////////////////////////////////////////////////////////////////////////////////
// Class:       MichelImageDrawerData
// Plugin Type: analyzer (art v2_11_03)
// File:        MichelImageDrawer_module.cc
//
// Generated at Tue Nov  6 14:17:47 2018 by Aidan Lewis Reynolds, using 
// cetskelgen from cetlib version v3_03_01.
//
// TODO: * Expand purity and efficiency tests 
// 
////////////////////////////////////////////////////////////////////////////////

// Art framework includes for dealing with art events
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

// LArSoft includes for dealing with reconstructed larsoft objects
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/ArtDataHelper/MVAReader.h"
#include "larreco/RecoAlg/ImagePatternAlgs/DataProvider/DataProviderAlg.h"

// Analysis helper tools for helpful analyser functions
#include "dune/Protodune/Analysis/ProtoDUNETrackUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNEShowerUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNETruthUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNEPFParticleUtils.h"

// ROOT includes
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"

// C++ includes 
#include <iostream>
#include <vector>

// Michel Analysis includes for my reusable Michel analysis functions
#include "MichelAnalysis.h"

namespace MichelAnalysis {

	class MichelImageDrawerData; 

	struct Config
	{

		using Name = fhicl::Name;

		// Debug config
		fhicl::Atom<bool> Debug { Name("Debug") };

		// Event selection tunable parameters
		fhicl::Atom<double> FidVolCut { Name("FidVolCut") };

		// Data product labels
		fhicl::Atom<std::string> CalorimetryTag { Name("CalorimetryTag") };
		fhicl::Atom<std::string> TrackTag       { Name("TrackTag") };
		fhicl::Atom<std::string> ShowerTag      { Name("ShowerTag") };
		fhicl::Atom<std::string> PFParticleTag  { Name("PFParticleTag") };
		fhicl::Atom<std::string> MCParticleTag  { Name("MCParticleTag") };
		fhicl::Atom<std::string> NNetTag        { Name("NNetTag") };

		fhicl::Atom<size_t> NWireBins        { Name("NWireBins") };

	};

}

using Parameters = art::EDAnalyzer::Table<MichelAnalysis::Config>;

class MichelAnalysis::MichelImageDrawerData : public art::EDAnalyzer 
{

public:

	explicit MichelImageDrawerData(const Parameters & config);

	// Plugins should not be copied or assigned.
	MichelImageDrawerData(MichelImageDrawerData const &) = delete;
	MichelImageDrawerData(MichelImageDrawerData &&) = delete;
	MichelImageDrawerData & operator = (MichelImageDrawerData const &) = delete;
	MichelImageDrawerData & operator = (MichelImageDrawerData &&) = delete;

	// Required functions.
	void analyze(art::Event const & event) override;

	// selected optional functions.
	void beginJob() override;
	void endJob() override;
	void reconfigure(fhicl::ParameterSet const & pset);

private:

	//////////////////////////////////////////////////////////////////////////////
	// Private variables
	//////////////////////////////////////////////////////////////////////////////

	// Debug and config options
	bool fDebug {false};

	// Image counting
	int nImagesDrawn {0};

	// Data product labels, read from fhicl
	std::string fCalorimetryTag;
	std::string fTrackTag;
	std::string fShowerTag;
	std::string fPFParticleTag;  
	std::string fMCParticleTag;  
	std::string fNNetTag;

	// Event selection tunable parameters
	double fFidVolCut;

	// Image data
	size_t fNWireBins;

	// Relevant data structs to be filled into TTrees
	MichelAnalysis::HitData            fHit;
	MichelAnalysis::DaughterData       fDaughter;
	MichelAnalysis::EventSelectionData fEventSelection;

	// Selected primary tracks
	std::vector<MichelAnalysis::TrackData> fSelectedTracks;

	// TTrees
	TTree * fHitTree            { nullptr };
	TTree * fDaughterTree       { nullptr };
	TTree * fEventSelectionTree { nullptr };

	// Canvas
	TCanvas * fCanvas = new TCanvas("canv", "canv");

	//////////////////////////////////////////////////////////////////////////////
	// Private Functions
	//////////////////////////////////////////////////////////////////////////////

	void InitialiseHitTree(art::ServiceHandle<art::TFileService> tfs);
	void InitialiseDaughterTree(art::ServiceHandle<art::TFileService> tfs);
	void InitialiseEventSelectionTree(art::ServiceHandle<art::TFileService> tfs);

};

// TODO: I could move back to default constructor now I got rid of the 
//       dataprovider
MichelAnalysis::MichelImageDrawerData::MichelImageDrawerData(Parameters const & config) 
	: EDAnalyzer(config), fDebug(config().Debug()), 
	  fCalorimetryTag(config().CalorimetryTag()), fTrackTag(config().TrackTag()), 
	  fShowerTag(config().ShowerTag()), fPFParticleTag(config().PFParticleTag()), 
	  fMCParticleTag(config().MCParticleTag()), fNNetTag(config().NNetTag()), 
	  fFidVolCut(config().FidVolCut()), fNWireBins(config().NWireBins())
{ 
}

void MichelAnalysis::MichelImageDrawerData::beginJob()
{

	std::cout << "------------------------------ Begin Michel Energy Image Job"
	          << " ------------------------------" << std::endl;

	art::ServiceHandle<art::TFileService> tfs;

	MichelAnalysis::MichelImageDrawerData::InitialiseHitTree(tfs);
	MichelAnalysis::MichelImageDrawerData::InitialiseDaughterTree(tfs);
	MichelAnalysis::MichelImageDrawerData::InitialiseEventSelectionTree(tfs);

}

void MichelAnalysis::MichelImageDrawerData::endJob()
{

	std::cout << "------------------------------ End Michel Energy Image Job"
	          << " ------------------------------" << std::endl;

}

// Defines the branch structure of the hit by hit tree
void MichelAnalysis::MichelImageDrawerData::InitialiseHitTree(
  art::ServiceHandle<art::TFileService> tfs)
{

	fHitTree = tfs -> make<TTree>("hit", "Hit info");

	fHitTree -> Branch("fHitPeakTime",        & fHit.PeakTime);
	fHitTree -> Branch("fHitIntegral",        & fHit.Integral);
	fHitTree -> Branch("fHitCNNScore_Michel", & fHit.CNNScore_Michel);
	fHitTree -> Branch("fHitCNNScore_EM",     & fHit.CNNScore_EM);

}

// Defines the branch structure of the daughter by daughter tree
void MichelAnalysis::MichelImageDrawerData::InitialiseDaughterTree(
  art::ServiceHandle<art::TFileService> tfs)
{

	fDaughterTree = tfs -> make<TTree>("daughter", "daughter info");

	fDaughterTree -> Branch("fDaughterNHits", 
	                        & fDaughter.NHits);
	fDaughterTree -> Branch("fDaughterNMichelHits", 
	                        & fDaughter.NMichelHits);
	fDaughterTree -> Branch("fDaughterFractionMichelHits", 
	                        & fDaughter.FractionMichelHits);

}

// Defines the branch structure of the event selection tree
void MichelAnalysis::MichelImageDrawerData::InitialiseEventSelectionTree(
  art::ServiceHandle<art::TFileService> tfs)
{

	fEventSelectionTree = tfs -> make<TTree>("eventselection", 
	                                         "event selection info");

	fEventSelectionTree -> Branch("fEventSelectionNTrue", 
	                              & fEventSelection.NTrue);
	fEventSelectionTree -> Branch("fEventSelectionNSelected", 
	                              & fEventSelection.NSelected);
	fEventSelectionTree -> Branch("fEventSelectionNCorrectlySelected", 
	                              & fEventSelection.NCorrectlySelected);
	fEventSelectionTree -> Branch("fEventSelectionNIncorrectlySelected", 
	                              & fEventSelection.NIncorrectlySelected);

}

// Main analyzer function
void MichelAnalysis::MichelImageDrawerData::analyze(art::Event const & event)
{

	art::ServiceHandle<art::TFileService> tfs;
	const size_t eventNumber  = event.id().event();
	const size_t runNumber    = event.id().run();
	std::cout << "Moving to event: " << eventNumber << " in run " 
	          << runNumber << std::endl;

	// reset event selection data
	fEventSelection = MichelAnalysis::ResetEventSelectionData();

	// Get Services
	std::cout << "-------------------------------- Getting Services Objects "
	          << "---------------------------------" <<std::endl;

	MichelAnalysis::ProtoDUNEUtils pdUtils = MichelAnalysis::GetProtoDUNEUtils();
	MichelAnalysis::DUNEUtils duneUtils = MichelAnalysis::GetDUNEUtils();

	// Get reconstructed objects
	std::cout << "------------------------------ Getting Reconstructed Objects "
	          << "------------------------------" <<std::endl;


	vh_vec_pfp_t pfparticles = event.getValidHandle<vec_pfp_t>(fPFParticleTag);
	mva_hit_t hitcnnscores(event, fNNetTag);

	// Event selection algorithm is finding Michels based on looking through the 
	// daughters of primary pfparticles
	for (auto const & pfparticle : *pfparticles) 
	{

		bool primarySelected { false };

		// Discard pfparticles that are not primary: we are looking for a muon at 
		// first
		if (!MichelAnalysis::IsPrimaryTrack(pfparticle, fPFParticleTag, fTrackTag, 
		                                    event, pdUtils)) 
		{
			continue;
		}

		const recob::Track * primarytrack = 
		  pdUtils.pfputil.GetPFParticleTrack(pfparticle, event, fPFParticleTag, 
		                                     fTrackTag);

		// Discard none track like primaries: muons should be track like
		if (!MichelAnalysis::IsGoodPrimaryTrack(primarytrack)) { continue; }

		// Now we have a muon we want to look through its daughters
		const std::vector<const recob::Shower *> daughterShowers =
		  pdUtils.pfputil.GetPFParticleDaughterShowers(pfparticle, event, 
		                                               fPFParticleTag, fShowerTag);

		for (auto const & daughterShower : daughterShowers) 
		{

			// Skipping if primary already selected
			// TODO: define a "best daughter" as the one with the highest michel hit 
			//       fraction 
			if (primarySelected) { break; }

			// The daughter is analysed based on the CNN scores associated to its hits
			//   - The fraction of hits with a "high" michel like CNN score is used to
			//     select events
			const vec_ptr_hit_t daughterShowerHits = 
			  MichelAnalysis::GetRecoShowerHits(* daughterShower, event, fShowerTag);

			MichelAnalysis::AnalyseDaughterShowerHits(fDaughter, fHit, duneUtils,
			                                          daughterShowerHits, 
			                                          hitcnnscores, fHitTree);

			fDaughterTree -> Fill();

			// Main event selection
			//   - Fraction michel hits
			//   - Inside fiducial volume
			if (!MichelAnalysis::EventSelection(*primarytrack, *daughterShower,
			                                    fDaughter)) 
			{ 
				continue; 
			}

			primarySelected = true;

			// Add the primary to the list of selected tracks for filling into the 
			// track tree
			fSelectedTracks.push_back(
			  MichelAnalysis::GetTrackData(primarytrack, pdUtils, event, fTrackTag));

			// Image drawing: name indicates if the event was correctly or incorrectly
			// selected
			std::ostringstream buffer;
			buffer << "event_" << event.id().event() << "primary_" << 
			            pfparticle.Self();
			std::string basename = buffer.str();

			TH2D patchHist = MichelAnalysis::DrawMichelAsRawImage(daughterShowerHits, 
			                                                      event, basename, 
			                                                      tfs, fNWireBins,
			                                                      duneUtils);

			if (patchHist.GetEntries() == 0) { continue; }

			nImagesDrawn += 1;

			std::stringstream histTitle;
			histTitle << "Raw Data Patch Around Selected Michel Electron. "
			          << "Event: " << eventNumber;
			patchHist.SetTitle(histTitle.str().c_str());

			patchHist.Draw("colz");
			fCanvas -> SaveAs(Form("data_patch_%i.png", nImagesDrawn));
			fCanvas -> Clear();

		}

	}

	std::cout << "Selected michel like particles in event: " 
	          << fEventSelection.NSelected << std::endl;

	fEventSelectionTree -> Fill();

}

DEFINE_ART_MODULE(MichelAnalysis::MichelImageDrawerData)
