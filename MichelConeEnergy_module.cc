////////////////////////////////////////////////////////////////////////////////
// Class:       MichelConeEnergy
// Plugin Type: analyzer (art v2_11_03)
// File:        MichelConeEnergy_module.cc
//
// Generated at Thu Nov 22 11:49:28 2018 by Aiden Reynolds using cetskelgen
// from cetlib version v3_03_01.
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
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "lardata/ArtDataHelper/MVAReader.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"

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

namespace MichelAnalysis 
{
	class MichelConeEnergy;
}

class MichelAnalysis::MichelConeEnergy : public art::EDAnalyzer 
{

public:
	explicit MichelConeEnergy(fhicl::ParameterSet const & p);

	// Plugins should not be copied or assigned.
	MichelConeEnergy(MichelConeEnergy const &) = delete;
	MichelConeEnergy(MichelConeEnergy &&) = delete;
	MichelConeEnergy & operator = (MichelConeEnergy const &) = delete;
	MichelConeEnergy & operator = (MichelConeEnergy &&) = delete;

	// Required functions.
	void analyze(art::Event const & event) override;

	// selected optional functions.
	void beginJob() override;
	void reconfigure(fhicl::ParameterSet const & p); 

private:

	//////////////////////////////////////////////////////////////////////////////
	// Private variables and objects
	//////////////////////////////////////////////////////////////////////////////

	// Data product labels, read from fhicl
	std::string fTrackTag;
	std::string fShowerTag;
	std::string fPFParticleTag;  
	std::string fMCParticleTag;  
	std::string fNNetTag;

	// Event selection tunable parameters
	double fFidVolCut;

	calo::CalorimetryAlg fCaloAlg;

	// Relevant data structs to be filled into TTrees
	// We have to do it this way so that the tree has the address of the relevant
	// data
	MichelAnalysis::HitData            fHit;
	MichelAnalysis::DaughterData       fDaughter;
	MichelAnalysis::TrackData          fPrimary;
	MichelAnalysis::MichelEvent        fMichelEvent;

	// TTrees
	TTree * fHitTree            { nullptr };
	TTree * fDaughterTree       { nullptr };
	TTree * fMichelEventTree    { nullptr };

	//////////////////////////////////////////////////////////////////////////////
	// Private Functions
	//////////////////////////////////////////////////////////////////////////////

	void InitialiseHitTree(
	       const art::ServiceHandle<art::TFileService> & tfs);
	void InitialiseDaughterTree(
	       const art::ServiceHandle<art::TFileService> & tfs);
	void InitialiseMichelEventTree(
	       const art::ServiceHandle<art::TFileService> & tfs);

};


MichelAnalysis::MichelConeEnergy::MichelConeEnergy(
	fhicl::ParameterSet const & p)
	:
	EDAnalyzer(p), 
	fCaloAlg (p.get<fhicl::ParameterSet>("CalorimetryTag"))
{
	reconfigure(p);
}

void MichelAnalysis::MichelConeEnergy::reconfigure(
       fhicl::ParameterSet const & p) 
{

	fTrackTag       = p.get<std::string>("TrackTag");
	fShowerTag      = p.get<std::string>("ShowerTag");
	fPFParticleTag  = p.get<std::string>("PFParticleTag");
	fNNetTag        = p.get<std::string>("NNetTag");
	fFidVolCut      = p.get<double>("FidVolCut");

	return;

}

void MichelAnalysis::MichelConeEnergy::beginJob()
{

	std::cout << "------------------------------ Begin Michel Energy Image Job"
	          << " ------------------------------" << std::endl;

	art::ServiceHandle<art::TFileService> tfs;
	MichelAnalysis::MichelConeEnergy::InitialiseHitTree(tfs);
	MichelAnalysis::MichelConeEnergy::InitialiseDaughterTree(tfs);
	MichelAnalysis::MichelConeEnergy::InitialiseMichelEventTree(tfs);

}

// Defines the branch structure of the hit by hit tree
void MichelAnalysis::MichelConeEnergy::InitialiseHitTree(
       const art::ServiceHandle<art::TFileService> & tfs)
{

	fHitTree = tfs->make<TTree>("hit", "Hit info");

	fHitTree->Branch("PeakTime",        & fHit.PeakTime);
	fHitTree->Branch("Integral",        & fHit.Integral);
	fHitTree->Branch("CNNScore_Michel", & fHit.CNNScore_Michel);
	fHitTree->Branch("CNNScore_EM",     & fHit.CNNScore_EM);

}

// Defines the branch structure of the daughter by daughter tree
void MichelAnalysis::MichelConeEnergy::InitialiseDaughterTree(
       const art::ServiceHandle<art::TFileService> & tfs)
{

	fDaughterTree = tfs->make<TTree>("daughter", "daughter info");

	fDaughterTree->Branch(
	 "NHits",               & fDaughter.NHits);
	fDaughterTree->Branch(
	  "NMichelHits",        & fDaughter.NMichelHits);
	fDaughterTree->Branch(
	  "FractionMichelHits", & fDaughter.FractionMichelHits);

	fDaughterTree->Branch(
	  "MichelScores",       & fDaughter.MichelScores);
	fDaughterTree->Branch(
	  "EMScores",           & fDaughter.EMScores);
	fDaughterTree->Branch(
	  "Integrals",          & fDaughter.Integrals);
	fDaughterTree->Branch(
	  "PeakTimes",          & fDaughter.PeakTimes);

	fDaughterTree->Branch(
	  "DistanceToPrimary",  & fDaughter.DistanceToPrimary);
	fDaughterTree->Branch(
	  "DistanceToPrimaryX", & fDaughter.DistanceToPrimaryX);
	fDaughterTree->Branch(
	  "DistanceToPrimaryY", & fDaughter.DistanceToPrimaryY);
	fDaughterTree->Branch(
	  "DistanceToPrimaryZ", & fDaughter.DistanceToPrimaryZ);

}

// defines the branch structure of the track by track tree
void MichelAnalysis::MichelConeEnergy::InitialiseMichelEventTree(
       const art::ServiceHandle<art::TFileService> & tfs)
{

	fMichelEventTree = tfs->make<TTree>("michelevent", "Michel event info");

	fMichelEventTree->Branch(
	  "TrackID",                    & fMichelEvent.Primary.ID);

	fMichelEventTree->Branch(
	  "TrackIsMuonDecaying",        & fMichelEvent.Primary.IsMuonDecaying);

	fMichelEventTree->Branch(
	  "TrackVertexX",               & fMichelEvent.Primary.VertexX);
	fMichelEventTree->Branch(
	  "TrackVertexY",               & fMichelEvent.Primary.VertexY);
	fMichelEventTree->Branch(
	  "TrackVertexZ",               & fMichelEvent.Primary.VertexZ);

	fMichelEventTree->Branch(
	  "TrackEndX",                  & fMichelEvent.Primary.EndX);
	fMichelEventTree->Branch(
	  "TrackEndY",                  & fMichelEvent.Primary.EndY);
	fMichelEventTree->Branch(
	  "TrackEndZ",                  & fMichelEvent.Primary.EndZ);

	fMichelEventTree->Branch(
	  "TrackHasT0",                 & fMichelEvent.Primary.HasT0);
	fMichelEventTree->Branch(
	  "TrackT0",                    & fMichelEvent.Primary.T0);

	fMichelEventTree->Branch(
	  "DaughterNHits",              & fMichelEvent.Daughter.NHits);
	fMichelEventTree->Branch(
	  "DaughterNMichelHits",        & fMichelEvent.Daughter.NMichelHits);
	fMichelEventTree->Branch(
	  "DaughterFractionMichelHits", & fMichelEvent.Daughter.FractionMichelHits);

	fMichelEventTree->Branch(
	  "DaughterMichelScores",       & fMichelEvent.Daughter.MichelScores);
	fMichelEventTree->Branch(
	  "DaughterEMScores",           & fMichelEvent.Daughter.EMScores);
	fMichelEventTree->Branch(
	  "DaughterIntegrals",          & fMichelEvent.Daughter.Integrals);
	fMichelEventTree->Branch(
	  "DaughterPeakTimes",          & fMichelEvent.Daughter.PeakTimes);

	fMichelEventTree->Branch(
	  "DaughterDistanceToPrimary",  & fMichelEvent.Daughter.DistanceToPrimary);
	fMichelEventTree->Branch(
	  "DaughterDistanceToPrimaryX", & fMichelEvent.Daughter.DistanceToPrimaryX);
	fMichelEventTree->Branch(
	  "DaughterDistanceToPrimaryY", & fMichelEvent.Daughter.DistanceToPrimaryY);
	fMichelEventTree->Branch(
	  "DaughterDistanceToPrimaryZ", & fMichelEvent.Daughter.DistanceToPrimaryZ);

}

void MichelAnalysis::MichelConeEnergy::analyze(art::Event const & event)
{

	std::cout << "Start of analyzer" << std::endl;

	// Get all the data products
	// Utilites
	MichelAnalysis::ProtoDUNEUtils pdUtils = MichelAnalysis::GetProtoDUNEUtils();
	MichelAnalysis::DUNEUtils duneUtils = MichelAnalysis::GetDUNEUtils();
	// Reco data products
	vh_vec_pfp_t pfparticles = event.getValidHandle<vec_pfp_t>(fPFParticleTag);
	mva_hit_t hitcnnscores(event, fNNetTag);

	// Loop over all primary pfparticles to find muons
	for (auto const & pfparticle : * pfparticles) 
	{

		// Cut to find primary tracks
		if (!MichelAnalysis::IsPrimaryTrack(pfparticle, fPFParticleTag, fTrackTag, 
		                                    event, pdUtils)) 
		{
			continue;
		}

		const recob::Track * primarytrack = 
		  pdUtils.pfputil.GetPFParticleTrack(pfparticle, event, fPFParticleTag, 
		                                     fTrackTag);
		if (!MichelAnalysis::IsGoodPrimaryTrack(primarytrack)) { continue; }

		// Get details of this primary
		fPrimary = MichelAnalysis::GetTrackData(primarytrack, pdUtils, event, 
		                                        fTrackTag );

		// Look for all daughter showers of this track to get michels
		const std::vector<const recob::Shower *> daughterShowers =
		  pdUtils.pfputil.GetPFParticleDaughterShowers(pfparticle, event, 
		                                               fPFParticleTag, fShowerTag);

		// Keep track of the most Michel like daughter of this track
		const recob::Shower * bestShower {nullptr}; 
		MichelAnalysis::DaughterData bestDaughter;
		float bestMichelFraction = 0.; 

		// For all daughters get hold of the Michel-like score
		for (auto const & daughterShower : daughterShowers) 
		{

			if (daughterShower == nullptr) { continue; }

			// Resest daughter data so we don't get old and new data mixed up
			fDaughter = MichelAnalysis::ResetDaughterData();

			// Get hold of the hits in the shower and analyse them to look for michel 
			// like hits
			const vec_ptr_hit_t daughterShowerHits = 
			  MichelAnalysis::GetRecoShowerHits(* daughterShower, event, fShowerTag);

			MichelAnalysis::AnalyseDaughterShowerHits(fDaughter, fHit, duneUtils,
			                                          daughterShowerHits, 
			                                          hitcnnscores, fPrimary.T0,
			                                          fCaloAlg, event, fHitTree);

			// Update best daughter if applicable
			if (fDaughter.FractionMichelHits > bestMichelFraction) 
			{
				bestMichelFraction = fDaughter.FractionMichelHits;
				bestShower = daughterShower;
				bestDaughter = fDaughter;
			}

			fDaughterTree->Fill();

		}

		if (bestShower == nullptr) { continue; }

		// Now take our best daughter and do the event selection based on the hits
		const vec_ptr_hit_t daughterShowerHits = 
		  MichelAnalysis::GetRecoShowerHits(* bestShower, event, fShowerTag);

		// Make a michel event object that has details of the daughter and the 
		// primary.
		fMichelEvent = { fPrimary, bestDaughter };
		fMichelEventTree->Fill();

	}

	std::cout << "End of analyzer" << std::endl;

}

DEFINE_ART_MODULE(MichelAnalysis::MichelConeEnergy)
