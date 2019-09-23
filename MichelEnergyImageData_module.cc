// Plugin Type: analyzer (art v2_11_03)
// File:        MichelEnergyImage_module.cc
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
	class MichelEnergyImage;
}

class MichelAnalysis::MichelEnergyImage : public art::EDAnalyzer 
{

public:
	explicit MichelEnergyImage(fhicl::ParameterSet const & p);

	// Plugins should not be copied or assigned.
	MichelEnergyImage(MichelEnergyImage const &) = delete;
	MichelEnergyImage(MichelEnergyImage &&) = delete;
	MichelEnergyImage & operator = (MichelEnergyImage const &) = delete;
	MichelEnergyImage & operator = (MichelEnergyImage &&) = delete;

	// Required functions.
	void analyze(art::Event const & event) override;

	// selected optional functions.
	void reconfigure(fhicl::ParameterSet const & p); 

private:

	//////////////////////////////////////////////////////////////////////////////
	// Private variables and objects
	//////////////////////////////////////////////////////////////////////////////

	// Data product labels, read from fhicl
	std::string fTrackTag;
	std::string fShowerTag;
	std::string fPFParticleTag;  
	std::string fNNetTag;
	std::string fCaloTag;
	std::string fHitTag;

	// Image size etc
	size_t fNWireBins;

	// Drawing and name bookeping
	TCanvas * fCanvas = new TCanvas("canv", "canv");
	int nImagesDrawn {0};

	// Event selection tunable parameters
	double fFidVolCut;

	calo::CalorimetryAlg fCaloAlg;

	// Relevant data structs to be filled into TTrees
	MichelAnalysis::HitData            fHit;
	MichelAnalysis::DaughterData       fDaughter;
	MichelAnalysis::TrackData          fPrimary;

};


MichelAnalysis::MichelEnergyImage::MichelEnergyImage(
	fhicl::ParameterSet const & p)
	:
	EDAnalyzer(p), 
	fCaloAlg (p.get<fhicl::ParameterSet>("CalorimetryTag"))
{
	reconfigure(p);
}

void MichelAnalysis::MichelEnergyImage::reconfigure(
       fhicl::ParameterSet const & p) 
{

	fTrackTag      = p.get<std::string>("TrackTag");
	fShowerTag     = p.get<std::string>("ShowerTag");
	fPFParticleTag = p.get<std::string>("PFParticleTag");
	fNNetTag       = p.get<std::string>("NNetTag");
	fCaloTag       = p.get<std::string>("CaloTag");
	fHitTag        = p.get<std::string>("HitTag");
	fFidVolCut     = p.get<double>("FidVolCut");
	fNWireBins     = p.get<size_t>("NWireBins");

	return;

}

void MichelAnalysis::MichelEnergyImage::analyze(art::Event const & event)
{

	art::ServiceHandle<art::TFileService> tfs;

	//////////////////////////////////////////////////////////////////////////////
	// Get all the data products

	// Utilites
	MichelAnalysis::ProtoDUNEUtils pdUtils = MichelAnalysis::GetProtoDUNEUtils();
	MichelAnalysis::DUNEUtils duneUtils    = MichelAnalysis::GetDUNEUtils();

	// Reco data products
	vh_vec_pfp_t pfparticles = event.getValidHandle<vec_pfp_t>(fPFParticleTag);
	mva_hit_t hitcnnscores(event, fNNetTag);

	art::Handle<std::vector<recob::Hit> > hitCollector;
	vec_ptr_hit_t allHits;
	if (event.getByLabel(fHitTag, hitCollector)) 
	{
		art::fill_ptr_vector(allHits, hitCollector);
	}

	auto allClusters = 
	  event.getValidHandle<std::vector<recob::Cluster>>(fPFParticleTag);
	//////////////////////////////////////////////////////////////////////////////

	// Loop over all primary pfparticles to find muons
	for (auto const & pfparticle : * pfparticles) 
	{

		// Cut to find primary tracks
		if (!MichelAnalysis::IsPrimaryTrack(pfparticle, fPFParticleTag, fTrackTag, 
		                                    event, pdUtils)) 
		{ 
			continue; 
		};

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
		const recob::Shower * bestShower { nullptr }; 
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
			                                          fCaloAlg);

			// Update best daughter if applicable
			if (fDaughter.FractionMichelHits > bestMichelFraction) 
			{
				bestMichelFraction = fDaughter.FractionMichelHits;
				bestShower         = daughterShower;
				bestDaughter       = fDaughter;
			}

		}

		if (bestShower == nullptr) { continue; }

		const vec_ptr_hit_t showerHits = 
		  MichelAnalysis::GetRecoShowerHits(* bestShower, event, fShowerTag);

		// Cut to ensure enough collection hits in shower
		size_t nColHits { 0 };
		for (auto const & hit : showerHits)
		{
			if (hit->WireID().Plane != 2) { continue; }
			nColHits += 1;
		}
		if (nColHits < 5) { continue; }

		if (!MichelAnalysis::EventSelection(* primarytrack, * bestShower, 
		                                    bestDaughter)) 
		{ continue; }

		// Calculate normalisation factor
		// Using calos make the 60cm 10 bin averaged dQ/dx constant
		size_t nHitsUsed { 0 };
		float avgDQDX { 0. };
		std::vector<anab::Calorimetry> calos = 
		  pdUtils.trackutil.GetRecoTrackCalorimetry(* primarytrack, event, 
		                                              fTrackTag, fCaloTag);

		for (size_t itcal = 0; itcal < calos.size(); itcal++) 
		{

			size_t const Nhits = calos[itcal].dEdx().size();

			for (size_t itHit = 0; itHit < Nhits; itHit++)  
			{

				if (nHitsUsed >= 10) { break; }
				auto dqdx = (calos[itcal].dQdx())[itHit];
				auto resrange = (calos[itcal].ResidualRange())[itHit];

				if (resrange > 60.f)
				{
					avgDQDX += dqdx;
					nHitsUsed += 1;
				}

			}

		}

		if (nHitsUsed < 10) { continue; }

		// Normalise to 300 avg dqdx ---> don't need to be too accurate but this is
		// a reasonable value
		avgDQDX /= static_cast<float>(nHitsUsed);
		float norm { 300.f / avgDQDX };

		// If norm is too significant ignore event
		if (norm < 0.33f || norm > 3.f) { continue; }

		std::cout << "Making image..." << std::flush;

		std::ostringstream buffer;
		buffer << "event_" << event.id().event() << "_primary" << pfparticle.Self();
		std::string basename = buffer.str();

		std::vector<TH2D> trainingData = 
		  DrawMichelTrainingData_Data(showerHits, allHits, * pfparticles, 
		                              fPFParticleTag, hitcnnscores, fPrimary.T0, 
		                              fCaloAlg, event, basename, tfs, fNWireBins, 
		                              norm, pdUtils, duneUtils);

		if (trainingData.size() == 0) { continue; } 
		if (trainingData[0].GetEntries() == 0) { continue; }

		TH2D rawImage = DrawMichelAsRawImage(showerHits, event, basename, tfs, 
		                                     fNWireBins, duneUtils);

		std::cout << "\n";

	}

}

DEFINE_ART_MODULE(MichelAnalysis::MichelEnergyImage)
