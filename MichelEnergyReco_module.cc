// Plugin Type: analyzer (art v2_11_03)
// File:        MichelEnergyReco_module.cc
//
// Generated at Wed June 19 11:49:28 2019 by Aiden Reynolds using cetskelgen
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
	class MichelEnergyReco;
}

class MichelAnalysis::MichelEnergyReco : public art::EDAnalyzer
{

////////////////////////////////////////////////////////////////////////////////
public:
	explicit MichelEnergyReco(fhicl::ParameterSet const & p);

	// Plugins should not be copied or assigned
	MichelEnergyReco(MichelEnergyReco const &) = delete;
	MichelEnergyReco(MichelEnergyReco &&) = delete;
	MichelEnergyReco & operator = (MichelEnergyReco const &) = delete;
	MichelEnergyReco & operator = (MichelEnergyReco &&) = delete;

	// Required functions
	void analyze(art::Event const & event) override;

	// Selected optional functions
	void reconfigure(fhicl::ParameterSet const &p);

////////////////////////////////////////////////////////////////////////////////
private:

	// Data product tags
	std::string fHitTag;
	std::string fTrackTag;
	std::string fShowerTag;
	std::string fPFParticleTag;  
	std::string fMCParticleTag;  
	std::string fNNetTag;

	// Services
	calo::CalorimetryAlg fCaloAlg;

	// Event selection tunable parameters
	double fFidVolCut;

	// Relevant personal data structures
	MichelAnalysis::HitData            fHit;
	MichelAnalysis::DaughterData       fDaughter;
	MichelAnalysis::TrackData          fPrimary;

};


MichelAnalysis::MichelEnergyReco::MichelEnergyReco(
	fhicl::ParameterSet const & p)
	:
	EDAnalyzer(p), 
	fCaloAlg (p.get<fhicl::ParameterSet>("CalorimetryTag"))
{
	reconfigure(p);
}

void MichelAnalysis::MichelEnergyReco::reconfigure(
       fhicl::ParameterSet const & p) 
{

	fHitTag        = p.get<std::string>("HitTag");
	fTrackTag      = p.get<std::string>("TrackTag");
	fShowerTag     = p.get<std::string>("ShowerTag");
	fPFParticleTag = p.get<std::string>("PFParticleTag");
	fMCParticleTag = p.get<std::string>("MCParticleTag");
	fNNetTag       = p.get<std::string>("NNetTag");
	fFidVolCut     = p.get<double>("FidVolCut");

	return;

}

void MichelAnalysis::MichelEnergyReco::analyze(art::Event const & event)
{

	art::ServiceHandle<art::TFileService> tfs;

	//////////////////////////////////////////////////////////////////////////////
	// Read relevant data from file

	// Utils
	MichelAnalysis::ProtoDUNEUtils pdUtils = MichelAnalysis::GetProtoDUNEUtils();
	MichelAnalysis::DUNEUtils duneUtils = MichelAnalysis::GetDUNEUtils();

	// Reco
	vh_vec_pfp_t pfparticles = event.getValidHandle<vec_pfp_t>(fPFParticleTag);
	mva_hit_t hitcnnscores(event, fNNetTag);

	// art::Handle<std::vector<recob::Hit> > hitCollector;
	// vec_ptr_hit_t allHits;
	// if (event.getByLabel(fHitTag, hitCollector)) 
	// {
	// 	art::fill_ptr_vector(allHits, hitCollector);
	// }

	// MC
	// vh_vec_mcp_t mcparticles = event.getValidHandle<vec_mcp_t>(fMCParticleTag);
	// pid_map_t    particleMap = MichelAnalysis::BuildParticleMap(mcparticles);

	//////////////////////////////////////////////////////////////////////////////

	// Use primary pandora tracks as start point for Michel electron search
	for (auto const & pfparticle : * pfparticles) 
	{
		bool isPT = MichelAnalysis::IsPrimaryTrack(pfparticle, fPFParticleTag, 
		                                           fTrackTag, event, pdUtils);
		if (!isPT) { continue; }

		const recob::Track * primarytrack = 
		  pdUtils.pfputil.GetPFParticleTrack(pfparticle, event, fPFParticleTag, 
		                                     fTrackTag);

		if (!MichelAnalysis::IsGoodPrimaryTrack(primarytrack)) { continue; }

		std::cout << "Found good track\n" << std::endl;

		fPrimary = MichelAnalysis::GetTrackData(primarytrack, pdUtils, event,
		                                        fTrackTag);

		const std::vector<const recob::Shower *> daughterShowers =
		  pdUtils.pfputil.GetPFParticleDaughterShowers(pfparticle, event, 
		                                               fPFParticleTag, fShowerTag);

		const recob::Shower * bestShower { nullptr } ;
		MichelAnalysis::DaughterData bestDaughter;
		float bestMichelFraction = 0.;

		for (auto const & daughterShower : daughterShowers) 
		{

			if (daughterShower == nullptr) { continue; }

			fDaughter = MichelAnalysis::ResetDaughterData();

			const vec_ptr_hit_t daughterHits = 
			  MichelAnalysis::GetRecoShowerHits(* daughterShower, event, fShowerTag);

			MichelAnalysis::AnalyseDaughterShowerHits(fDaughter, fHit, duneUtils,
			                                          daughterHits, hitcnnscores,
			                                          fPrimary.T0, fCaloAlg);

			if (fDaughter.FractionMichelHits > bestMichelFraction)
			{
				bestMichelFraction = fDaughter.FractionMichelHits;
				bestShower         = daughterShower;
				bestDaughter       = fDaughter;
			}

		}

		if (bestShower == nullptr) { continue; }

		bool eventSelected = MichelAnalysis::EventSelection(* primarytrack, 
		                                                    * bestShower,
		                                                    bestDaughter);

		if (!eventSelected) { continue; }

		std::cout << "Event selected, reconstructing hits" << std::endl;

		const vec_ptr_hit_t showerHits = MichelAnalysis::GetRecoShowerHits(
		                                   * bestShower, event, fShowerTag);

		// Get details of one of the hits and use it as the centre of the image
		// try to use plane 2 if possible since that is the collection view
		unsigned int centralWire     { 0 };
		unsigned int centralTick     { 0 };
		unsigned int nCollectionHits { 0 };

		art::Ptr<recob::Hit>  firstHit;
		for (auto const & showerHit : showerHits)
		{

			if (showerHit->WireID().Plane != 2) { continue; }

			if (firstHit.isNull()) { firstHit = showerHit; }

			nCollectionHits += 1;
			centralWire     += showerHit->WireID().Wire;
			centralTick     += showerHit->PeakTime();

		}

		if (nCollectionHits < 1) { continue; }

		centralWire /= nCollectionHits;
		centralTick /= nCollectionHits;

		// Below commented olut for now due to unused vars
		// const unsigned int plane    = firstHit->WireID().Plane;
		// const unsigned int tpc      = firstHit->WireID().TPC;
		// const unsigned int cryostat = firstHit->WireID().Cryostat;

		// double driftVelocity   = duneUtils.detprop->DriftVelocity();
		// double clockFrequency  = duneUtils.detclock->TPCClock().Frequency();
		// 
		// double wireBinSize   = duneUtils.geom->WirePitch(plane);
		// double timeBinSize   = driftVelocity / clockFrequency;
		// 
		// double aspectRatio = wireBinSize / timeBinSize;

		// constexpr size_t nWireBins = 164;
		// const size_t nTimeBins = static_cast<size_t>(
		//   static_cast<double>(nWireBins) * aspectRatio);

		// const vec_ptr_hit_t selectedHits;
		// for (auto const & hit : allHits) 
		// {

		// 	if (hit->WireID().Plane != plane)       { continue; }
		// 	if (hit->WireID().TPC != tpc)           { continue; }
		// 	if (hit->WireID().Cryostat != cryostat) { continue; }

		// 	const unsigned int wire = hit->WireID().Wire;
		// 	const int wireDiff = static_cast<int> (centralWire) - 
		// 	                     static_cast<int> (wire);
		// 	if (abs(wireDiff) > static_cast<int>(nWireBins) / 2) { continue; }

		// 	const unsigned int time = hit->PeakTime();
		// 	const int timeDiff = static_cast<int> (centralTick) - 
		// 	                     static_cast<int> (time);
		// 	if (abs(timeDiff) > static_cast<int>(nTimeBins) / 2) { continue; }

		// 	const HitData & hitData = GetHitData(hit, hitcnnscores);
		// 	const bool hitIsMichel  = HitIsMichel(hit, duneUtils);

		// }

		// TODO: Select all extra hits
		// TODO: Purity and Completeness measures

	}

}
