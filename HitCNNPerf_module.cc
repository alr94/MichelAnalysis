////////////////////////////////////////////////////////////////////////
// Class:       HitCNNPerf
// Plugin Type: analyzer (art v2_07_03)
// File:        HitCNNPerf_module.cc
//
// Generated at Mon Sep  4 06:55:33 2017 by Leigh Whitehead using cetskelgen
// from cetlib version v3_00_01.
//
// This module is designed to show some usage examples of the analysis
// tools that I have been producing for protoDUNE. The aim has been to
// simplify the associations between the different objects to make
// some of the low-level art features more transparent
//
// The code is split into a few sections with each focusing on a different
// type of initial object
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/ArtDataHelper/MVAWrapperBase.h"
#include "lardata/ArtDataHelper/MVAReader.h"

#include "dune/Protodune/Analysis/ProtoDUNETrackUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNEShowerUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNETruthUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNEPFParticleUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNEDataUtils.h"

#include "TTree.h"

namespace nnet {
	class HitCNNPerf;
}


class nnet::HitCNNPerf : public art::EDAnalyzer {
public:

	explicit HitCNNPerf(fhicl::ParameterSet const & p);
	// The compiler-generated destructor is fine for non-base
	// classes without bare pointers or other resource use.

	// Plugins should not be copied or assigned.
	HitCNNPerf(HitCNNPerf const &) = delete;
	HitCNNPerf(HitCNNPerf &&) = delete;
	HitCNNPerf & operator = (HitCNNPerf const &) = delete;
	HitCNNPerf & operator = (HitCNNPerf &&) = delete;

	virtual void beginJob() override;
	virtual void endJob() override;

	// Required functions.
	void analyze(art::Event const & e) override;

private:

	TTree* fHitTree;

	bool fIsTrackHit;
	double fTrackScore;
	double fShowerScore;
	int fEvent;
};


nnet::HitCNNPerf::HitCNNPerf(fhicl::ParameterSet const & p)
	:
	EDAnalyzer(p)
{

}

void nnet::HitCNNPerf::beginJob()
{

	// Make an output tre to store information on pandora hits
	art::ServiceHandle<art::TFileService> tfs;

	fHitTree = tfs->make<TTree>("HitTree", "Hit Tree");
	fHitTree->Branch("isTrackHit",&fIsTrackHit,"isTrackHit/O");
	fHitTree->Branch("trackScore",&fTrackScore,"trackScore/D");
	fHitTree->Branch("showerScore",&fShowerScore,"showerScore/D");
	fHitTree->Branch("event",&fEvent,"event/I");
}

void nnet::HitCNNPerf::analyze(art::Event const & evt)
{

	fEvent = evt.id().event();

	// ========
	// Get results and products from the CNN

	auto cluResults = 
	  anab::MVAReader< recob::Cluster, 4 >::create(evt, 
	                                               "emtrkmichelid:emtrkmichel");
 
	if (!cluResults) { return; }
	
	int trkLikeIdx = cluResults->getIndex("track");
	int emLikeIdx = cluResults->getIndex("em");
	if ((trkLikeIdx < 0) || (emLikeIdx < 0)) { return; }
	
	const art::FindManyP< recob::Hit > hitsFromClusters(cluResults->dataHandle(), 
	                                                    evt, 
	                                                    cluResults->dataTag());
	const auto & cnnOuts = cluResults->outputs();

	// ========
	// Get the pandora showers

	auto recoShowers = 
	  evt.getValidHandle<std::vector<recob::Shower>>("pandoraShower");
	const art::FindManyP<recob::Hit> hitsFromShowers(recoShowers, evt,
	                                                 "pandoraShower");

	for(unsigned int s = 0; s < recoShowers->size(); ++s ){

		std::vector<art::Ptr<recob::Hit>> hits = hitsFromShowers.at(s);

		for(auto const h : hits){

			double shwscore = -999;
			double trkscore = -999;
			fIsTrackHit = false;

			for(unsigned int c = 0; c < cnnOuts.size(); ++c){

				double trkOrEm = cnnOuts[c][trkLikeIdx] + cnnOuts[c][emLikeIdx];
				if (trkOrEm > 0) { 
					shwscore = cnnOuts[c][emLikeIdx] / trkOrEm;
					trkscore = cnnOuts[c][trkLikeIdx] / trkOrEm;
				}
 
				for(auto const ch : hitsFromClusters.at(c)){
					if(h.key() == ch.key()){
						fShowerScore = shwscore;
						fTrackScore  = trkscore;
					}
				}
			}
			fHitTree->Fill();
		}
	}

	// ========
	// Get the pandora tracks

	auto recoTracks = 
		evt.getValidHandle<std::vector<recob::Track>>("pandoraTrack");
	const art::FindManyP<recob::Hit> hitsFromTracks(recoTracks, evt,
	                                                "pandoraTrack");

	for(unsigned int t = 0; t < recoTracks->size(); ++t ){

		std::vector<art::Ptr<recob::Hit>> hits = hitsFromTracks.at(t);

		for(auto const h : hits){

			double shwscore = -999;
			double trkscore = -999;
			fIsTrackHit = true;

			for(unsigned int c = 0; c < cnnOuts.size(); ++c){

				double trkOrEm = cnnOuts[c][trkLikeIdx] + cnnOuts[c][emLikeIdx];
				if (trkOrEm > 0) { 
					shwscore = cnnOuts[c][emLikeIdx] / trkOrEm;
					trkscore = cnnOuts[c][trkLikeIdx] / trkOrEm;
				}

				for(auto const ch : hitsFromClusters.at(c)){
					if(h.key() == ch.key()){
						fShowerScore = shwscore;
						fTrackScore  = trkscore;
					}
				}
			}
			fHitTree->Fill();
		}
	}

}

void nnet::HitCNNPerf::endJob()
{

}

DEFINE_ART_MODULE(nnet::HitCNNPerf)

