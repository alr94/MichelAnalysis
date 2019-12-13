#ifndef MICHELANALYSIS_H
#define MICHELANALYSIS_H

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

// larsoft includes for dealing with reconstructed larsoft objects
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "lardata/ArtDataHelper/MVAReader.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Simulation/LArG4Parameters.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"

// Analysis helper tools for helpful analyser functions
#include "dune/Protodune/Analysis/ProtoDUNETrackUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNEShowerUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNETruthUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNEPFParticleUtils.h"

// ROOT includes
#include "TTree.h"
#include "TH2D.h"

// C++ includes
#include <iostream>
#include <vector>
#include <stack>

// Typedefs for readability
using vh_vec_pfp_t   = art::ValidHandle<std::vector<recob::PFParticle>>;
using vec_pfp_t      = std::vector<recob::PFParticle>;
using pfp_t          = recob::PFParticle;

using vh_vec_mcp_t   = art::ValidHandle<std::vector<simb::MCParticle>>;
using vec_mcp_t      = std::vector<simb::MCParticle>;
using mcp_t          = simb::MCParticle;

using vec_ptr_sp_t  = std::vector<art::Ptr<recob::SpacePoint>>;

using vec_ptr_hit_t = std::vector<art::Ptr<recob::Hit>>;
using mva_hit_t     = anab::MVAReader<recob::Hit, 4>;

using vh_vec_clu_t  = art::ValidHandle<std::vector<recob::Cluster>>;
using vec_ptr_clu_t = std::vector<art::Ptr<recob::Cluster>>;
using mva_clu_t     = anab::MVAReader<recob::Cluster, 4>;

using pid_map_t      = std::unordered_map<int, const simb::MCParticle *>; 

namespace MichelAnalysis 
{
	//////////////////////////////////////////////////////////////////////////////
	// ProtoDUNE utilites
	//////////////////////////////////////////////////////////////////////////////

	struct ProtoDUNEUtils
	{
		protoana::ProtoDUNEPFParticleUtils &  pfputil;
		protoana::ProtoDUNETrackUtils &       trackutil;
		protoana::ProtoDUNETruthUtils &       truthutil;
	};

	ProtoDUNEUtils GetProtoDUNEUtils ();

	struct DUNEUtils
	{
		const geo::Geometry *                     geom;
		const detinfo::DetectorProperties *       detprop;
		const detinfo::DetectorClocks *           detclock;
		const cheat::BackTracker *                bt;
		const cheat::ParticleInventory *          pis;
		const sim::LArG4Parameters *              larg4;
	};

	DUNEUtils GetDUNEUtils ();

	//////////////////////////////////////////////////////////////////////////////
	// Hit data manipulation
	//////////////////////////////////////////////////////////////////////////////

	struct HitData
	{
		float PeakTime;
		float Integral;
		float Energy;

		double X;
		double Y;
		double Z;

		unsigned int CNNIndex_Michel;
		unsigned int CNNIndex_EM;

		float CNNScore_Michel;
		float CNNScore_EM;

		float CLUScore_Michel;
		float CLUScore_EM;

	};

	HitData GetHitData(const art::Ptr<recob::Hit> hit, 
	                   const mva_hit_t & hitcnnscores);

	HitData GetHitData(const art::Ptr<recob::Hit> hit, 
	                   double t0,
	                   const calo::CalorimetryAlg & caloAlg,
	                   const DUNEUtils & duneUtils,
	                   const art::Event & event,
	                   const mva_hit_t & hitcnnscores);

	HitData GetHitData(const art::Ptr<recob::Hit> hit, 
	                   double t0,
	                   const calo::CalorimetryAlg & caloAlg,
	                   const MichelAnalysis::DUNEUtils & duneUtils,
	                   const art::Event & event,
	                   const vec_ptr_sp_t spacePoint,
	                   const mva_hit_t & hitcnnscores);

	double HitToEnergy(art::Ptr<recob::Hit> hit, 
	                   double t0, 
	                   const calo::CalorimetryAlg & caloAlg,
	                   const DUNEUtils & duneUtils); 

	double HitToEnergy(art::Ptr<recob::Hit> hit, 
	                   double X, double Y, double Z,
	                   double t0, 
	                   const calo::CalorimetryAlg & caloAlg,
	                   const art::Event & event,
	                   const DUNEUtils & duneUtils); 

	bool HitsAreMichel(std::vector<art::Ptr<recob::Hit>> michelHits, 
                     const DUNEUtils & duneUtils);

	double HitToTrueEnergy(art::Ptr<recob::Hit> hit, const DUNEUtils & duneUtils);

	const simb::MCParticle * HitToMCP(art::Ptr<recob::Hit> hit, 
	                                  const DUNEUtils & duneUtils);
	bool HitIsMichel(const art::Ptr<recob::Hit> hit, 
	                 const DUNEUtils & duneUtils);

	bool MCPIsMichel(const simb::MCParticle * mcParticle);

	const simb::MCParticle * HitsToMcParticle(const vec_ptr_hit_t & hits, 
	                                          const DUNEUtils & duneUtils);

	//////////////////////////////////////////////////////////////////////////////
	// Track data manipulation
	//////////////////////////////////////////////////////////////////////////////

	struct TrackData
	{
		int     ID;

		bool    IsMuonDecaying;

		double  VertexX;
		double  VertexY;
		double  VertexZ;

		double  EndX;
		double  EndY;
		double  EndZ;

		bool    HasT0;
		double  T0;
	};

	void PrintTrackData(const TrackData & trackData);

	TrackData GetTrackData(const recob::Track * track, 
	                       const ProtoDUNEUtils & pdUtils, 
	                       const art::Event & event,
	                       const std::string trackModule);

	bool TrackIsDecayingMuon(const recob::Track * track, 
	                         const art::Event & event,
	                         const ProtoDUNEUtils & pdUtils,  
	                         const std::string trackModule, 
	                         const pid_map_t & particleMap);

	//////////////////////////////////////////////////////////////////////////////
	// Daughter data manipulation
	//////////////////////////////////////////////////////////////////////////////

	struct DaughterData
	{

		unsigned int       NHits;
		unsigned int       NMichelHits;
		float              FractionMichelHits; 

		// Data for this daughter
		std::vector<float> MichelScores;
		std::vector<float> EMScores;
		std::vector<float> Integrals;
		std::vector<float> Energies;
		std::vector<float> PeakTimes;

		float AvgMichelScore;
		float AvgEMScore;

		float DistanceToPrimary;
		float DistanceToPrimaryX;
		float DistanceToPrimaryY;
		float DistanceToPrimaryZ;

	};

	DaughterData ResetDaughterData();

	void AnalyseDaughterShowerHits(DaughterData & daughterData, 
	                               HitData & hitData,
	                               const DUNEUtils & duneUtils, 
	                               const vec_ptr_hit_t & daughterShowerHits,
	                               const mva_hit_t & hitcnnscores, 
	                               TTree * hitTree);

	void AnalyseDaughterShowerHits(DaughterData & daughterData, 
	                               HitData & hitData,
	                               const DUNEUtils & duneUtils, 
	                               const vec_ptr_hit_t & daughterShowerHits,
	                               const mva_hit_t & hitcnnscores);

	void AnalyseDaughterShowerHits(DaughterData & daughterData, 
	                               HitData & hitData, 
	                               const DUNEUtils & duneUtils,
	                               const vec_ptr_hit_t & daughterShowerHits,
	                               const mva_hit_t & hitcnnscores, 
	                               double t0,
	                               const calo::CalorimetryAlg & caloAlg,
	                               const art::Event & event,
	                               TTree * hitTree); 

	void AnalyseDaughterShowerHits(DaughterData & daughterData, 
	                               HitData & hitData, 
	                               const DUNEUtils & duneUtils,
	                               const vec_ptr_hit_t & daughterShowerHits,
	                               const mva_hit_t & hitcnnscores, 
	                               double t0,
	                               const calo::CalorimetryAlg & caloAlg,
	                               const art::Event & event); 

	//////////////////////////////////////////////////////////////////////////////
	// Data obtainers
	//////////////////////////////////////////////////////////////////////////////

	const vec_ptr_hit_t GetRecoTrackHits(const recob::Track & track,  
	                                     const art::Event & evt, 
	                                     const std::string trackModule);

	const vec_ptr_hit_t GetRecoShowerHits(const recob::Shower & shower, 
	                                      const art::Event & evt, 
                                        const std::string showerModule);

	void ResizeWireDriftData(std::vector<std::vector<float>> & data, 
	                         size_t nWires, 
	                         size_t nDrifts);

	void SetWireData(std::vector<std::vector<float>> & data, 
	                 const std::vector<float> & adcValues, 
	                 size_t wireIndex);

	void GetWireDriftData(std::vector<std::vector<float>> & data, 
	                     const std::vector<recob::Wire> & wires, 
	                     const DUNEUtils & duneUtils, 
	                     size_t plane, 
	                     size_t tpc, 
	                     size_t cryo);

	bool FillPatch (std::vector<std::vector<float>> & patch, 
	                const std::vector<std::vector<float>> & data, 
	                size_t wire, 
	                float drift, 
	                size_t patchSizeWire, 
	                size_t patchSizeDrift);

	//////////////////////////////////////////////////////////////////////////////
	// PFP manipulators
	//////////////////////////////////////////////////////////////////////////////

	bool IsPrimaryTrack(const recob::PFParticle & pfparticle, 
	                    const std::string particleLabel,
	                    const std::string trackLabel,
	                    const art::Event & event,
	                    const ProtoDUNEUtils & pDUtils);

	bool IsGoodPrimaryTrack(const recob::Track * track);

	const vec_ptr_hit_t GetPFParticleHits(const recob::PFParticle & particle, 
	                                      art::Event const &evt, 
	                                      const ProtoDUNEUtils & pdUtils,
	                                      const std::string particleLabel); 

	std::vector<int> GetPFPDaughtersRecursively(
	                   const recob::PFParticle & pfparticle, 
	                   const vec_pfp_t & pfparticles);

	//////////////////////////////////////////////////////////////////////////////
	// True Particle Manipulators
	//////////////////////////////////////////////////////////////////////////////

	int CountDecayingMuons(const vec_mcp_t & particles, 
	                       const pid_map_t & particleMap, 
	                       const DUNEUtils & duneUtils, 
	                       double fidVolCut);

	bool IsMuonDecaying(const mcp_t & particle, 
	                    const pid_map_t & particleMap);

	bool InsideFidVol(const double position[3], 
	                  const DUNEUtils & duneUtils, 
	                  double fidVolCut);

	pid_map_t BuildParticleMap(const vh_vec_mcp_t & particles);

	//////////////////////////////////////////////////////////////////////////////
	// Event Selection 
	//////////////////////////////////////////////////////////////////////////////

	struct EventSelectionData
	{
		int NTrue;
		int NSelected;
		int NCorrectlySelected;
		int NIncorrectlySelected;
	};

	EventSelectionData ResetEventSelectionData();

	bool EventSelection(const DaughterData & daughterData);

	bool EventSelection(const recob::Track & primaryTrack, 
	                    const recob::Shower & daughterShower,
	                    DaughterData & daughterData);

	bool EventSelection(const recob::Track & primaryTrack, 
	                    const recob::Track & daughterTrack,
	                    DaughterData & daughterData);

	void EventSelectionTests(EventSelectionData & eventSelectionData, 
	                         TrackData & primaryData, 
	                         const ProtoDUNEUtils & pdUtils, 
	                         const recob::Track * primarytrack, 
	                         const pid_map_t & particleMap, 
	                         const art::Event & event, 
	                         const std::string trackTag); 

	//////////////////////////////////////////////////////////////////////////////
	// Plotting Functions
	//////////////////////////////////////////////////////////////////////////////

	const float GetADCNorm(const ProtoDUNEUtils & pdUtils, 
	                       const recob::Track * track,
	                       const art::Event & event,
	                       const std::string trackLabel,
	                       const std::string caloLabel);

	TH2D DrawPatchAsHist(const std::vector<std::vector<float>> & patch, 
	                     const std::string basename, 
	                     art::ServiceHandle<art::TFileService> tfs);

	TH2D DrawMichelAsRawImage(const vec_ptr_hit_t & hits, 
	                          const art::Event & event,  
	                          const std::string basename, 
	                          art::ServiceHandle<art::TFileService> tfs, 
	                          const size_t nWireBins,
	                          const DUNEUtils & duneUtils,
	                          const std::string wireLabel);

	std::vector<TH2D> DrawMichelTrainingData(
	  const vec_ptr_hit_t & daughterHits,
	  const vec_pfp_t & pfparticles,
	  const std::string particleLabel,
	  const mva_hit_t & hitcnnscores,
	  const art::Event & event,
	  const std::string basename,
	  art::ServiceHandle<art::TFileService> tfs,
	  const size_t nWireBins,
	  const float norm,
	  const MichelAnalysis::ProtoDUNEUtils & pdUtils,
	  const MichelAnalysis::DUNEUtils & duneUtils);

	std::vector<TH2D> DrawMichelTrainingData(
	  const vec_ptr_hit_t & daughterHits, 
	  const vh_vec_clu_t & allClusters,
	  const std::string particleLabel,
	  const mva_hit_t & hitcnnscores,
	  const art::Event & event,
	  const std::string basename,
	  art::ServiceHandle<art::TFileService> tfs,
	  const size_t nWireBins,
	  const float norm,
	  const MichelAnalysis::ProtoDUNEUtils & pdUtils,
	  const MichelAnalysis::DUNEUtils & duneUtils);

	std::vector<TH2D> DrawMichelTrainingData(
	  const vec_ptr_hit_t & daughterHits,
	  const vec_ptr_hit_t & allHits,
	  const std::string particleLabel,
	  const mva_hit_t & hitcnnscores,
	  const art::Event & event,
	  const std::string basename,
	  art::ServiceHandle<art::TFileService> tfs,
	  const size_t nWireBins,
	  const float norm,
	  const MichelAnalysis::ProtoDUNEUtils & pdUtils,
	  const MichelAnalysis::DUNEUtils & duneUtils);

void DrawMichelTrainingData(
       const vec_ptr_hit_t & daughterHits,
       const vec_ptr_hit_t & allHits,
       const vec_pfp_t & pfparticles,
       const std::string particleLabel,
       const mva_hit_t & hitcnnscores,
       const double t0,
       const simb::MCParticle & mcParticle,
       const calo::CalorimetryAlg & caloAlg, 
       const art::Event & event,
       const std::string basename,
       art::ServiceHandle<art::TFileService> tfs,
       const size_t nWireBins,
       const float norm,
       const MichelAnalysis::ProtoDUNEUtils & pdUtils,
       const MichelAnalysis::DUNEUtils & duneUtils);

std::vector<TH2D> DrawMichelTrainingData_Data(
       const vec_ptr_hit_t & daughterHits,
       const vec_ptr_hit_t & allHits,
       const vec_pfp_t & pfparticles,
       const std::string particleLabel,
       const mva_hit_t & hitcnnscores,
       const double t0,
       const calo::CalorimetryAlg & caloAlg, 
       const art::Event & event,
       const std::string basename,
       art::ServiceHandle<art::TFileService> tfs,
       const size_t nWireBins,
       const float norm,
       const MichelAnalysis::ProtoDUNEUtils & pdUtils,
       const MichelAnalysis::DUNEUtils & duneUtils);
	//////////////////////////////////////////////////////////////////////////////
	// Analysis data
	//////////////////////////////////////////////////////////////////////////////

	struct MichelEvent
	{
		TrackData    Primary;
		DaughterData Daughter;
	};

} // namespace MichelAnalysis

#endif
