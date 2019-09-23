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
#include "lardataobj/RecoBase/Wire.h"
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
#include <iterator>
#include <stack>

// Michel analysis includes
#include "MichelAnalysis.h"

//////////////////////////////////////////////////////////////////////////////
// ProtoDUNE utilites
//////////////////////////////////////////////////////////////////////////////

// Returns a struct containing general protodune utilities
MichelAnalysis::ProtoDUNEUtils MichelAnalysis::GetProtoDUNEUtils ()
{

	protoana::ProtoDUNEPFParticleUtils pfputil;
	protoana::ProtoDUNETrackUtils      trackutil;
	protoana::ProtoDUNETruthUtils      truthutil;

	MichelAnalysis::ProtoDUNEUtils pdUtils = { pfputil, trackutil, truthutil };

	return pdUtils;

}

// Returns a struct containing general dune utilities
MichelAnalysis::DUNEUtils MichelAnalysis::GetDUNEUtils() 
{

	const geo::Geometry * geom = &* (art::ServiceHandle<geo::Geometry>());

	const sim::LArG4Parameters * larg4 = 
	  &* (art::ServiceHandle<sim::LArG4Parameters>());

	const detinfo::DetectorProperties * detprop = 
	  lar::providerFrom<detinfo::DetectorPropertiesService>();

	const detinfo::DetectorClocks * detclock = 
	  lar::providerFrom<detinfo::DetectorClocksService>();

	const cheat::BackTracker * bt = 
	  lar::providerFrom<cheat::BackTrackerService>();

	const cheat::ParticleInventory * pis = 
	  lar::providerFrom<cheat::ParticleInventoryService>();

	MichelAnalysis::DUNEUtils duneUtils = 
	  {geom, detprop, detclock, bt, pis, larg4};

	return duneUtils;

}

//////////////////////////////////////////////////////////////////////////////
// Hit data manipulation
//////////////////////////////////////////////////////////////////////////////

// Fills hit data struct with hit details
MichelAnalysis::HitData MichelAnalysis::GetHitData(
                          const art::Ptr<recob::Hit> hit, 
                          const mva_hit_t & hitcnnscores)
{

	MichelAnalysis::HitData hitData;

	hitData.PeakTime          = hit->PeakTime();
	hitData.Integral          = hit->Integral();
	hitData.Energy            = 0; 

	hitData.CNNIndex_Michel   = hitcnnscores.getIndex("michel");
	hitData.CNNIndex_EM       = hitcnnscores.getIndex("em");

	const std::array<float, 4> hitcnnscore = hitcnnscores.getOutput(hit);

	hitData.CNNScore_Michel   = hitcnnscore[hitData.CNNIndex_Michel];
	hitData.CNNScore_EM       = hitcnnscore[hitData.CNNIndex_EM];

	return hitData;

}

MichelAnalysis::HitData MichelAnalysis::GetHitData(
                          const art::Ptr<recob::Hit> hit, 
                          double t0,
                          const calo::CalorimetryAlg & caloAlg,
                          const MichelAnalysis::DUNEUtils & duneUtils,
                          const mva_hit_t & hitcnnscores)
{

	MichelAnalysis::HitData hitData = 
	  MichelAnalysis::GetHitData(hit, hitcnnscores);

	hitData.Energy = MichelAnalysis::HitToEnergy(hit, t0, caloAlg, duneUtils); 

	return hitData;

}

const simb::MCParticle * MichelAnalysis::HitsToMcParticle(
                           const vec_ptr_hit_t & hits, 
                           const MichelAnalysis::DUNEUtils & duneUtils)
{

	std::unordered_map<int, double> trackIDE;
	for (auto hit : hits)
	{
		auto const & ides = duneUtils.bt->HitToTrackIDEs(hit);
		for (auto const & ide : ides) { trackIDE[ide.trackID] += ide.energy; }
	}

	int bestID {0}; double totE {0}, maxE {0};
	for (auto const & contrib : trackIDE)
	{
		totE += contrib.second;
		if (contrib.second <= maxE) { continue; }
		maxE   = contrib.second;
		bestID = contrib.first;
	}

	const simb::MCParticle * mcParticle { nullptr };
	if ((maxE > 0) && (totE > 0))
	{
		if (bestID < 0) { bestID = - bestID; }
		mcParticle = duneUtils.pis->TrackIdToParticle_P(bestID);
	}

	return mcParticle;

}

bool MichelAnalysis::MCPIsMichel(const simb::MCParticle * mcParticle)
{

	if (mcParticle != nullptr) 
	{
		return abs(mcParticle -> PdgCode()) == 11 && (mcParticle -> Process() == 
			"Decay" || mcParticle -> Process() == "muMinusCaptureAtRest"); 
	}
	else { return false; }

}

double MichelAnalysis::HitToEnergy(art::Ptr<recob::Hit> hit, 
                                   double t0, 
                                   const calo::CalorimetryAlg & caloAlg,
                                   const MichelAnalysis::DUNEUtils & duneUtils) 
{

	constexpr int collectionView {2};

	if (hit->View() != collectionView) { return 0.; }

	double energy { 0. };

	double lifetime_correction = caloAlg.LifetimeCorrection(hit->PeakTime(), t0);

	double lifetime_corrected_integral { hit->Integral() * lifetime_correction };

	double electrons { caloAlg.ElectronsFromADCArea(lifetime_corrected_integral,
	                                                collectionView) }; 

	double electronsToMeV { 1000. / duneUtils.larg4->GeVToElectrons() };

	energy = electrons * electronsToMeV;

	double recombination_factor = 0.69; // TODO: This parameter needs more thought

	energy /= recombination_factor;

	return energy;

}

bool MichelAnalysis::HitsAreMichel(std::vector<art::Ptr<recob::Hit>> michelHits, 
                                   const DUNEUtils & duneUtils)
{

	std::unordered_map<int,double> trkIDE;
	for (auto const & hit : michelHits) {
		for (auto const & ide: duneUtils.bt->HitToTrackIDEs(hit)) {
			trkIDE[ide.trackID] += ide.energy;
		}
	}

	int best_id = 0; double tot_e = 0, max_e = 0;
	for (auto const & contrib : trkIDE) {
		tot_e += contrib.second;
		if (contrib.second <= max_e) {continue;}
		max_e = contrib.second;
		best_id = contrib.first;
	}
	
	const simb::MCParticle * mcParticle = nullptr;
	if (max_e > 0 && tot_e > 0) {
		if (best_id < 0) {
			best_id = -best_id;
		}
		mcParticle = duneUtils.pis->TrackIdToParticle_P(best_id);
	}

	if (mcParticle != nullptr) {
		return abs(mcParticle->PdgCode()) == 11 && (mcParticle->Process() == 
			"Decay" || mcParticle->Process() == "muMinusCaptureAtRest"); 
	}
	else { return false; }

}

double MichelAnalysis::HitToTrueEnergy(art::Ptr<recob::Hit> hit, 
                                       const DUNEUtils & duneUtils)
{

	double trueE = 0.;

	auto const & hitIDEs = duneUtils.bt->HitToAvgSimIDEs(hit);
	for (auto const & hitIDE : hitIDEs) { trueE += hitIDE.energy; }

	return trueE;

}

const simb::MCParticle * MichelAnalysis::HitToMCP(art::Ptr<recob::Hit> hit, 
                                                  const DUNEUtils & duneUtils)
{

	std::unordered_map<int,double> trkIDE;
	for (auto const & ide: duneUtils.bt->HitToTrackIDEs(hit)) {
		trkIDE[ide.trackID] += ide.energy;
	}

	int best_id = 0; double tot_e = 0, max_e = 0;
	for (auto const & contrib : trkIDE) {
		tot_e += contrib.second;
		if (contrib.second <= max_e) {continue;}
		max_e = contrib.second;
		best_id = contrib.first;
	}
	
	const simb::MCParticle * mcParticle = 0;
	if (max_e > 0 && tot_e > 0) {
		if (best_id < 0) { best_id = -best_id; }
		mcParticle = duneUtils.pis->TrackIdToParticle_P(best_id);
	}

	return mcParticle;

}

bool MichelAnalysis::HitIsMichel(art::Ptr<recob::Hit> hit, 
                                 const DUNEUtils & duneUtils)
{

	std::unordered_map<int,double> trkIDE;
	for (auto const & ide: duneUtils.bt->HitToTrackIDEs(hit)) {
		trkIDE[ide.trackID] += ide.energy;
	}

	int best_id = 0; double tot_e = 0, max_e = 0;
	for (auto const & contrib : trkIDE) {
		tot_e += contrib.second;
		if (contrib.second <= max_e) {continue;}
		max_e = contrib.second;
		best_id = contrib.first;
	}
	
	const simb::MCParticle * mcParticle = 0;
	if (max_e > 0 && tot_e > 0) {
		if (best_id < 0) { best_id = -best_id; }
		mcParticle = duneUtils.pis->TrackIdToParticle_P(best_id);
	}

	if (mcParticle != 0) {
		return abs(mcParticle->PdgCode()) == 11 && (mcParticle->Process() == 
			"Decay" || mcParticle->Process() == "muMinusCaptureAtRest"); 
	}
	else { return false; }

}

//////////////////////////////////////////////////////////////////////////////
// Track data manipulation
//////////////////////////////////////////////////////////////////////////////

// Prints data from track information struct
void MichelAnalysis::PrintTrackData(const MichelAnalysis::TrackData & trackData)
{

	std::cout << trackData.ID << std::endl; 

	std::cout << trackData.VertexX << ", " << trackData.VertexY << ", "
	          << trackData.VertexZ << std::endl;

	std::cout << trackData.EndX << ", " << trackData.EndY << ", "
	          << trackData.EndZ << std::endl;
}

// Fills track data struct with track details
MichelAnalysis::TrackData MichelAnalysis::GetTrackData(
                            const recob::Track * track,
                            const MichelAnalysis::ProtoDUNEUtils & pdUtils,
                            const art::Event & event,
                            const std::string trackModule)
{

	MichelAnalysis::TrackData trackData;

	trackData.ID      = track->ID();

	trackData.VertexX = track->Vertex().X();
	trackData.VertexY = track->Vertex().Y();
	trackData.VertexZ = track->Vertex().Z();

	trackData.EndX    = track->End().X();
	trackData.EndY    = track->End().Y();
	trackData.EndZ    = track->End().Z();

	std::vector<anab::T0> t0s = pdUtils.trackutil.GetRecoTrackT0(* track, event, 
	                                                             trackModule);

	trackData.HasT0 = (t0s.size() > 0);
	if (trackData.HasT0) { trackData.T0 = t0s[0].Time(); }
	else                 { trackData.T0 = 0; }

	return trackData;

}

bool MichelAnalysis::TrackIsDecayingMuon(
       const recob::Track * track, 
       const art::Event & event, 
       const MichelAnalysis::ProtoDUNEUtils & pdUtils, 
       const std::string trackModule, 
       const pid_map_t & particleMap)
{

	const simb::MCParticle * mcParticle = 
	  pdUtils.truthutil.GetMCParticleFromRecoTrack(* track, event, trackModule);

	if (mcParticle == nullptr) { return false; }
	return MichelAnalysis::IsMuonDecaying(* mcParticle, particleMap); 

}

////////////////////////////////////////////////////////////////////////////////
// Daughter data manipulation
////////////////////////////////////////////////////////////////////////////////

MichelAnalysis::DaughterData MichelAnalysis::ResetDaughterData()
{
	MichelAnalysis::DaughterData daughterData;
	return daughterData;
}

// Analyses the daughter hits from showers
//   - Fills daughter data struct
//   - Fills hit information into hit tree
void MichelAnalysis::AnalyseDaughterShowerHits(
       MichelAnalysis::DaughterData & daughterData, 
       MichelAnalysis::HitData & hitData, 
       const MichelAnalysis::DUNEUtils & duneUtils,
       const vec_ptr_hit_t & daughterShowerHits,
       const mva_hit_t & hitcnnscores, 
       TTree * hitTree) 
{

	daughterData.NMichelHits  = 0;
	daughterData.NHits        = daughterShowerHits.size();

	daughterData.AvgMichelScore = 0.;
	daughterData.AvgEMScore     = 0.;
	for (auto const & daughterShowerHit : daughterShowerHits) 
	{

		hitData = MichelAnalysis::GetHitData(daughterShowerHit, hitcnnscores);

		daughterData.MichelScores.push_back(hitData.CNNScore_Michel);
		daughterData.EMScores.push_back(hitData.CNNScore_EM);
		daughterData.Integrals.push_back(hitData.Integral);
		daughterData.PeakTimes.push_back(hitData.PeakTime);
		daughterData.Energies.push_back(hitData.Energy);

		daughterData.AvgMichelScore += hitData.CNNScore_Michel;
		daughterData.AvgEMScore += hitData.CNNScore_EM;

		hitTree ->Fill();

		if (hitData.CNNScore_Michel > 0.9) { daughterData.NMichelHits += 1; }

	}

	if (daughterData.NHits != 0) 
	{  
		daughterData.AvgMichelScore /= daughterData.NHits;
		daughterData.AvgEMScore     /= daughterData.NHits;
	}

	daughterData.FractionMichelHits = 
	  static_cast<float>(daughterData.NMichelHits) /
	  static_cast<float>(daughterData.NHits);

}

void MichelAnalysis::AnalyseDaughterShowerHits(
       MichelAnalysis::DaughterData & daughterData, 
       MichelAnalysis::HitData & hitData, 
       const MichelAnalysis::DUNEUtils & duneUtils,
       const vec_ptr_hit_t & daughterShowerHits,
       const mva_hit_t & hitcnnscores) 
{

	daughterData.NMichelHits  = 0;
	daughterData.NHits        = daughterShowerHits.size();

	for (auto const & daughterShowerHit : daughterShowerHits) 
	{

		hitData = MichelAnalysis::GetHitData(daughterShowerHit, hitcnnscores);

		daughterData.MichelScores.push_back(hitData.CNNScore_Michel);
		daughterData.EMScores.push_back(hitData.CNNScore_EM);
		daughterData.Integrals.push_back(hitData.Integral);
		daughterData.PeakTimes.push_back(hitData.PeakTime);
		daughterData.Energies.push_back(hitData.Energy);

		if (hitData.CNNScore_Michel > 0.9) { daughterData.NMichelHits += 1; }

	}

	daughterData.FractionMichelHits = 
	  static_cast<float>(daughterData.NMichelHits) /
	  static_cast<float>(daughterData.NHits);

}

void MichelAnalysis::AnalyseDaughterShowerHits(
       MichelAnalysis::DaughterData & daughterData, 
       MichelAnalysis::HitData & hitData, 
       const MichelAnalysis::DUNEUtils & duneUtils,
       const vec_ptr_hit_t & daughterShowerHits,
       const mva_hit_t & hitcnnscores, 
       double t0,
       const calo::CalorimetryAlg & caloAlg,
       TTree * hitTree) 
{

	daughterData.NMichelHits  = 0;
	daughterData.NHits        = daughterShowerHits.size();

	for (auto const & daughterShowerHit : daughterShowerHits) 
	{

		hitData = MichelAnalysis::GetHitData(daughterShowerHit, t0, caloAlg,
		                                     duneUtils, hitcnnscores);

		daughterData.MichelScores.push_back(hitData.CNNScore_Michel);
		daughterData.EMScores.push_back(hitData.CNNScore_EM);
		daughterData.Integrals.push_back(hitData.Integral);
		daughterData.PeakTimes.push_back(hitData.PeakTime);
		daughterData.Energies.push_back(hitData.Energy);

		hitTree ->Fill();

		if (hitData.CNNScore_Michel > 0.9) { daughterData.NMichelHits += 1; }

	}

	daughterData.FractionMichelHits = 
	  static_cast<float>(daughterData.NMichelHits) /
	  static_cast<float>(daughterData.NHits);

}

void MichelAnalysis::AnalyseDaughterShowerHits(
       MichelAnalysis::DaughterData & daughterData, 
       MichelAnalysis::HitData & hitData, 
       const MichelAnalysis::DUNEUtils & duneUtils,
       const vec_ptr_hit_t & daughterShowerHits,
       const mva_hit_t & hitcnnscores, 
       double t0,
       const calo::CalorimetryAlg & caloAlg) 
{

	daughterData.NMichelHits  = 0;
	daughterData.NHits        = daughterShowerHits.size();

	for (auto const & daughterShowerHit : daughterShowerHits) 
	{

		hitData = MichelAnalysis::GetHitData(daughterShowerHit, t0, caloAlg,
		                                     duneUtils, hitcnnscores);

		daughterData.MichelScores.push_back(hitData.CNNScore_Michel);
		daughterData.EMScores.push_back(hitData.CNNScore_EM);
		daughterData.Integrals.push_back(hitData.Integral);
		daughterData.PeakTimes.push_back(hitData.PeakTime);
		daughterData.Energies.push_back(hitData.Energy);

		if (hitData.CNNScore_Michel > 0.9) { daughterData.NMichelHits += 1; }

	}

	daughterData.FractionMichelHits = 
	  static_cast<float>(daughterData.NMichelHits) /
	  static_cast<float>(daughterData.NHits);

}

////////////////////////////////////////////////////////////////////////////////
// Data obtainers
////////////////////////////////////////////////////////////////////////////////

// Gets hits from reconstructed track
const vec_ptr_hit_t MichelAnalysis::GetRecoTrackHits(
                      const recob::Track & track, 
                      const art::Event & event, 
                      const std::string trackModule) 
{

	auto recoTracks =  
	  event.getValidHandle<std::vector<recob::Track>>(trackModule);

	art::FindManyP<recob::Hit> findHits(recoTracks, event, trackModule);
	vec_ptr_hit_t inputHits = findHits.at(track.ID());

	return inputHits;

}

// Gets hits from reconstructed shower
const vec_ptr_hit_t MichelAnalysis::GetRecoShowerHits(
                      const recob::Shower & shower, 
                      const art::Event & event, 
                      const std::string showerModule) 
{

	auto recoShowers = 
	  event.getValidHandle<std::vector<recob::Shower>>(showerModule);

	art::FindManyP<recob::Hit> findHits(recoShowers, event, showerModule);

	int actualIndex = shower.ID();
	if(shower.ID() < 0)
	{
		for(unsigned int s = 0; s < recoShowers->size(); ++s)
		{
			const recob::Shower thisShower = (* recoShowers)[s];
			if(fabs(thisShower.Length() - shower.Length()) < 1e-5)
			{ 
				actualIndex = s; 
				continue; 
			}
		}
	}

	vec_ptr_hit_t inputHits = findHits.at(actualIndex);

	return inputHits;

}

// Sets the patch size for wire drift data
void MichelAnalysis::ResizeWireDriftData(std::vector<std::vector<float>> & data,
                                         size_t nWires, 
                                         size_t nDrifts)
{

	data.resize(nWires);

	for (auto & wire : data)
	{

		wire.resize(nDrifts);
		std::fill(wire.begin(), wire.end(), 0);

	}

}

// Fills the patch with the wire drift data
void MichelAnalysis::SetWireData(std::vector<std::vector<float>> & data, 
                                 const std::vector<float> & adcValues,
                                 size_t wireIndex) 
{

	if (wireIndex >= data.size()) { return; }
	if (adcValues.empty()) { return; }

	auto & wireData = data[wireIndex];

	if (adcValues.size() <= wireData.size())
	{
		std::copy(adcValues.begin(), adcValues.end(), wireData.begin());
	}
	else
	{
		std::copy(adcValues.begin(), adcValues.begin() + wireData.size(), 
		          wireData.begin());
	}

}

// Gets hold of the wire drift data for a given plane, tpc and cryostat
void MichelAnalysis::GetWireDriftData(
       std::vector<std::vector<float>> & data, 
       const std::vector<recob::Wire> & wires, 
       const MichelAnalysis::DUNEUtils & duneUtils, 
       size_t plane, 
       size_t tpc, 
       size_t cryo)
{

	size_t nWires   = duneUtils.geom->Nwires(plane, tpc, cryo);
	size_t nDrifts  = duneUtils.detprop->NumberTimeSamples(); 

	MichelAnalysis::ResizeWireDriftData(data, nWires, nDrifts);

	for (auto const & wire : wires)
	{
		auto wireChannel = wire.Channel();
		size_t wireIndex = 0; 
		for (auto const & wireID : duneUtils.geom->ChannelToWire(wireChannel))
		{
			if ((wireID.Plane == plane) && (wireID.TPC == tpc) && 
			    (wireID.Cryostat == cryo))
			{
				wireIndex = wireID.Wire;
				const auto adcValues = wire.Signal();
				MichelAnalysis::SetWireData(data, adcValues, wireIndex);
			}
		}
	}

}

// Fills the patch with the wire drift data
//   - Zero padding on areas outside of initial plane
//   - Aspect ratio based on drift vs wire bin sizes
//   - TODO: implement scaling based on drift binning vs wire binning
bool MichelAnalysis::FillPatch(std::vector<std::vector<float>> & patch, 
                               const std::vector<std::vector<float>> & data,
                               size_t wire, 
                               float drift, 
                               size_t patchSizeWire,
                               size_t patchSizeDrift)
{

	patch.resize(patchSizeWire);

	int patchHalfSizeWire  = patchSizeWire / 2;
	int patchHalfSizeDrift = patchSizeDrift / 2;

	int wireLow   = wire - patchHalfSizeWire;
	int wireHigh  = wire + patchHalfSizeWire;

	int driftLow   = static_cast<int>(drift) - patchHalfSizeDrift;
	int driftHigh  = static_cast<int>(drift) + patchHalfSizeDrift;

	int wireSize = data.size();

	float maxADC = 0;

	for (int wireIndex = wireLow; wireIndex < wireHigh; ++wireIndex)
	{

		int wirePatchIndex = wireIndex - wireLow;

		auto & patchDriftData = patch[wirePatchIndex];
		patchDriftData.resize(patchSizeDrift); 

		if ((wireIndex >= 0) && (wireIndex < wireSize))
		{

			auto & sourceDriftData = data[wireIndex];

			int driftSize = sourceDriftData.size();

			for (int driftIndex = driftLow; driftIndex < driftHigh; ++driftIndex)
			{

				int driftPatchIndex = driftIndex - driftLow;

				if ((driftIndex >= 0) && (driftIndex < driftSize))
				{
					auto adc = sourceDriftData[driftIndex];
					if (adc > maxADC) { maxADC = adc; }
					patchDriftData[driftPatchIndex] = adc;
				}
				else
				{
					patchDriftData[driftPatchIndex] = 0;
				}

			}

		}
		else
		{
			std::fill(patchDriftData.begin(), patchDriftData.end(), 0);
		}

	}

	if (maxADC < 1E-9) { return false; }
	return true;

}

////////////////////////////////////////////////////////////////////////////////
// PFP manipulators
////////////////////////////////////////////////////////////////////////////////

// Checks if a given pfparticle is:
//   - Track like
//   - Primary
bool MichelAnalysis::IsPrimaryTrack(
       const recob::PFParticle & pfparticle, 
       const std::string particleLabel,
       const std::string trackLabel,
       const art::Event & event,
       const MichelAnalysis::ProtoDUNEUtils & pdUtils)
{

	bool primaryTrack {true};

	if (!pfparticle.IsPrimary()) { primaryTrack = false; }
	if (!pdUtils.pfputil.IsPFParticleTracklike(pfparticle, event, particleLabel, 
	                                           trackLabel)) 
	{ 
		primaryTrack = false; 
	}

	return primaryTrack;

}

// Checks if a given track is good
bool MichelAnalysis::IsGoodPrimaryTrack(const recob::Track * track)
{
	bool goodTrack {true};
	if (track == nullptr) { goodTrack = false; }
	return goodTrack;
}

const vec_ptr_hit_t MichelAnalysis::GetPFParticleHits(
                      const recob::PFParticle & particle, 
                      art::Event const & evt, 
                      const MichelAnalysis::ProtoDUNEUtils & pdUtils,
                      const std::string particleLabel) 
{

	const std::vector<const recob::Cluster*> pfpClusters = 
	  pdUtils.pfputil.GetPFParticleClusters(particle, evt, particleLabel);

	auto allClusters = 
	  evt.getValidHandle<std::vector<recob::Cluster>>(particleLabel);

	const art::FindManyP<recob::Hit> findHits(allClusters, evt, particleLabel);

	vec_ptr_hit_t pfpHits;
	for(auto const & cluster : pfpClusters){

		const vec_ptr_hit_t & clusterHits = findHits.at(cluster->ID());

		pfpHits.reserve(pfpHits.size() + clusterHits.size());
		pfpHits.insert(pfpHits.end(), 
		               std::make_move_iterator(clusterHits.begin()), 
		               std::make_move_iterator(clusterHits.end()));

	}

	return pfpHits;

	// const std::vector<const recob::SpacePoint*> spacePoints = 
	//   pdUtils.pfputil.GetPFParticleSpacePoints(particle, evt, particleLabel);

	// auto allSpacePoints = 
	//   evt.getValidHandle<std::vector<recob::SpacePoint>>(particleLabel);

	// const art::FindManyP<recob::Hit> findHits(allSpacePoints, evt, particleLabel);

	// vec_ptr_hit_t pfpHits;
	// for(auto const & sp : spacePoints)
	// {
	// 	const vec_ptr_hit_t & spacePointHits = findHits.at(sp->ID());
	// 	
	// 	pfpHits.reserve(pfpHits.size() + spacePointHits.size());
	// 	pfpHits.insert(pfpHits.end(), 
	// 	               std::make_move_iterator(spacePointHits.begin()), 
	// 	               std::make_move_iterator(spacePointHits.end()));
	// }

	// return pfpHits;

}

std::vector<int> MichelAnalysis::GetPFPDaughtersRecursively(
                   const recob::PFParticle & pfparticle, 
                   const vec_pfp_t & pfparticles)
{

	std::vector<int> ids;

	std::stack<int> daughterIdStack;
	for (auto const & daughterID : pfparticle.Daughters())
	{
		daughterIdStack.push(daughterID);
	}

	while (!daughterIdStack.empty())
	{

		auto const & daughter = pfparticles[daughterIdStack.top()];
		daughterIdStack.pop();

		ids.push_back(daughter.Self());

		for (auto const & subDaughterID : daughter.Daughters())
		{
			daughterIdStack.push(subDaughterID);
		}

	}

	return ids;

}


////////////////////////////////////////////////////////////////////////////////
// True Particle Manipulators
////////////////////////////////////////////////////////////////////////////////

// Counts the number of true decaying muons in the given mc particles
//   - Map is required for significant efficiency improvement
int MichelAnalysis::CountDecayingMuons(
      const vec_mcp_t & particles, 
      const pid_map_t & particleMap, 
      const MichelAnalysis::DUNEUtils & duneUtils, 
      double fidVolCut)
{
	int nMichel = 0;

	for ( auto const & particle : particles ) 
	{

		// Only consider muons and check if the decay
		if (abs(particle.PdgCode()) == 13) 
		{

			// Use the particles end point and check if it is inside the fiducial 
			// volume
			double mcEnd[3] = { particle.EndX(), particle.EndY(), particle.EndZ() };

			if (MichelAnalysis::IsMuonDecaying(particle, particleMap) && 
			    MichelAnalysis::InsideFidVol(mcEnd, duneUtils, fidVolCut)) 
			{
				nMichel += 1;
			}

		}

	}

	return nMichel;

}

// Checks if a given mc particle is a decaying muon
bool MichelAnalysis::IsMuonDecaying(const mcp_t & particle, 
                                    const pid_map_t & particleMap)
{

	if ( abs(particle.PdgCode()) != 13) { return false; }

	if (particle.EndProcess() != "FastScintillation") { return false; }
	
	bool hasElectron(false), hasNuMu(false), hasNuE(false);
	
	size_t n_daughters = particle.NumberDaughters();
	for (size_t d = 0; d < n_daughters; d++) {

		auto daughter_search = particleMap.find(particle.Daughter(d));
		if (daughter_search == particleMap.end()) { continue; }

		auto const & daughter = * ((* daughter_search).second);
		int daughterPDG = abs(daughter.PdgCode());

		if (daughterPDG == 11) { hasElectron = true; }
		if (daughterPDG == 14) { hasNuMu = true; }
		if (daughterPDG == 12) { hasNuE = true; }

	}

	return (hasElectron && hasNuMu && hasNuE);

}

// Checks if a given 3d point is inside the fiducial volume
bool MichelAnalysis::InsideFidVol(const double position [3], 
                                  const DUNEUtils & duneUtils,
                                  double fidVolCut)
{

	bool inside = false;

	geo::TPCID idtpc = duneUtils.geom->FindTPCAtPosition(position);

	if (duneUtils.geom->HasTPC(idtpc)) 
	{

		const geo::TPCGeo & tpcgeo = duneUtils.geom->GetElement(idtpc);
		double minx = tpcgeo.MinX(); double maxx = tpcgeo.MaxX();
		double miny = tpcgeo.MinY(); double maxy = tpcgeo.MaxY();
		double minz = tpcgeo.MinZ(); double maxz = tpcgeo.MaxZ();
		
		for (size_t cryo = 0; cryo < duneUtils.geom->Ncryostats(); cryo++) 
		{

			const geo::CryostatGeo & cryostat = duneUtils.geom->Cryostat(cryo);

			for (size_t tpc = 0; tpc < cryostat.NTPC(); tpc++) 
			{

				const geo::TPCGeo & tpcg = cryostat.TPC(tpc);

				if (tpcg.MinX() < minx) { minx = tpcg.MinX(); }
				if (tpcg.MaxX() > maxx) { maxx = tpcg.MaxX(); }

				if (tpcg.MinY() < miny) { miny = tpcg.MinY(); }
				if (tpcg.MaxY() > maxy) { maxy = tpcg.MaxY(); }

				if (tpcg.MinZ() < minz) { minz = tpcg.MinZ(); }
				if (tpcg.MaxZ() > maxz) { maxz = tpcg.MaxZ(); }

			}
		}

		double dista = fabs(minx - position[0]), distb = fabs(position[0] - maxx); 
		if ((position[0] > minx) && (position[0] < maxx) && 
				(dista > fidVolCut) && (distb > fidVolCut)) 
		{ 
			inside = true; 
		}

		dista = fabs(miny - position[1]); distb = fabs(position[1]-maxy);
		if (inside && (position[1] > miny) && (position[1] < maxy) &&
				(dista > fidVolCut) && (distb > fidVolCut)) 
		{
			inside = true;
		}
		else { inside = false; }

		dista = fabs(minz - position[2]); distb = fabs(position[2] - maxz);
		if (inside && (position[2] > minz) && (position[2] < maxz) && 
		    (dista > fidVolCut) && (distb > fidVolCut)) 
		{ 
			inside = true; 
		}
		else { inside = false; }

	}

	return inside;

}

// Builds a map between particle id and mc particle
//   - Provides significant efficiency improvement when doing particle searches 
//     for truth matching
pid_map_t MichelAnalysis::BuildParticleMap(const vh_vec_mcp_t & particles)
{

	pid_map_t particleMap;

	for (auto const & particle : * particles) 
	{
		particleMap[particle.TrackId()] = &particle;
	}
	
	return particleMap;

}

//////////////////////////////////////////////////////////////////////////////
// Event Selection Functions
//////////////////////////////////////////////////////////////////////////////

// Sets the event selection information in the current event to empty
MichelAnalysis::EventSelectionData MichelAnalysis::ResetEventSelectionData()
{
	MichelAnalysis::EventSelectionData data = {0, 0, 0, 0};
	return data;
}

// Performs events selection based on the data from the analysis of a daughter 
// particle
bool MichelAnalysis::EventSelection(const DaughterData & daughterData)
{
	bool eventSelected = { true };
	if (daughterData.FractionMichelHits < 0.8) { eventSelected = false; }
	return eventSelected;
}

// Event selection:
//   - Daughter particle is michel like
//   - Daughter is close to primary track end
bool MichelAnalysis::EventSelection(const recob::Track & primaryTrack,
                                    const recob::Track & daughterTrack,
                                    DaughterData & daughterData)
{

	bool eventSelected = { MichelAnalysis::EventSelection(daughterData) };

	if (eventSelected == false) { return eventSelected; }

	// I had to construct these from scratch because the underlying object for 
	// primaryTrack().Vertex() got changed to a ROOT::Math::PositionVector3D
	const TVector3 primaryTrackVertex  (primaryTrack.Vertex().X(), 
			                                primaryTrack.Vertex().Y(),
			                                primaryTrack.Vertex().Z());

	const TVector3 primaryTrackEnd     (primaryTrack.End().X(), 
			                                primaryTrack.End().Y(),
			                                primaryTrack.End().Z());

	const TVector3 daughterTrackVertex (daughterTrack.Vertex().X(),
	                                    daughterTrack.Vertex().Y(),
	                                    daughterTrack.Vertex().Z());

	auto vertexVec = primaryTrackVertex - daughterTrackVertex;
	auto endVec    = primaryTrackEnd - daughterTrackVertex;

	auto vertexDist = vertexVec.Mag2();
	auto endDist    = endVec.Mag2();

	bool endClosest = (endDist < vertexDist);
	if (endClosest) 
	{ 
		daughterData.DistanceToPrimary  = endVec.Mag2(); 
		daughterData.DistanceToPrimaryX = endVec.X(); 
		daughterData.DistanceToPrimaryY = endVec.Y(); 
		daughterData.DistanceToPrimaryZ = endVec.Z(); 

		// TODO: this track end cut should be in a fcl
		if (endDist > 5) { eventSelected = false; }
	}
	else 
	{ 
		daughterData.DistanceToPrimary  = vertexVec.Mag2(); 
		daughterData.DistanceToPrimaryX = vertexVec.X(); 
		daughterData.DistanceToPrimaryY = vertexVec.Y(); 
		daughterData.DistanceToPrimaryZ = vertexVec.Z(); 

		// TODO: this track end cut should be in a fcl
		if (vertexDist > 5) { eventSelected = false; }
	}

	return eventSelected;

}

// Event selection:
//   - Daughter particle is michel like
//   - Daughter is close to primary track end
bool MichelAnalysis::EventSelection(const recob::Track & primaryTrack,
                                    const recob::Shower & daughterShower,
                                    DaughterData & daughterData)
{

	bool eventSelected = { MichelAnalysis::EventSelection(daughterData) };

	if (eventSelected == false) { return eventSelected; }


	// I had to construct these from scratch because the underlying object for 
	// primaryTrack().Vertex() got changed to a ROOT::Math::PositionVector3D
	const TVector3 primaryTrackVertex (primaryTrack.Vertex().X(), 
			                               primaryTrack.Vertex().Y(),
			                               primaryTrack.Vertex().Z());
	const TVector3 primaryTrackEnd    (primaryTrack.End().X(), 
			                            	 primaryTrack.End().Y(),
			                            	 primaryTrack.End().Z());

	const TVector3 daughterShowerStart (daughterShower.ShowerStart().X(), 
			                                daughterShower.ShowerStart().Y(),
			                                daughterShower.ShowerStart().Z());

	auto vertexVec = primaryTrackVertex - daughterShowerStart;
	auto endVec    = primaryTrackEnd - daughterShowerStart;

	auto vertexDist = vertexVec.Mag2();
	auto endDist    = endVec.Mag2();

	bool endClosest = (endDist < vertexDist);
	if (endClosest) 
	{ 
		daughterData.DistanceToPrimary  = endVec.Mag2(); 
		daughterData.DistanceToPrimaryX = endVec.X(); 
		daughterData.DistanceToPrimaryY = endVec.Y(); 
		daughterData.DistanceToPrimaryZ = endVec.Z(); 

		// TODO: this track end cut should be in a fcl
		if (endDist > 5) { eventSelected = false; }
	}
	else 
	{ 
		daughterData.DistanceToPrimary  = vertexVec.Mag2(); 
		daughterData.DistanceToPrimaryX = vertexVec.X(); 
		daughterData.DistanceToPrimaryY = vertexVec.Y(); 
		daughterData.DistanceToPrimaryZ = vertexVec.Z(); 

		// TODO: this track end cut should be in a fcl
		if (vertexDist > 5) { eventSelected = false; }
	}

	return eventSelected;

}

// Efficiency and purity tests on selected michel electrons
void MichelAnalysis::EventSelectionTests(
  MichelAnalysis::EventSelectionData & eventSelectionData, 
  MichelAnalysis::TrackData & primaryData,
  const MichelAnalysis::ProtoDUNEUtils & pdUtils, 
  const recob::Track * primarytrack, 
  const pid_map_t & particleMap, 
  const art::Event & event, 
  std::string trackTag)  
{

	eventSelectionData.NSelected += 1;

	// TODO: I need a better event selection test for purity
	primaryData.IsMuonDecaying = MichelAnalysis::TrackIsDecayingMuon(primarytrack,
	                               event, pdUtils, trackTag, particleMap);

	if (primaryData.IsMuonDecaying) 
	{ 
		eventSelectionData.NCorrectlySelected += 1; 
	}
	else { eventSelectionData.NIncorrectlySelected += 1; }

}

//////////////////////////////////////////////////////////////////////////////
// Plotting Functions
//////////////////////////////////////////////////////////////////////////////

// Returns a 2D histogram of the raw data in a patch
TH2D MichelAnalysis::DrawPatchAsHist(
  const std::vector<std::vector<float>> & patch, 
  const std::string basename, 
  art::ServiceHandle<art::TFileService> tfs)
{

	size_t patchSizeWire = patch.size();
	if (patchSizeWire <= 0) 
	{
		TH2D patchHist = TH2D();
		return patchHist;
	}

	size_t patchSizeDrift = patch[0].size();

	TH2D * patchHist = tfs->make<TH2D>((basename + "_raw").c_str(), 
	                                   "Raw Data Image", patchSizeWire, 0, 
	                                   patchSizeWire, patchSizeDrift, 0, 
	                                   patchSizeDrift);

			for (size_t i = 0; i < patchSizeWire; ++i)
			{
				for (size_t j = 0; j < patch.at(i).size(); ++j) 
				{ 
					auto adc = patch[i][j];
					patchHist->Fill(static_cast<double>(i), static_cast<double>(j), adc);
				}
			}

			return *patchHist;

}

// Draws a given set of hits as a 2D histogram
TH2D MichelAnalysis::DrawMichelAsRawImage(
       const vec_ptr_hit_t & hits,
       const art::Event & event,
       const std::string basename, 
       art::ServiceHandle<art::TFileService> tfs,
       const size_t nWireBins,
       const MichelAnalysis::DUNEUtils & duneUtils)
{

	// Get details of one of the hits and use it as the centre of the image
	// try to use plane 2 if possible since that is the collection view
	art::Ptr<recob::Hit>  hit = hits[0];
	unsigned int plane = hit->WireID().Plane;
	size_t i = 0;
	while (plane != 2)
	{
		i += 1;
		if (i >= hits.size()) { break; }
		hit = hits[i];
		plane = hit->WireID().Plane;
	}

	const unsigned int wire     = hit->WireID().Wire;
	plane    = hit->WireID().Plane;
	const unsigned int tpc      = hit->WireID().TPC;
	const unsigned int cryostat = hit->WireID().Cryostat;
	const float        peakTime = hit->PeakTime();

	art::Handle<std::vector<recob::Wire>> wireHandle;
	event.getByLabel("caldata", wireHandle);

	std::vector<std::vector<float>> rawData;
	MichelAnalysis::GetWireDriftData(rawData, * wireHandle, duneUtils, plane,
	                                 tpc, cryostat);

	std::vector<std::vector<float>> patch;

	double driftVelocity   = duneUtils.detprop->DriftVelocity();
	double clockFrequency  = duneUtils.detclock->TPCClock().Frequency();

	double wireBinSize   = duneUtils.geom->WirePitch(plane);
	double timeBinSize   = driftVelocity / clockFrequency;

	double aspectRatio = wireBinSize / timeBinSize;

	// TODO: read the nBins from the fhicl file
	const size_t nTimeBins = static_cast<size_t>(
	  static_cast<double>(nWireBins) * aspectRatio);

	if (!MichelAnalysis::FillPatch(patch, rawData, wire, peakTime, nWireBins, 
	                          nTimeBins))
	{
		TH2D patchHist = TH2D();
		return patchHist;
	}

	TH2D patchHist = MichelAnalysis::DrawPatchAsHist(patch, basename, tfs);

	return patchHist;

}

std::vector<TH2D> MichelAnalysis::DrawMichelTrainingData(
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
       const MichelAnalysis::DUNEUtils & duneUtils)
{


	// Get details of one of the hits and use it as the centre of the image
	// try to use plane 2 if possible since that is the collection view
	unsigned int centralWire     { 0 };
	unsigned int centralTick     { 0 };
	unsigned int nCollectionHits { 0 };

	art::Ptr<recob::Hit>  hit;
	for (auto const & daughterHit : daughterHits)
	{

		if (daughterHit->WireID().Plane != 2) { continue; }

		if (hit.isNull()) { hit = daughterHit; }

		nCollectionHits += 1;
		centralWire     += daughterHit->WireID().Wire;
		centralTick     += daughterHit->PeakTime();

	}

	if (nCollectionHits < 1) 
	{ 
		std::cout << "No collection plane hits\n";
		std::vector<TH2D> trainingData;

		TH2D dataHist      = TH2D();
		TH2D energyHist    = TH2D();
		TH2D cnnEMHist     = TH2D();
		TH2D cnnMichelHist = TH2D();
		TH2D cluEMHist     = TH2D();
		TH2D cluMichelHist = TH2D();
		TH2D truthHist     = TH2D();

		trainingData.push_back(dataHist);
		trainingData.push_back(energyHist);
		trainingData.push_back(cnnEMHist);
		trainingData.push_back(cnnMichelHist);
		trainingData.push_back(cluEMHist);
		trainingData.push_back(cluMichelHist);
		trainingData.push_back(truthHist);

		return trainingData;

	}

	centralWire /= nCollectionHits;
	centralTick /= nCollectionHits;

	const unsigned int plane    { hit->WireID().Plane };
	const unsigned int tpc      { hit->WireID().TPC };
	const unsigned int cryostat { hit->WireID().Cryostat };

	double driftVelocity  { duneUtils.detprop->DriftVelocity() };
	double clockFrequency { duneUtils.detclock->TPCClock().Frequency() };
	double wireBinSize    { duneUtils.geom->WirePitch(plane) };
	double timeBinSize    { driftVelocity / clockFrequency };
	double aspectRatio    { wireBinSize / timeBinSize };
	
	// TODO: read the nBins from the fhicl file
	const size_t nTimeBins = static_cast<size_t>(
	  static_cast<double>(nWireBins) * aspectRatio);


	TH2D * dataHist      = tfs->make<TH2D>((basename + "_wire").c_str(), 
	                                     "Wire Time Image", nWireBins, 0, 
	                                     nWireBins, nWireBins, 0, nTimeBins);
	TH2D * energyHist    = tfs->make<TH2D>((basename + "_energy").c_str(), 
	                                     "Energy Time Image", nWireBins, 0, 
	                                     nWireBins, nWireBins, 0, nTimeBins);
	TH2D * cnnEMHist     = tfs->make<TH2D>((basename + "_cnnem").c_str(), 
	                                       "EM Score Image", nWireBins, 0, 
	                                       nWireBins, nWireBins, 0, nTimeBins);
	TH2D * cnnMichelHist = tfs->make<TH2D>((basename + "_cnnmichel").c_str(), 
	                                     "Michel Score Image", nWireBins, 0, 
	                                     nWireBins, nWireBins, 0, nTimeBins);
	TH2D * cluEMHist     = tfs->make<TH2D>((basename + "_cluem").c_str(), 
	                                       "EM Score Image", nWireBins, 0, 
	                                       nWireBins, nWireBins, 0, nTimeBins);
	TH2D * cluMichelHist = tfs->make<TH2D>((basename + "_clumichel").c_str(), 
	                                     "Michel Score Image", nWireBins, 0, 
	                                     nWireBins, nWireBins, 0, nTimeBins);
	TH2D * truthHist     = tfs->make<TH2D>((basename + "_truth").c_str(),
	                                     "Truth Image", nWireBins, 0, nWireBins,
	                                     nWireBins, 0, nTimeBins);

	for (auto const & pfparticle : pfparticles)
	{

		const vec_ptr_hit_t pfphits = 
		  MichelAnalysis::GetPFParticleHits(pfparticle, event, pdUtils, 
		                                    particleLabel);

		double cluScoreEM     = 0.;
		double cluScoreMichel = 0.;
		for (auto const & hit : pfphits) 
		{
			const HitData & hitData = GetHitData(hit, hitcnnscores);
			cluScoreEM     += hitData.CNNScore_EM;
			cluScoreMichel += hitData.CNNScore_Michel;
		}

		if (pfphits.size() != 0)
		{
			cluScoreEM     /= pfphits.size();
			cluScoreMichel /= pfphits.size();
		}

		for (auto const & hit : pfphits) 
		{

			if (hit->WireID().Plane != plane)       { continue; }
			if (hit->WireID().TPC != tpc)           { continue; }
			if (hit->WireID().Cryostat != cryostat) { continue; }

			const unsigned int wire = hit->WireID().Wire;
			const unsigned int time = hit->PeakTime();

			const int wireDiff = static_cast<int> (centralWire) - 
			                     static_cast<int> (wire);
			const int timeDiff = static_cast<int> (centralTick) - 
			                     static_cast<int> (time);

			if (abs(wireDiff) > static_cast<int>(nWireBins) / 2 || 
			    abs(timeDiff) > static_cast<int>(nTimeBins) / 2) 
			{ 
				continue; 
			}

			const int wireBin = (static_cast<int>(nWireBins) / 2) - wireDiff;
			const int timeBin = (static_cast<int>(nTimeBins) / 2) - timeDiff;

			HitData hitData         = GetHitData(hit, hitcnnscores);
			hitData.CLUScore_EM     = cluScoreEM;
			hitData.CLUScore_Michel = cluScoreMichel;

			const bool hitIsMichel  = HitIsMichel(hit, duneUtils);

			dataHist      -> Fill(wireBin, timeBin, hitData.Integral * norm);
			energyHist    -> Fill(wireBin, timeBin, hitData.Energy * norm);
			cnnEMHist     -> Fill(wireBin, timeBin, hitData.CNNScore_EM);
			cnnMichelHist -> Fill(wireBin, timeBin, hitData.CNNScore_Michel);
			cluEMHist     -> Fill(wireBin, timeBin, hitData.CLUScore_EM);
			cluMichelHist -> Fill(wireBin, timeBin, hitData.CLUScore_Michel);
			truthHist     -> Fill(wireBin, timeBin, static_cast<double>(hitIsMichel));

		}

	}

	std::vector<TH2D> trainingData;
	trainingData.push_back(* dataHist);
	trainingData.push_back(* energyHist);
	trainingData.push_back(* cnnEMHist);
	trainingData.push_back(* cnnMichelHist);
	trainingData.push_back(* cluEMHist);
	trainingData.push_back(* cluMichelHist);
	trainingData.push_back(* truthHist);

	return trainingData;

}

std::vector<TH2D> MichelAnalysis::DrawMichelTrainingData(
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
       const MichelAnalysis::DUNEUtils & duneUtils)
{


	// Get details of one of the hits and use it as the centre of the image
	// try to use plane 2 if possible since that is the collection view
	unsigned int centralWire     { 0 };
	unsigned int centralTick     { 0 };
	unsigned int nCollectionHits { 0 };

	art::Ptr<recob::Hit>  hit;
	for (auto const & daughterHit : daughterHits)
	{

		if (daughterHit->WireID().Plane != 2) { continue; }

		if (hit.isNull()) { hit = daughterHit; }

		nCollectionHits += 1;
		centralWire     += daughterHit->WireID().Wire;
		centralTick     += daughterHit->PeakTime();

	}

	if (nCollectionHits < 1) 
	{ 
		std::cout << "No collection plane hits\n";
		std::vector<TH2D> trainingData;

		TH2D dataHist      = TH2D();
		TH2D energyHist    = TH2D();
		TH2D cnnEMHist     = TH2D();
		TH2D cnnMichelHist = TH2D();
		TH2D cluEMHist     = TH2D();
		TH2D cluMichelHist = TH2D();
		TH2D truthHist     = TH2D();

		trainingData.push_back(dataHist);
		trainingData.push_back(energyHist);
		trainingData.push_back(cnnEMHist);
		trainingData.push_back(cnnMichelHist);
		trainingData.push_back(cluEMHist);
		trainingData.push_back(cluMichelHist);
		trainingData.push_back(truthHist);

		return trainingData;

	}

	centralWire /= nCollectionHits;
	centralTick /= nCollectionHits;

	const unsigned int plane    { hit->WireID().Plane };
	const unsigned int tpc      { hit->WireID().TPC };
	const unsigned int cryostat { hit->WireID().Cryostat };

	double driftVelocity  { duneUtils.detprop->DriftVelocity() };
	double clockFrequency { duneUtils.detclock->TPCClock().Frequency() };
	double wireBinSize    { duneUtils.geom->WirePitch(plane) };
	double timeBinSize    { driftVelocity / clockFrequency };
	double aspectRatio    { wireBinSize / timeBinSize };
	
	// TODO: read the nBins from the fhicl file
	const size_t nTimeBins = static_cast<size_t>(
	  static_cast<double>(nWireBins) * aspectRatio);


	TH2D * dataHist      = tfs->make<TH2D>((basename + "_wire").c_str(), 
	                                     "Wire Time Image", nWireBins, 0, 
	                                     nWireBins, nWireBins, 0, nTimeBins);
	TH2D * energyHist    = tfs->make<TH2D>((basename + "_energy").c_str(), 
	                                     "Energy Time Image", nWireBins, 0, 
	                                     nWireBins, nWireBins, 0, nTimeBins);
	TH2D * cnnEMHist     = tfs->make<TH2D>((basename + "_cnnem").c_str(), 
	                                       "EM Score Image", nWireBins, 0, 
	                                       nWireBins, nWireBins, 0, nTimeBins);
	TH2D * cnnMichelHist = tfs->make<TH2D>((basename + "_cnnmichel").c_str(), 
	                                     "Michel Score Image", nWireBins, 0, 
	                                     nWireBins, nWireBins, 0, nTimeBins);
	TH2D * cluEMHist     = tfs->make<TH2D>((basename + "_cluem").c_str(), 
	                                       "EM Score Image", nWireBins, 0, 
	                                       nWireBins, nWireBins, 0, nTimeBins);
	TH2D * cluMichelHist = tfs->make<TH2D>((basename + "_clumichel").c_str(), 
	                                     "Michel Score Image", nWireBins, 0, 
	                                     nWireBins, nWireBins, 0, nTimeBins);
	TH2D * truthHist     = tfs->make<TH2D>((basename + "_truth").c_str(),
	                                     "Truth Image", nWireBins, 0, nWireBins,
	                                     nWireBins, 0, nTimeBins);

	const art::FindManyP<recob::Hit> findHits(allClusters, event, particleLabel);
	for (auto const & cluster : * allClusters)
	{

		const vec_ptr_hit_t & clusterHits = findHits.at(cluster.ID());

		double cluScoreEM     = 0.;
		double cluScoreMichel = 0.;
		for (auto const & hit : clusterHits) 
		{
			const HitData & hitData = GetHitData(hit, hitcnnscores);
			cluScoreEM     += hitData.CNNScore_EM;
			cluScoreMichel += hitData.CNNScore_Michel;
		}

		if (clusterHits.size() != 0)
		{
			cluScoreEM     /= clusterHits.size();
			cluScoreMichel /= clusterHits.size();
		}

		for (auto const & hit : clusterHits) 
		{

			if (hit->WireID().Plane != plane)       { continue; }
			if (hit->WireID().TPC != tpc)           { continue; }
			if (hit->WireID().Cryostat != cryostat) { continue; }

			const unsigned int wire = hit->WireID().Wire;
			const unsigned int time = hit->PeakTime();

			const int wireDiff = static_cast<int> (centralWire) - 
			                     static_cast<int> (wire);
			const int timeDiff = static_cast<int> (centralTick) - 
			                     static_cast<int> (time);

			if (abs(wireDiff) > static_cast<int>(nWireBins) / 2 || 
			    abs(timeDiff) > static_cast<int>(nTimeBins) / 2) 
			{ 
				continue; 
			}

			const int wireBin = (static_cast<int>(nWireBins) / 2) - wireDiff;
			const int timeBin = (static_cast<int>(nTimeBins) / 2) - timeDiff;

			HitData hitData         = GetHitData(hit, hitcnnscores);
			hitData.CLUScore_EM     = cluScoreEM;
			hitData.CLUScore_Michel = cluScoreMichel;

			const bool hitIsMichel  = HitIsMichel(hit, duneUtils);

			dataHist      -> Fill(wireBin, timeBin, hitData.Integral * norm);
			energyHist    -> Fill(wireBin, timeBin, hitData.Energy * norm);
			cnnEMHist     -> Fill(wireBin, timeBin, hitData.CNNScore_EM);
			cnnMichelHist -> Fill(wireBin, timeBin, hitData.CNNScore_Michel);
			cluEMHist     -> Fill(wireBin, timeBin, hitData.CLUScore_EM);
			cluMichelHist -> Fill(wireBin, timeBin, hitData.CLUScore_Michel);
			truthHist     -> Fill(wireBin, timeBin, static_cast<double>(hitIsMichel));

		}

	}

	std::vector<TH2D> trainingData;
	trainingData.push_back(* dataHist);
	trainingData.push_back(* energyHist);
	trainingData.push_back(* cnnEMHist);
	trainingData.push_back(* cnnMichelHist);
	trainingData.push_back(* cluEMHist);
	trainingData.push_back(* cluMichelHist);
	trainingData.push_back(* truthHist);

	return trainingData;

}

std::vector<TH2D> MichelAnalysis::DrawMichelTrainingData(
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
       const MichelAnalysis::DUNEUtils & duneUtils)
{


	// Get details of one of the hits and use it as the centre of the image
	// try to use plane 2 if possible since that is the collection view
	unsigned int centralWire     { 0 };
	unsigned int centralTick     { 0 };
	unsigned int nCollectionHits { 0 };

	art::Ptr<recob::Hit>  hit;
	for (auto const & daughterHit : daughterHits)
	{

		if (daughterHit->WireID().Plane != 2) { continue; }

		if (hit.isNull()) { hit = daughterHit; }

		nCollectionHits += 1;
		centralWire     += daughterHit->WireID().Wire;
		centralTick     += daughterHit->PeakTime();

	}

	if (nCollectionHits < 1) 
	{ 
		std::cout << "No collection plane hits\n";
		std::vector<TH2D> trainingData;

		TH2D dataHist      = TH2D();
		TH2D energyHist    = TH2D();
		TH2D cnnEMHist     = TH2D();
		TH2D cnnMichelHist = TH2D();
		TH2D truthHist     = TH2D();

		trainingData.push_back(dataHist);
		trainingData.push_back(energyHist);
		trainingData.push_back(cnnEMHist);
		trainingData.push_back(cnnMichelHist);
		trainingData.push_back(truthHist);

		return trainingData;

	}

	centralWire /= nCollectionHits;
	centralTick /= nCollectionHits;

	const unsigned int plane    { hit->WireID().Plane };
	const unsigned int tpc      { hit->WireID().TPC };
	const unsigned int cryostat { hit->WireID().Cryostat };

	double driftVelocity  { duneUtils.detprop->DriftVelocity() };
	double clockFrequency { duneUtils.detclock->TPCClock().Frequency() };
	double wireBinSize    { duneUtils.geom->WirePitch(plane) };
	double timeBinSize    { driftVelocity / clockFrequency };
	double aspectRatio    { wireBinSize / timeBinSize };
	
	// TODO: read the nBins from the fhicl file
	const size_t nTimeBins = static_cast<size_t>(
	  static_cast<double>(nWireBins) * aspectRatio);


	TH2D * dataHist      = tfs->make<TH2D>((basename + "_wire").c_str(), 
	                                     "Wire Time Image", nWireBins, 0, 
	                                     nWireBins, nWireBins, 0, nTimeBins);
	TH2D * energyHist    = tfs->make<TH2D>((basename + "_energy").c_str(), 
	                                     "Energy Time Image", nWireBins, 0, 
	                                     nWireBins, nWireBins, 0, nTimeBins);
	TH2D * cnnEMHist     = tfs->make<TH2D>((basename + "_cnnem").c_str(), 
	                                       "EM Score Image", nWireBins, 0, 
	                                       nWireBins, nWireBins, 0, nTimeBins);
	TH2D * cnnMichelHist = tfs->make<TH2D>((basename + "_cnnmichel").c_str(), 
	                                     "Michel Score Image", nWireBins, 0, 
	                                     nWireBins, nWireBins, 0, nTimeBins);
	TH2D * truthHist     = tfs->make<TH2D>((basename + "_truth").c_str(),
	                                     "Truth Image", nWireBins, 0, nWireBins,
	                                     nWireBins, 0, nTimeBins);

	for (auto const & hit : allHits)
	{

		if (hit->WireID().Plane != plane)       { continue; }
		if (hit->WireID().TPC != tpc)           { continue; }
		if (hit->WireID().Cryostat != cryostat) { continue; }

		const unsigned int wire = hit->WireID().Wire;
		const unsigned int time = hit->PeakTime();

		const int wireDiff = static_cast<int> (centralWire) - 
		                     static_cast<int> (wire);
		const int timeDiff = static_cast<int> (centralTick) - 
		                     static_cast<int> (time);

		if (abs(wireDiff) > static_cast<int>(nWireBins) / 2 || 
		    abs(timeDiff) > static_cast<int>(nTimeBins) / 2) 
		{ 
			continue; 
		}

		const int wireBin = (static_cast<int>(nWireBins) / 2) - wireDiff;
		const int timeBin = (static_cast<int>(nTimeBins) / 2) - timeDiff;

		HitData hitData         = GetHitData(hit, hitcnnscores);

		const bool hitIsMichel  = HitIsMichel(hit, duneUtils);

		dataHist      -> Fill(wireBin, timeBin, hitData.Integral * norm);
		energyHist    -> Fill(wireBin, timeBin, hitData.Energy * norm);
		cnnEMHist     -> Fill(wireBin, timeBin, hitData.CNNScore_EM);
		cnnMichelHist -> Fill(wireBin, timeBin, hitData.CNNScore_Michel);
		truthHist     -> Fill(wireBin, timeBin, static_cast<double>(hitIsMichel));

	}

	std::vector<TH2D> trainingData;
	trainingData.push_back(* dataHist);
	trainingData.push_back(* energyHist);
	trainingData.push_back(* cnnEMHist);
	trainingData.push_back(* cnnMichelHist);
	trainingData.push_back(* truthHist);

	return trainingData;

}

std::vector<TH2D> MichelAnalysis::DrawMichelTrainingData(
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
       const MichelAnalysis::DUNEUtils & duneUtils)
{

	// Get details of one of the hits and use it as the centre of the image
	// try to use plane 2 if possible since that is the collection view
	unsigned int centralWire     { 0 };
	unsigned int centralTick     { 0 };
	unsigned int nCollectionHits { 0 };

	art::Ptr<recob::Hit>  hit;
	for (auto const & daughterHit : daughterHits)
	{
		if (daughterHit->WireID().Plane != 2) { continue; }
		if (hit.isNull())                     { hit = daughterHit; }

		nCollectionHits += 1;
		centralWire     += daughterHit->WireID().Wire;
		centralTick     += daughterHit->PeakTime();
	}

	if (nCollectionHits < 1) 
	{ 
		std::cout << "No collection plane hits\n";
		std::vector<TH2D> trainingData;
		return trainingData;
	}

	centralWire /= nCollectionHits;
	centralTick /= nCollectionHits;

	const unsigned int plane    { hit->WireID().Plane };
	const unsigned int tpc      { hit->WireID().TPC };
	const unsigned int cryostat { hit->WireID().Cryostat };

	double driftVelocity  { duneUtils.detprop->DriftVelocity() };
	double clockFrequency { duneUtils.detclock->TPCClock().Frequency() };
	double wireBinSize    { duneUtils.geom->WirePitch(plane) };
	double timeBinSize    { driftVelocity / clockFrequency };

	double aspectRatio    { wireBinSize / timeBinSize };
	const size_t nTimeBins = static_cast<size_t>(
	  static_cast<double>(nWireBins) * aspectRatio);

	TH2D * dataHist       = tfs->make<TH2D>((basename + "_wire").c_str(), 
	                                        "Wire Time Image", nWireBins, 0, 
	                                        nWireBins, nWireBins, 0, nTimeBins);
	TH2D * energyHist     = tfs->make<TH2D>((basename + "_energy").c_str(), 
	                                        "Energy Time Image", nWireBins, 0, 
	                                        nWireBins, nWireBins, 0, nTimeBins);
	TH2D * cnnEMHist      = tfs->make<TH2D>((basename + "_cnnem").c_str(), 
	                                        "EM Score Image", nWireBins, 0, 
	                                        nWireBins, nWireBins, 0, nTimeBins);
	TH2D * cnnMichelHist  = tfs->make<TH2D>((basename + "_cnnmichel").c_str(), 
	                                        "Michel Score Image", nWireBins, 0, 
	                                        nWireBins, nWireBins, 0, nTimeBins);
	TH2D * cluEMHist      = tfs->make<TH2D>((basename + "_cluem").c_str(), 
	                                       "EM Score Image", nWireBins, 0, 
	                                       nWireBins, nWireBins, 0, nTimeBins);
	TH2D * cluMichelHist  = tfs->make<TH2D>((basename + "_clumichel").c_str(), 
	                                        "Michel Score Image", nWireBins, 0, 
	                                        nWireBins, nWireBins, 0, nTimeBins);
	TH2D * trueEnergyHist = tfs->make<TH2D>((basename + "_trueEnergy").c_str(),
	                                        "True Energy Deposits", nWireBins, 0, 
	                                        nWireBins, nWireBins, 0, nTimeBins);
	TH2D * truthHist      = tfs->make<TH2D>((basename + "_truth").c_str(),
	                                        "Truth Image", nWireBins, 0, 
	                                        nWireBins, nWireBins, 0, nTimeBins);
	TH2D * nTrueHist      = tfs->make<TH2D>((basename + "_nTrue").c_str(),
	                                        "NumberOfTrueHits", nWireBins, 0, 
	                                        nWireBins, nWireBins, 0, nTimeBins);

	// Keep track of used hits to not double fill later
	// - We want to use clusters to get more consistent cnn scores
	// - We don't quite get all hits with clusters so need to keep track
	//   for when we get the last hits later
	std::unordered_set<size_t> usedHitKeys;
	int nTrueHits {0};
	for (auto const & pfparticle : pfparticles)
	{

		const vec_ptr_hit_t pfphits = 
		  MichelAnalysis::GetPFParticleHits(pfparticle, event, pdUtils, 
		                                    particleLabel);

		double cluScoreEM     = 0.;
		double cluScoreMichel = 0.;
		for (auto const & hit : pfphits) 
		{
			usedHitKeys.insert(hit.key());
			HitData hitData = GetHitData(hit, t0, caloAlg, duneUtils,
			                                     hitcnnscores);
			cluScoreEM     += hitData.CNNScore_EM;
			cluScoreMichel += hitData.CNNScore_Michel;
		}

		if (pfphits.size() != 0)
		{
			cluScoreEM     /= pfphits.size();
			cluScoreMichel /= pfphits.size();
		}

		for (auto const & hit : pfphits) 
		{

			if (hit->WireID().Plane    != plane)    { continue; }
			if (hit->WireID().TPC      != tpc)      { continue; }
			if (hit->WireID().Cryostat != cryostat) { continue; }

			const simb::MCParticle * mcp = HitToMCP(hit, duneUtils);
			if (mcp != nullptr) 
			{
				if (mcp -> TrackId() == mcParticle.TrackId()) { nTrueHits++; }
			}

			const unsigned int wire = hit->WireID().Wire;
			const unsigned int time = hit->PeakTime();

			const int wireDiff = static_cast<int> (centralWire) - 
			                     static_cast<int> (wire);
			const int timeDiff = static_cast<int> (centralTick) - 
			                     static_cast<int> (time);

			if (abs(wireDiff) > static_cast<int>(nWireBins) / 2 || 
			    abs(timeDiff) > static_cast<int>(nTimeBins) / 2) 
			{ 
				continue; 
			}

			const int wireBin = (static_cast<int>(nWireBins) / 2) - wireDiff;
			const int timeBin = (static_cast<int>(nTimeBins) / 2) - timeDiff;

			HitData hitData = GetHitData(hit, t0, caloAlg, duneUtils,
			                             hitcnnscores);

			hitData.CLUScore_EM     = cluScoreEM;
			hitData.CLUScore_Michel = cluScoreMichel;

			const double trueEnergy = HitToTrueEnergy(hit, duneUtils);
			const bool hitIsMichel  = HitIsMichel(hit, duneUtils);

			dataHist       -> Fill(wireBin, timeBin, hitData.Integral * norm);
			energyHist     -> Fill(wireBin, timeBin, hitData.Energy * norm);
			cnnEMHist      -> Fill(wireBin, timeBin, hitData.CNNScore_EM);
			cnnMichelHist  -> Fill(wireBin, timeBin, hitData.CNNScore_Michel);
			cluEMHist      -> Fill(wireBin, timeBin, hitData.CLUScore_EM);
			cluMichelHist  -> Fill(wireBin, timeBin, hitData.CLUScore_Michel);
			trueEnergyHist -> Fill(wireBin, timeBin, trueEnergy);
			if (hitIsMichel) 
			{
				if (mcp != nullptr) 
				{
					if (mcp -> TrackId() == mcParticle.TrackId()) 
					{ 
						truthHist -> Fill(wireBin, timeBin, 1.); 
					}
				}
			}
		}

	}

	for (auto const & hit : allHits)
	{

		std::unordered_set<size_t>::const_iterator getHit = 
		  usedHitKeys.find(hit.key());

		if (getHit != usedHitKeys.end())        { continue; }
		if (hit->WireID().Plane != plane)       { continue; }
		if (hit->WireID().TPC != tpc)           { continue; }
		if (hit->WireID().Cryostat != cryostat) { continue; }

		const simb::MCParticle * mcp = HitToMCP(hit, duneUtils);
		if (mcp != nullptr) 
		{
			if (mcp -> TrackId() == mcParticle.TrackId()) { nTrueHits++; }
		}

		const unsigned int wire = hit -> WireID().Wire;
		const unsigned int time = hit -> PeakTime();

		const int wireDiff = static_cast<int> (centralWire) - 
		                     static_cast<int> (wire);
		const int timeDiff = static_cast<int> (centralTick) - 
		                     static_cast<int> (time);

		if (abs(wireDiff) > static_cast<int>(nWireBins) / 2 || 
		    abs(timeDiff) > static_cast<int>(nTimeBins) / 2) 
		{ 
			continue; 
		}

		const int wireBin = (static_cast<int>(nWireBins) / 2) - wireDiff;
		const int timeBin = (static_cast<int>(nTimeBins) / 2) - timeDiff;

		HitData hitData = GetHitData(hit, t0, caloAlg, duneUtils,
		                             hitcnnscores);

		const double trueEnergy = HitToTrueEnergy(hit, duneUtils);
		const bool hitIsMichel  = HitIsMichel(hit, duneUtils);

		dataHist       -> Fill(wireBin, timeBin, hitData.Integral * norm);
		energyHist     -> Fill(wireBin, timeBin, hitData.Energy * norm);
		cnnEMHist      -> Fill(wireBin, timeBin, hitData.CNNScore_EM);
		cnnMichelHist  -> Fill(wireBin, timeBin, hitData.CNNScore_Michel);
		cluEMHist      -> Fill(wireBin, timeBin, hitData.CNNScore_EM);
		cluMichelHist  -> Fill(wireBin, timeBin, hitData.CNNScore_Michel);
		trueEnergyHist -> Fill(wireBin, timeBin, trueEnergy);
		if (hitIsMichel) 
		{ 
			if (mcp != nullptr) 
			{
				if (mcp -> TrackId() == mcParticle.TrackId()) 
				{ 
					truthHist -> Fill(wireBin, timeBin, 1.); 
				}
			}
		} 

	}

	for (size_t i = 0; i < (nWireBins + 2)*(nTimeBins + 2); i++)
	{
		nTrueHist -> SetBinContent(i, nTrueHits);
	}

	std::vector<TH2D> trainingData;
	trainingData.push_back(* dataHist);
	trainingData.push_back(* energyHist);
	trainingData.push_back(* cnnEMHist);
	trainingData.push_back(* cnnMichelHist);
	trainingData.push_back(* cluEMHist);
	trainingData.push_back(* cluMichelHist);
	trainingData.push_back(* trueEnergyHist);
	trainingData.push_back(* truthHist);
	trainingData.push_back(* nTrueHist);

	return trainingData;

}

std::vector<TH2D> MichelAnalysis::DrawMichelTrainingData_Data(
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
       const MichelAnalysis::DUNEUtils & duneUtils)
{

	// Get details of one of the hits and use it as the centre of the image
	// try to use plane 2 if possible since that is the collection view
	unsigned int centralWire     { 0 };
	unsigned int centralTick     { 0 };
	unsigned int nCollectionHits { 0 };

	art::Ptr<recob::Hit>  hit;
	for (auto const & daughterHit : daughterHits)
	{

		if (daughterHit->WireID().Plane != 2) { continue; }

		if (hit.isNull()) { hit = daughterHit; }

		nCollectionHits += 1;
		centralWire     += daughterHit->WireID().Wire;
		centralTick     += daughterHit->PeakTime();

	}

	if (nCollectionHits < 1) 
	{ 
		std::cout << "No collection plane hits\n";
		std::vector<TH2D> trainingData;

		TH2D dataHist       = TH2D();
		TH2D energyHist     = TH2D();
		TH2D cnnEMHist      = TH2D();
		TH2D cnnMichelHist  = TH2D();
		TH2D cluEMHist      = TH2D();
		TH2D cluMichelHist  = TH2D();

		trainingData.push_back(dataHist);
		trainingData.push_back(energyHist);
		trainingData.push_back(cnnEMHist);
		trainingData.push_back(cnnMichelHist);
		trainingData.push_back(cluEMHist);
		trainingData.push_back(cluMichelHist);

		return trainingData;

	}

	centralWire /= nCollectionHits;
	centralTick /= nCollectionHits;

	const unsigned int plane    { hit->WireID().Plane };
	const unsigned int tpc      { hit->WireID().TPC };
	const unsigned int cryostat { hit->WireID().Cryostat };

	double driftVelocity  { duneUtils.detprop->DriftVelocity() };
	double clockFrequency { duneUtils.detclock->TPCClock().Frequency() };
	double wireBinSize    { duneUtils.geom->WirePitch(plane) };
	double timeBinSize    { driftVelocity / clockFrequency };
	double aspectRatio    { wireBinSize / timeBinSize };
	
	// TODO: read the nBins from the fhicl file
	const size_t nTimeBins = static_cast<size_t>(
	  static_cast<double>(nWireBins) * aspectRatio);


	TH2D * dataHist      = tfs->make<TH2D>((basename + "_wire").c_str(), 
	                                     "Wire Time Image", nWireBins, 0, 
	                                     nWireBins, nWireBins, 0, nTimeBins);
	TH2D * energyHist    = tfs->make<TH2D>((basename + "_energy").c_str(), 
	                                     "Energy Time Image", nWireBins, 0, 
	                                     nWireBins, nWireBins, 0, nTimeBins);
	TH2D * cnnEMHist     = tfs->make<TH2D>((basename + "_cnnem").c_str(), 
	                                       "EM Score Image", nWireBins, 0, 
	                                       nWireBins, nWireBins, 0, nTimeBins);
	TH2D * cnnMichelHist = tfs->make<TH2D>((basename + "_cnnmichel").c_str(), 
	                                     "Michel Score Image", nWireBins, 0, 
	                                     nWireBins, nWireBins, 0, nTimeBins);
	TH2D * cluEMHist     = tfs->make<TH2D>((basename + "_cluem").c_str(), 
	                                       "EM Score Image", nWireBins, 0, 
	                                       nWireBins, nWireBins, 0, nTimeBins);
	TH2D * cluMichelHist = tfs->make<TH2D>((basename + "_clumichel").c_str(), 
	                                     "Michel Score Image", nWireBins, 0, 
	                                     nWireBins, nWireBins, 0, nTimeBins);

	// Keep track of used hits to not double fill later
	// - We want to use clusters to get more consistent cnn scores
	// - We don't quite get all hits with clusters so need to keep track
	//   for when we get the last hits later
	std::unordered_set<size_t> usedHitKeys;
	for (auto const & pfparticle : pfparticles)
	{

		const vec_ptr_hit_t pfphits = 
		  MichelAnalysis::GetPFParticleHits(pfparticle, event, pdUtils, 
		                                    particleLabel);

		double cluScoreEM     = 0.;
		double cluScoreMichel = 0.;
		for (auto const & hit : pfphits) 
		{
			usedHitKeys.insert(hit.key());
			HitData hitData = GetHitData(hit, t0, caloAlg, duneUtils, hitcnnscores);
			cluScoreEM     += hitData.CNNScore_EM;
			cluScoreMichel += hitData.CNNScore_Michel;
		}

		if (pfphits.size() != 0)
		{
			cluScoreEM     /= pfphits.size();
			cluScoreMichel /= pfphits.size();
		}

		for (auto const & hit : pfphits) 
		{

			if (hit->WireID().Plane != plane)       { continue; }
			if (hit->WireID().TPC != tpc)           { continue; }
			if (hit->WireID().Cryostat != cryostat) { continue; }

			const unsigned int wire = hit->WireID().Wire;
			const unsigned int time = hit->PeakTime();

			const int wireDiff = static_cast<int> (centralWire) - 
			                     static_cast<int> (wire);
			const int timeDiff = static_cast<int> (centralTick) - 
			                     static_cast<int> (time);

			if (abs(wireDiff) > static_cast<int>(nWireBins) / 2 || 
			    abs(timeDiff) > static_cast<int>(nTimeBins) / 2) 
			{ 
				continue; 
			}

			const int wireBin = (static_cast<int>(nWireBins) / 2) - wireDiff;
			const int timeBin = (static_cast<int>(nTimeBins) / 2) - timeDiff;

			HitData hitData         = GetHitData(hit, t0, caloAlg, duneUtils, hitcnnscores);
			hitData.CLUScore_EM     = cluScoreEM;
			hitData.CLUScore_Michel = cluScoreMichel;

			dataHist       -> Fill(wireBin, timeBin, hitData.Integral * norm);
			energyHist     -> Fill(wireBin, timeBin, hitData.Energy * norm);
			cnnEMHist      -> Fill(wireBin, timeBin, hitData.CNNScore_EM);
			cnnMichelHist  -> Fill(wireBin, timeBin, hitData.CNNScore_Michel);
			cluEMHist      -> Fill(wireBin, timeBin, hitData.CLUScore_EM);
			cluMichelHist  -> Fill(wireBin, timeBin, hitData.CLUScore_Michel);

		}

	}

	for (auto const & hit : allHits)
	{

		std::unordered_set<size_t>::const_iterator getHit = 
		  usedHitKeys.find(hit.key());

		if (getHit != usedHitKeys.end())        { continue; }
		if (hit->WireID().Plane != plane)       { continue; }
		if (hit->WireID().TPC != tpc)           { continue; }
		if (hit->WireID().Cryostat != cryostat) { continue; }

		const unsigned int wire = hit -> WireID().Wire;
		const unsigned int time = hit -> PeakTime();

		const int wireDiff = static_cast<int> (centralWire) - 
		                     static_cast<int> (wire);
		const int timeDiff = static_cast<int> (centralTick) - 
		                     static_cast<int> (time);

		if (abs(wireDiff) > static_cast<int>(nWireBins) / 2 || 
		    abs(timeDiff) > static_cast<int>(nTimeBins) / 2) 
		{ 
			continue; 
		}

		const int wireBin = (static_cast<int>(nWireBins) / 2) - wireDiff;
		const int timeBin = (static_cast<int>(nTimeBins) / 2) - timeDiff;

		HitData hitData = GetHitData(hit, t0, caloAlg, duneUtils, hitcnnscores);

		dataHist       -> Fill(wireBin, timeBin, hitData.Integral * norm);
		energyHist     -> Fill(wireBin, timeBin, hitData.Energy * norm);
		cnnEMHist      -> Fill(wireBin, timeBin, hitData.CNNScore_EM);
		cnnMichelHist  -> Fill(wireBin, timeBin, hitData.CNNScore_Michel);
		cluEMHist      -> Fill(wireBin, timeBin, hitData.CNNScore_EM);
		cluMichelHist  -> Fill(wireBin, timeBin, hitData.CNNScore_Michel);

	}

	std::vector<TH2D> trainingData;
	trainingData.push_back(* dataHist);
	trainingData.push_back(* energyHist);
	trainingData.push_back(* cnnEMHist);
	trainingData.push_back(* cnnMichelHist);
	trainingData.push_back(* cluEMHist);
	trainingData.push_back(* cluMichelHist);

	return trainingData;

}
