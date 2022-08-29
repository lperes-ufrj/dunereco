/**
*  @file   dunereco/AtmosphericAna/Atmospheric.h
*
*  @brief  Header file for the neutrino atmospheric algorithm. 
*
*  $Log: $
*/
#ifndef DUNE_NEUTRINO_ATMOSPHERIC_ALG_H
#define DUNE_NEUTRINO_ATMOSPHERIC_ALG_H


// Art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <fstream>
#include <array>
#include <iterator>

// some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "Math/GenVector/LorentzVector.h"
#include "TVector3.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

// "art" includes (canvas, and gallery)
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
//#include "gallery/Event.h"
//#include "gallery/ValidHandle.h"

//#include "fhiclcpp/ParameterSet.h"

// LArSoft, nutools includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/Shower.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/MVAPIDResult.h"
//#include "larana/ParticleIdentification/Chi2PIDAlg.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/TrackingTypes.h"

// Geometry includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/TestUtils/geometry_unit_test_base.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// DUNE libraries
#include "dunecore/Geometry/DuneApaChannelMapAlg.h"
#include "dunecore/Geometry/GeoObjectSorterAPA.h"
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"
#include "dunereco/FDSensOpt/NeutrinoEnergyRecoAlg/NeutrinoEnergyRecoAlg.h"
#include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"


namespace atm
{
  class Atmospheric;
}
class atm::Atmospheric : public art::EDAnalyzer
{
public:
  typedef art::Handle<std::vector<recob::PFParticle>> PFParticleHandle;
  typedef art::Handle<std::vector<recob::Track>> TrackHandle;
  typedef art::Handle<std::vector<recob::Shower>> ShowerHandle;
  typedef art::Handle<std::vector<recob::SpacePoint>> SpacePointHandle;
  typedef art::Handle<std::vector<recob::Hit>> HitHandle;
  typedef std::map<size_t, art::Ptr<recob::PFParticle>> PFParticleIdMap;

  typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<Coord_t>, ROOT::Math::GlobalCoordinateSystemTag> vtxPos;


  explicit Atmospheric(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Atmospheric(Atmospheric const &) = delete;
  Atmospheric(Atmospheric &&) = delete;
  Atmospheric &operator=(Atmospheric const &) = delete;
  Atmospheric &operator=(Atmospheric &&) = delete;

  // Required functions.
  void analyze(art::Event const &e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:
  // Declare member data here.

  // Root Tree
  TTree *m_AtmTree; ///<
  TTree *m_AllEvents;

  // Implemented functions
  bool IsVisibleParticle(int PdgCode, std::string process);
  void ResetCounters();
  void GeoLimits(art::ServiceHandle<geo::Geometry const> &geom, float fFidVolCutX, float fFidVolCutY, float fFidVolCutZ);
  bool insideFV(geo::Point_t const &vertex);
  double CalcPIDA( std::vector<art::Ptr<anab::Calorimetry>> calos, bool IsFlipped);
  double PIDAFunc( double dedx, double resrng );
  void IsTrackBackwardsAndFlip( std::vector<art::Ptr<anab::Calorimetry>> calos, TVector3 trkStart, TVector3 trkEnd, bool IsFlipped);

  //TVector3 RandomSunPosition;
  //std::ifstream SunPositions;

  // Variables to save in the tree
  int m_run;   ///<
  int m_event; ///<
  int fnGeantParticles;

  std::vector<int> fNPrimaryDaughters;
  int fNPrimaries;

  std::vector<int> fDaughterTrackID;

  std::vector<float> fTrackStartX, fTrackStartY, fTrackStartZ;
  std::vector<float> fTrackEndX, fTrackEndY, fTrackEndZ;
  std::vector<bool> fIsTrackFlipped;

  std::vector<float> fminChi2value;
  std::vector<int> fminChi2PDG;

  bool InvertTrack;

  std::vector<int> fPrimaryPDGReco;
  std::vector<std::vector<double>> fPrimaryRecoVertex;
  std::vector<double> fTotalMomRecoRangeUnitVect;
  std::vector<double> fTotalMomRecoCalVectUnit;
  double fMCCosAzimuthNu;
  std::vector<int> fShowerID;
  std::vector<float> fShowerDirectionX;
  std::vector<float> fShowerDirectionY;
  std::vector<float> fShowerDirectionZ;
  std::vector<std::vector<double>> fShowerDirectionErr;
  std::vector<std::vector<double>> fShowerShowerStart;
  std::vector<std::vector<double>> fShowerShowerStartErr;
  std::vector<int> fShowerbest_plane;
  std::vector<double> fShowerLength;
  std::vector<double> fShowerOpenAngle;
  std::vector<std::vector<double>> fShowerEnergy;
  std::vector<double> fPIDANoFlip;
  std::vector<double> fPIDAwithFlipMaybe;

  float fDistVertex = -1;

  //Reco BDT Variables
  double fCosThetaDetTotalMom;
  double fCosPhiDetTotalMom;
  double fTotalMomentumP;
  int fnTracks;
  int fnShowers;
  int fnSpacePoints;
  double fPIDALongestTrack;
  double fLongestTrack;
  double fHighestTrackSummedADC;
  double fHighestShowerSummedADC;
  double fLargeShowerOpenAngle;
  double fLongestShower;
  float fCVN_NCScore;
  double fFracTotalChargeLongTrack;
  double fAvarageTrackLength;
  float fEventRecoEnergy;

  bool fIsNC_CVNPred;

  //Check Direction 

  double fDiffCosAngleTotalMom;
  double fDiffCosAngleLongestTrack;


  //MVA bits - Just one for event -> Takes the longest track
  double fRecoTrackMVAEvalRatio;
  double fRecoTrackMVAConcentration;
  double fRecoTrackMVACoreHaloRatio;
  double fRecoTrackMVAConicalness;
  double fRecoTrackMVAdEdxStart;
  double fRecoTrackMVAdEdxEnd;
  double fRecoTrackMVAdEdxEndRatio;
  double fRecoTrackMVAElectron;
  double fRecoTrackMVAPion;
  double fRecoTrackMVAMuon;
  double fRecoTrackMVAProton;
  double fRecoTrackMVAPhoton;


  //BackTracking
  std::vector<int> fDaughterTrackTruePDG;
  std::vector<int> fRecoTrackTruePDG;
  std::vector<bool> fRecoTrackTruePrimary;
  std::vector<double> fRecoTrackTrueMomX;
  std::vector<double> fRecoTrackTrueMomY;
  std::vector<double> fRecoTrackTrueMomZ;
  std::vector<double> fRecoTrackTrueMomT;
  std::vector<double> fRecoTrackTrueStartX;
  std::vector<double> fRecoTrackTrueStartY;
  std::vector<double> fRecoTrackTrueStartZ;
  std::vector<double> fRecoTrackTrueStartT;
  std::vector<double> fRecoTrackTrueEndX;
  std::vector<double> fRecoTrackTrueEndY;
  std::vector<double> fRecoTrackTrueEndZ;
  std::vector<double> fRecoTrackTrueEndT;
 
  //Truth information
  std::vector<int> fCCNC;
  std::vector<double> fThetaNuLepton;
  std::vector<int> fMCPrimaryNuPDG;
  std::vector<std::vector<double>> fMCInitialPositionNu;
  std::vector<std::vector<double>> fMCNuMomentum;
  std::vector<int> fMCTrackId;
  std::vector<int> fMCPdgCode;
  std::vector<int> fMCMother;
  std::vector<int> fMCNumberDaughters;
  std::vector<std::string> fMCProcess;
  std::vector<std::string> fMCEndProcess;
  std::vector<std::vector<double>> fMCPosition;
  std::vector<std::vector<double>> fMCEndPosition;
  std::vector<std::vector<double>> fMCMomentum;
  std::vector<std::vector<double>> fMCEndMomentum;
  std::vector<double> fMCStartEnergy;
  std::vector<double> fMCEndEnergy;
  std::vector<int> fMCStatusCode;
  std::vector<bool> fisMCinside;
  std::vector<double> fTotalMomentumUnitVect;
  std::vector<double> fEnergyShowerLinearlyCorrected;

  // Module labels and parameters
  std::string fGeneratorModuleLabel;
  std::string fGeantModuleLabel;
  std::string fShowerModuleLabel;
  std::string fTrackModuleLabel;
  std::string fPFParticleModuleLabel;
  std::string fCaloModuleLabel;
  std::string fPIDModuleLabel;
  std::string fSpacePointModuleLabel;
  std::string fHitModuleLabel;
  std::string fCVNModuleLabel;
  std::string fMVAPIDModuleLabel;
  bool fSaveGeantInfo;
  bool fCheatVertex;
  bool fShowerRecoSave;
  std::string fWireModuleLabel;

  trkf::TrackMomentumCalculator trkm;
  dune::NeutrinoEnergyRecoAlg fNeutrinoEnergyRecoAlg;
  calo::CalorimetryAlg fCalorimetryAlg;                    ///< the calorimetry algorithm

};
#endif //DUNE_NEUTRINO_ATMOSPHERIC_ALG_H