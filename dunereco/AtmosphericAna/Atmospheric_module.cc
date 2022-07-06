
////////////////////////////////////////////////////////////////////////
// Class:       Atmospheric
// Plugin Type: analyzer (Unknown Unknown)
// File:        Atmospheric_module.cc
//
// Generated at Mon Jan 17 19:47:24 2022 by Leonardo Peres using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

//Some PIDA and flip functions are adaptted from ProtonIdentification_module.cc 

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



// Detector Limits =========================
float fFidVolXmin = 0;
float fFidVolXmax = 0;
float fFidVolYmin = 0;
float fFidVolYmax = 0;
float fFidVolZmin = 0;
float fFidVolZmax = 0;

const int kMaxHits = 250;

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

  //std::vector<std::vector<float>> fDaughterTrackdEdx_Plane0;
  //std::vector<std::vector<float>> fDaughterTrackdEdx_Plane1;
  //std::vector<std::vector<float>> fDaughterTrackdEdx_Plane2;
  //std::vector<std::vector<float>> fDaughterTrackResidualRange_Plane0;
  //std::vector<std::vector<float>> fDaughterTrackResidualRange_Plane1;
  //std::vector<std::vector<float>> fDaughterTrackResidualRange_Plane2;
  //std::vector<float> fDaughterKE_Plane0;
  //std::vector<float> fDaughterKE_Plane1;
  //std::vector<float> fDaughterKE_Plane2;

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
  int    CaloPlane0, CaloPlane1, CaloPlane2;
  double ResRPlane0[kMaxHits], ResRPlane1[kMaxHits], ResRPlane2[kMaxHits];
  double dEdxPlane0[kMaxHits], dEdxPlane1[kMaxHits], dEdxPlane2[kMaxHits];

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
  float fCVN_NCProbability;
  double fFracTotalChargeLongTrack;
  double fAvarageTrackLength;

  //BackTracking
  std::vector<int> fDaughterTrackTruePDG;
 
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
  bool fSaveGeantInfo;
  bool fCheatVertex;
  bool fShowerRecoSave;

  trkf::TrackMomentumCalculator trkm;

};

atm::Atmospheric::Atmospheric(fhicl::ParameterSet const &p)
    : EDAnalyzer{p}, // More initializers here.
      fGeneratorModuleLabel(p.get<std::string>("GeneratorModuleLabel")),
      fGeantModuleLabel(p.get<std::string>("GeantModuleLabel")),
      fShowerModuleLabel(p.get<std::string>("ShowerModuleLabel")),
      fTrackModuleLabel(p.get<std::string>("TrackModuleLabel")),
      fPFParticleModuleLabel(p.get<std::string>("PFParticleModuleLabel")),
      fCaloModuleLabel(p.get<std::string>("CaloModuleLabel")),
      fPIDModuleLabel(p.get<std::string>("PIDModuleLabel")),
      fSpacePointModuleLabel(p.get<std::string>("SpacePointModuleLabel")),
      fHitModuleLabel(p.get<std::string>("HitModuleLabel")),
      fCVNModuleLabel(p.get<std::string>("CVNModuleLabel")),
      fSaveGeantInfo(p.get<bool>("SaveGeantInfo")),
      fCheatVertex(p.get<bool>("CheatVertex")),
      fShowerRecoSave(p.get<bool>("ShowerRecoSave"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void atm::Atmospheric::ResetCounters()
{
 // std::cout << "Reseting counters... " << std::endl;

  fCosThetaDetTotalMom = -2;
  fCosPhiDetTotalMom = -2;
  fnTracks = 0;
  fnShowers = 0;
  fnSpacePoints = 0;
  fTotalMomentumP = 0;
  fPIDALongestTrack = 0;
  fHighestTrackSummedADC = 0;
  fLongestTrack = 0;
  fHighestShowerSummedADC = 0;
  fLargeShowerOpenAngle = 0;
  fLongestShower = 0;
  fCVN_NCProbability = -1;
  fFracTotalChargeLongTrack = 0;
  fAvarageTrackLength = 0;


  fDistVertex = -2;


  fCCNC.clear();
  fNPrimaryDaughters.clear();
  fPrimaryPDGReco.clear();
  fNPrimaries = 0;
  InvertTrack = false;

  fDaughterTrackID.clear();
  fTrackStartX.clear();
  fTrackStartY.clear();
  fTrackStartZ.clear();
  fTrackEndX.clear();
  fTrackEndY.clear();
  fTrackEndZ.clear();
  fIsTrackFlipped.clear();
  fPIDAwithFlipMaybe.clear();
  fPIDANoFlip.clear();
  fPrimaryRecoVertex.clear();
  

  //fDaughterTrackdEdx_Plane0.clear();
  //fDaughterTrackdEdx_Plane1.clear();
  //fDaughterTrackdEdx_Plane2.clear();
  //fDaughterTrackResidualRange_Plane0.clear();
  //fDaughterTrackResidualRange_Plane1.clear();
  //fDaughterTrackResidualRange_Plane2.clear();  

  fminChi2PDG.clear();
  fminChi2value.clear();
  fTotalMomRecoRangeUnitVect.clear();
  fTotalMomRecoCalVectUnit.clear();
  fDaughterTrackTruePDG.clear();

  fnGeantParticles = -1;
  fThetaNuLepton.clear();
  fMCNuMomentum.clear();
  fMCPrimaryNuPDG.clear();
  fMCInitialPositionNu.clear();
  fMCCosAzimuthNu = -2;
  fMCTrackId.clear();
  fMCPdgCode.clear();
  fMCMother.clear();
  fMCProcess.clear();
  fMCEndProcess.clear();
  fMCNumberDaughters.clear();
  fMCPosition.clear();
  fMCEndPosition.clear();
  fMCMomentum.clear();
  fMCStartEnergy.clear();
  fMCEndMomentum.clear();
  fMCEndEnergy.clear();
  fMCStatusCode.clear();
  fisMCinside.clear();
  fTotalMomentumUnitVect.clear();

  fShowerID.clear();
  fShowerDirectionX.clear();
  fShowerDirectionY.clear();
  fShowerDirectionZ.clear();
  fShowerDirectionErr.clear();
  fShowerShowerStart.clear();
  fShowerShowerStartErr.clear();
  fShowerbest_plane.clear();
  fShowerLength.clear();
  fShowerOpenAngle.clear();
  fShowerEnergy.clear();
}

void atm::Atmospheric::analyze(art::Event const &evt)
{
  // pid::Chi2PIDAlg fChiAlg;
  //std::cout << "Test 1 " << std::endl;
  ResetCounters();
  //std::cout << "Test 2 " << std::endl;
  TVector3 TotalMomentumRecoRange;
  TVector3 TotalMomentumRecoCal;
  TLorentzVector TotalMomentumTrue;
  // Get Geometry
  art::ServiceHandle<geo::Geometry const> geom;
  GeoLimits(geom, 10, 10, 10);

  //we need the particle inventory service
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

  // and the clock data for event
  // auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  // channel quality
  // lariov::ChannelStatusProvider const& channelStatus = art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider();

  size_t nplanes = geom->Nplanes();
  std::vector<std::vector<unsigned int>> hits(nplanes);

  m_run = evt.run();
  m_event = evt.id().event();

 // std::cout << "  Run: " << m_run << std::endl;
 // std::cout << "  Event: " << m_event << std::endl;

    // Get CVN results
  art::Handle<std::vector<cvn::Result>> cvnResult;
  evt.getByLabel(fCVNModuleLabel, "cvnresult", cvnResult); //not sure if i need an instance name?? 

  if(!cvnResult->empty()) 
  {
      fCVN_NCProbability = (*cvnResult)[0].GetNCProbability();
  }


  std::map<int,float> PDGtoMass;
  PDGtoMass.insert(std::pair<int,float>(2212, 0.938272));
  PDGtoMass.insert(std::pair<int,float>(211, 0.13957));
  PDGtoMass.insert(std::pair<int,float>(321, 0.493677));
  PDGtoMass.insert(std::pair<int,float>(13, 0.105658));

  TVector3 vertical(0,1,0);

  auto const &MCTruthHandle = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorModuleLabel);
  auto const &MCTruthObjs = *MCTruthHandle;


  for (size_t i = 0; i < MCTruthObjs.size(); i++)
  {
    simb::MCTruth MCTruthObj = MCTruthObjs[i];
    fnGeantParticles = MCTruthObj.NParticles();

    const simb::MCNeutrino &MCNeutrino = MCTruthObj.GetNeutrino();
    fCCNC.push_back(MCNeutrino.CCNC());
    fThetaNuLepton.push_back(MCNeutrino.Theta());
    const simb::MCParticle &nu = MCNeutrino.Nu();
    fMCPrimaryNuPDG.push_back(nu.PdgCode());
    std::vector<double> tmp_NuMomentum = {nu.Momentum().X(), nu.Momentum().Y(), nu.Momentum().Z()};
    fMCNuMomentum.push_back(tmp_NuMomentum);
    std::vector<double> tmp_fMCInitialPositionNu = {nu.Vx(), nu.Vy(), nu.Vz()};
    fMCInitialPositionNu.push_back(tmp_fMCInitialPositionNu);
    
    fMCCosAzimuthNu =  vertical*nu.Momentum().Vect().Unit() ;
  
  }

  //std::cout << "Number of CC and/or NC interactions: " << fCCNC.size() << std::endl;
  //std::cout << "Nu Interaction (0=CC 1=NC): " << fCCNC[0] << std::endl;

  // Truth information (Save geant4 info)
  if (fSaveGeantInfo)
  {

    auto const &MCParticleHandle = evt.getValidHandle<std::vector<simb::MCParticle>>(fGeantModuleLabel);
    auto const &MCParticleObjs = *MCParticleHandle;

    fnGeantParticles = MCParticleObjs.size();
      for (size_t i = 0; i < MCParticleObjs.size(); i++)
      {
      

        const simb::MCParticle MCParticleObj = MCParticleObjs[i];
        fMCTrackId.push_back(MCParticleObj.TrackId());
        fMCPdgCode.push_back(MCParticleObj.PdgCode());
        fMCMother.push_back(MCParticleObj.Mother());
        fMCProcess.push_back(MCParticleObj.Process());
        fMCEndProcess.push_back(MCParticleObj.EndProcess());
        fMCNumberDaughters.push_back(MCParticleObj.NumberDaughters());
        fMCPosition.push_back({MCParticleObj.Position().X(), MCParticleObj.Position().Y(), MCParticleObj.Position().Z()});
        fMCEndPosition.push_back({MCParticleObj.EndPosition().X(), MCParticleObj.EndPosition().Y(), MCParticleObj.EndPosition().Z()});
        fMCMomentum.push_back({MCParticleObj.Px(), MCParticleObj.Py(), MCParticleObj.Pz()});
        fMCStartEnergy.push_back(MCParticleObj.E(0));
        fMCEndMomentum.push_back({MCParticleObj.EndPx(), MCParticleObj.EndPy(), MCParticleObj.EndPz()});
        fMCEndEnergy.push_back(MCParticleObj.EndE());
        fMCStatusCode.push_back(MCParticleObj.StatusCode());
        geo::Point_t EndMCPos(MCParticleObj.EndPosition().X(), MCParticleObj.EndPosition().Y(), MCParticleObj.EndPosition().Z());
        bool isMCinside_tmp = insideFV(EndMCPos);
        fisMCinside.push_back(isMCinside_tmp);

        // True total momentum -> it has to be a visible particle,
        //               it has to be stopped inside the detector,
        //               it has to be a primary particle, 
        //               it has to be a stable final state particle.

        if (IsVisibleParticle(MCParticleObj.PdgCode(), MCParticleObj.Process()) && isMCinside_tmp && MCParticleObj.StatusCode() == 1)
        {
          TotalMomentumTrue += MCParticleObj.Momentum();
        }
      }
    
    TVector3 TotalMomTrue = TotalMomentumTrue.Vect().Unit();
    std::vector<double> TotalMomTrueXYZ = {TotalMomTrue.X(), TotalMomTrue.Y(), TotalMomTrue.Z()};
    fTotalMomentumUnitVect = TotalMomTrueXYZ;
  }

    m_AllEvents->Fill();

  // Collect the PFParticles from the event
  PFParticleHandle pfParticleHandle;
  evt.getByLabel(fPFParticleModuleLabel, pfParticleHandle);

  TrackHandle trackHandle;
  evt.getByLabel(fTrackModuleLabel, trackHandle);

  ShowerHandle showerHandle;
  evt.getByLabel(fShowerModuleLabel, showerHandle);

  SpacePointHandle spacepointHandle;
  std::vector< art::Ptr<recob::SpacePoint> > SpacePointVect;
  if(evt.getByLabel(fPFParticleModuleLabel,spacepointHandle)) art::fill_ptr_vector(SpacePointVect, spacepointHandle);

  HitHandle hitHandle;
  std::vector< art::Ptr<recob::Hit> > HitVect;
  if(evt.getByLabel(fHitModuleLabel,hitHandle)) art::fill_ptr_vector(HitVect, hitHandle);



  if (!pfParticleHandle.isValid())
  {
    mf::LogDebug("Atmospheric") << "  Failed to find the PFParticles." << std::endl;
    return;
  }

  // Get the associations between PFParticles and tracks/showers from the event and track from calo
  art::FindManyP<recob::Track> pfPartToTrackAssoc(pfParticleHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Shower> pfPartToShowerAssoc(pfParticleHandle, evt, fShowerModuleLabel);
  art::FindManyP<anab::Calorimetry> TrackToCaloAssoc(trackHandle, evt, fCaloModuleLabel);
  art::FindManyP<recob::Hit> TrackToHitsAssoc(trackHandle, evt, fTrackModuleLabel);
  art::FindManyP<anab::ParticleID> TrackToPIDAssoc(trackHandle, evt, fPIDModuleLabel);
  art::FindManyP<recob::SpacePoint> pfpToSpacePoint(pfParticleHandle, evt, fSpacePointModuleLabel);
  art::FindManyP<recob::Vertex> pfPartToVertex(pfParticleHandle, evt,fPFParticleModuleLabel);
  art::FindManyP<recob::Hit> ShowerToHitAssoc(showerHandle, evt, fShowerModuleLabel);
  //art::FindManyP<recob::Hit> TrackToHitsAssoc(trackHandle, evt, fHitModuleLabel);


  const std::vector<art::Ptr<recob::PFParticle>> pfparticleVect = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt, fPFParticleModuleLabel);

  size_t neutrinoID = 99999;
 // Double_t VertexXYZ[3] = {};

  // Block adapted from ConsolidatedPFParticleAnalysisTemplate should work fine!
    for (size_t iPfp = 0; iPfp < pfparticleVect.size(); iPfp++)
    {

      //const std::vector<art::Ptr<recob::SpacePoint>> &associatedSpacePoints = pfpToSpacePoint.at(iPfp);
      const std::vector<art::Ptr<recob::Vertex>> &associatedVertex = pfPartToVertex.at(iPfp);
      
      if(!(pfparticleVect[iPfp]->IsPrimary() && (pfparticleVect[iPfp]->PdgCode() == 14 || pfparticleVect[iPfp]->PdgCode() == 12 )) ) continue;

    
      // const int pdg(pParticle->PdgCode());
      // std::cout << "We have Primary Particles! Yep" << std::endl;
      //std::cout << "associatedVertex.size() = " << associatedVertex.size() << std::endl;

      for (const art::Ptr<recob::Vertex> &vtx : associatedVertex)
      {
        if(!(vtx->isValid())) continue;
        vtxPos vertex = vtx->position();
       // vtx->XYZ(VertexXYZ);

        fPrimaryRecoVertex.push_back({vertex.X(),vertex.Y(),vertex.Z()});
        fDistVertex = pow(pow((fMCInitialPositionNu.at(0).at(0)-vertex.X()),2)+pow((fMCInitialPositionNu.at(0).at(1)-vertex.Y()),2)+pow((fMCInitialPositionNu.at(0).at(2)-vertex.Z()),2),0.5);
      }

      neutrinoID = pfparticleVect[iPfp]->Self();
      fNPrimaryDaughters.push_back(pfparticleVect[iPfp]->NumDaughters());
      fPrimaryPDGReco.push_back(pfparticleVect[iPfp]->PdgCode());
      
      fNPrimaries++;
    }

    if (neutrinoID == 99999) return;

    fnSpacePoints+=SpacePointVect.size();

    double TotalHitsADC = 0;

    for (size_t iHit = 0; iHit < HitVect.size(); iHit++){

      TotalHitsADC += HitVect[iHit]->SummedADC();

    }


    double LongestTrack = 1e-10;
    double HighestTrackSummedADC = 1e-10;
    long unsigned int LongestTrackID = -2;
    double PIDALongestTrack = 0;
    double AllLengthTracksSummed = 0;

  for (size_t iPfp = 0; iPfp < pfparticleVect.size(); iPfp++){

      if(pfparticleVect[iPfp]->Parent() != neutrinoID) continue;
      const std::vector<art::Ptr<recob::Track>> &associatedTracks = pfPartToTrackAssoc.at(iPfp);

      float minChi2value;
      int minChi2PDG;
      fnTracks+=associatedTracks.size();

      if (!associatedTracks.empty())
      {

  
        // std::cout << "We have Tracks! Yep" << std::endl;
        for (const art::Ptr<recob::Track> &trk : associatedTracks)
        {
          std::vector<art::Ptr<recob::Hit>> trackHits = TrackToHitsAssoc.at(trk.key());

          if(trk->NumberTrajectoryPoints() < 2) continue; 

          if(!insideFV(trk->End())) continue;
         // int trkidtruth = TruthMatchUtils::TrueParticleIDFromTotalTrueEnergy(clockData, trackHits, true);
         // const simb::MCParticle *particle = pi_serv->TrackIdToParticle_P(trkidtruth);

          fDaughterTrackID.push_back(trk.key());

         // fDaughterTrackTruePDG.push_back(particle->PdgCode());
          float trackADC = 0;

          for(const art::Ptr<recob::Hit> &hit : trackHits){
            
            trackADC += hit->SummedADC();
          }
          if(trackADC > fHighestTrackSummedADC) HighestTrackSummedADC = trackADC;
          
          AllLengthTracksSummed += trk->Length();
          if(trk->Length() > LongestTrack){
              LongestTrack = trk->Length();
              LongestTrackID = trk.key();
          } 
          
          std::vector<art::Ptr<anab::ParticleID>> trackPID = TrackToPIDAssoc.at(trk.key());
          std::map<int,float> PDGtoChi2;

          for (size_t i = 0; i < trackPID.size(); i++)
          {

            std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(i)->ParticleIDAlgScores();
            minChi2value = 9999;
            minChi2PDG = 9999;

            // Loop through AlgScoresVec and find the variables we want
            for (size_t i_algscore = 0; i_algscore < AlgScoresVec.size(); i_algscore++)
            {
              anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
              // if (AlgScore.fPlaneMask[2] != 1)  continue; //Only collection plane
              if(AlgScore.fAlgName == "Chi2") {
                //Sum Chi2 hyphotesis all planes 
                PDGtoChi2[AlgScore.fAssumedPdg] = AlgScore.fValue;


              }
            }
          }

          //std::cout << "\nThe map PDGtoChi2 is : \n";
          //std::cout << "\tKEY\tELEMENT\n";
          std::map<int, float>::iterator it;
          for( it = PDGtoChi2.begin(); it != PDGtoChi2.end(); ++it){
           // std::cout << '\t' << it->first << '\t' << it->second << '\n';
            if(it->second < minChi2value){
              minChi2value = it->second;
              minChi2PDG = it->first;
            }
          }
         // std::cout << "min Ch2 hypothesis : " << minChi2PDG << " with " << minChi2value << ". \n";
              fminChi2value.push_back(minChi2value);
              fminChi2PDG.push_back(minChi2PDG);
          // std::min(chi2_pion_backward, std::min(chi2_pion_forward, std::min(chi2_proton_backward, chi2_proton_forward)));

          const std::vector<art::Ptr<anab::Calorimetry>> associatedCalo = TrackToCaloAssoc.at(trk.key());
          if (associatedCalo.empty())
            continue;
          float KE = 0;
          
          
          for (const art::Ptr<anab::Calorimetry> &cal : associatedCalo)
          {
            if (!cal->PlaneID().isValid)
              continue;
            int planenum = cal->PlaneID().Plane;
            // std::cout << "pid: " << pidout.ParticleIDAlgScores.at(0) << std::endl;
            std::vector<float> temp_dEdx = cal->dEdx();
           // std::cout<< "temp_dEdx.size() = "<< temp_dEdx.size() << std::endl;
            //if (InvertTrack)
            //{
            //  //std::cout << "temp_dEdx.at(0)= " << temp_dEdx.at(0) << std::endl;
            //  std::reverse(temp_dEdx.begin(), temp_dEdx.end());
            //  //std::cout << "temp_dEdx.at(final)= " << temp_dEdx.at(temp_dEdx.size() - 1) << std::endl;
            //}

            //if (planenum == 0)
            //{
            //  fDaughterTrackdEdx_Plane0.push_back(temp_dEdx);
            //  fDaughterTrackResidualRange_Plane0.push_back(cal->ResidualRange());
            //  fDaughterKE_Plane0.push_back(cal->KineticEnergy());
            //}
//
            //if (planenum == 1)
            //{
            //  fDaughterTrackdEdx_Plane1.push_back(temp_dEdx);
            //  fDaughterTrackResidualRange_Plane1.push_back(cal->ResidualRange());
            //  fDaughterKE_Plane1.push_back(cal->KineticEnergy());
            //}

            if (planenum == 2)
            {
              //fDaughterTrackdEdx_Plane2.push_back(temp_dEdx);
              //fDaughterTrackResidualRange_Plane2.push_back(cal->ResidualRange());
              //fDaughterKE_Plane2.push_back(cal->KineticEnergy());
              KE = cal->KineticEnergy();
            }

            temp_dEdx.clear();
          } // Calo

        TVector3 TrackDirection(trk->StartDirection().X(), trk->StartDirection().Y(), trk->StartDirection().Z());

        TVector3 DaughterStartPoint(trk->Start().X(), trk->Start().Y(), trk->Start().Z());
        TVector3 DaughterEndPoint(trk->End().X(), trk->End().Y(), trk->End().Z());

        fTrackStartX.push_back(DaughterStartPoint.X());fTrackStartY.push_back(DaughterStartPoint.Y()); fTrackStartZ.push_back(DaughterStartPoint.Z());
        fTrackEndX.push_back(DaughterEndPoint.X()); fTrackEndY.push_back(DaughterEndPoint.Y()); fTrackEndZ.push_back(DaughterEndPoint.Z());

          // Cheat Vertex, all track directions based on the true vertex
        double PIDANoFlip = CalcPIDA(associatedCalo, InvertTrack);
        fPIDANoFlip.push_back(PIDANoFlip); 

        //IsTrackBackwardsAndFlip(associatedCalo, DaughterStartPoint, DaughterEndPoint, InvertTrack);

        if (InvertTrack) TrackDirection = -TrackDirection;
   
        double PIDAwithFlipMaybe = CalcPIDA(associatedCalo, InvertTrack);
        //std::cout << "InvertTrack = " << InvertTrack << "." << std::endl;
         //  std::cout << "PIDA = " << PIDA << std::endl;
        fPIDAwithFlipMaybe.push_back(PIDAwithFlipMaybe);
        fIsTrackFlipped.push_back(InvertTrack);

        if(trk.key() == LongestTrackID) PIDALongestTrack = PIDAwithFlipMaybe;

        if(minChi2PDG == 13 || minChi2PDG == 211 ){
          TotalMomentumRecoRange += trkm.GetTrackMomentum(trk->Length(), 13) * TrackDirection;
        //  std::cout << "TotalMomentumRecoRange.Mag() = " << TotalMomentumRecoRange.Mag() << std::endl;
        }

        if(minChi2PDG == 2212 || minChi2PDG == 321 ){
          TotalMomentumRecoRange += trkm.GetTrackMomentum(trk->Length(), 2212) * TrackDirection;
        //  std::cout << "TotalMomentumRecoRange.Mag() = " << TotalMomentumRecoRange.Mag() << std::endl;
        }

          double Pcal = sqrt((KE + PDGtoMass[minChi2PDG])*(KE + PDGtoMass[minChi2PDG])-PDGtoMass[minChi2PDG] * PDGtoMass[minChi2PDG]);
          TotalMomentumRecoCal += Pcal * TrackDirection;

      } // Tracks- loop

      fAvarageTrackLength =  AllLengthTracksSummed/fnTracks;
      
    }   // if there is a track
  }

    for (size_t iPfp = 0; iPfp < pfparticleVect.size(); iPfp++)
    {
      if (fShowerRecoSave)
      {

        if(pfparticleVect[iPfp]->Parent() != neutrinoID) continue;
        // SHOWERS RECO INFO ====================================================================
        const std::vector<art::Ptr<recob::Shower>> &associatedShowers = pfPartToShowerAssoc.at(iPfp);
        fnShowers += associatedShowers.size();

        if (!associatedShowers.empty())
        {
          for (const art::Ptr<recob::Shower> &Shower : associatedShowers)
          {
            fShowerID.push_back(Shower->ID());
            fShowerDirectionX.push_back(Shower->Direction().X());
            fShowerDirectionY.push_back(Shower->Direction().Y());
            fShowerDirectionZ.push_back(Shower->Direction().Z());
            fShowerDirectionErr.push_back({Shower->DirectionErr().X(), Shower->DirectionErr().Y(), Shower->DirectionErr().Z()});
            fShowerShowerStart.push_back({Shower->ShowerStart().X(), Shower->ShowerStart().Y(), Shower->ShowerStart().Z()});
            fShowerShowerStartErr.push_back({Shower->ShowerStartErr().X(), Shower->ShowerStartErr().Y(), Shower->ShowerStartErr().Z()});
            fShowerbest_plane.push_back(Shower->best_plane());
            fShowerEnergy.push_back(Shower->Energy());
            fShowerLength.push_back(Shower->Length());
            fShowerOpenAngle.push_back(Shower->OpenAngle()); 
            
            std::vector<art::Ptr<recob::Hit>> ShowerHits = ShowerToHitAssoc.at(Shower.key());
            
            if (Shower->Length() > fLongestShower) fLongestShower = Shower->Length();
            if (Shower->OpenAngle() > fLargeShowerOpenAngle) fLargeShowerOpenAngle = Shower->OpenAngle();

            
            float showerADC = 0;

            for(const art::Ptr<recob::Hit> &hit : ShowerHits){
            
              showerADC+= hit->SummedADC();

            }
            if(showerADC > fHighestShowerSummedADC) fHighestShowerSummedADC = showerADC;


          }
        }
      }
    }

    //TVector3 z(0,0,1);
    //std::cout << "TotalMomentumRecoRange.Mag() = " << TotalMomentumRecoRange.Mag() << std::endl; 
    if (TotalMomentumRecoRange.Mag() > 0.0)
    {
    
    fHighestTrackSummedADC = HighestTrackSummedADC;
    fLongestTrack = LongestTrack;
    fPIDALongestTrack = PIDALongestTrack;
    fFracTotalChargeLongTrack = HighestTrackSummedADC/TotalHitsADC;

    fTotalMomentumP = TotalMomentumRecoRange.Mag();
    fCosThetaDetTotalMom = TotalMomentumRecoRange.Unit().CosTheta();
    fCosPhiDetTotalMom = cos(TotalMomentumRecoRange.Unit().Phi());
    fTotalMomRecoRangeUnitVect = {TotalMomentumRecoRange.Unit().X(), TotalMomentumRecoRange.Unit().Y(), TotalMomentumRecoRange.Unit().Z()};
    //std::cout << "fTotalMomRecoRangeUnitVect =" << fTotalMomRecoRangeUnitVect <<std::endl;
    fTotalMomRecoCalVectUnit = {TotalMomentumRecoCal.Unit().X(), TotalMomentumRecoCal.Unit().Y(), TotalMomentumRecoCal.Unit().Z()}; 
    //std::cout << "fTotalMomRecoCalVectUnit =" << fTotalMomRecoCalVectUnit <<std::endl;
    }
  

  //Fill the Tree just for a NC event, and if there is at least one track per event in the BDT variables
  if(TotalMomentumRecoRange.Mag() > 0) m_AtmTree->Fill();


  
}

//Evaluate the charge deposition in the beggining and at the end of the track
void atm::Atmospheric::IsTrackBackwardsAndFlip( std::vector<art::Ptr<anab::Calorimetry>> calos, TVector3 trkStart, TVector3 trkEnd, bool IsFlipped ){

  // I need to decide which way to go through the hits...
  double SumdEdx_St = 0., SumdEdx_En = 0.; //ResRng_St = 0, ResRng_En = 0;
  unsigned int AvHits = 5;
 // if (calos[2]->dEdx().size()) {
 //   ResRng_St  = calos[2]->ResidualRange()[0];
 //   ResRng_En  = calos[2]->ResidualRange()[calos[2]->dEdx().size()-1];
 // }
  // If don't have 2*AvHits (10) hits on the collection plane then use mid point + 1 hit
  if ( calos[2]->dEdx().size() < (2*AvHits) ) {
    AvHits = 1 + (0.5 * calos[2]->dEdx().size());
  }
  for ( unsigned int PlHit=0; PlHit < calos[2]->dEdx().size(); ++PlHit ) { // loop through hits on the collection plane
    if ( PlHit <= AvHits ) {
      SumdEdx_St += calos[2]->dEdx()[PlHit];
    }
    if ( calos[2]->dEdx().size() - PlHit <= AvHits ) {
      SumdEdx_En += calos[2]->dEdx()[PlHit];
    }
    //std::cout << "Looking at hit " << PlHit << " of " << (int)calos[2]->dEdx().size() << "...SumdEdx_St = " << SumdEdx_St << ", and SumdEdx_En = " << SumdEdx_En << std::endl;
  }
  double AvdEdx_St = SumdEdx_St / AvHits;
  double AvdEdx_En = SumdEdx_En / AvHits;
  // The dEdx at the start of the track should be less than that at the end...
  bool LowdEdxSt = false;
  std::cout << "AvdEdx_St = " << AvdEdx_St << "\t" <<  "AvdEdx_En = " << AvdEdx_En << "." << std::endl;
  if ( AvdEdx_St < AvdEdx_En )
    LowdEdxSt = true;
  
  if ( !LowdEdxSt) {
    std::cout << "Track backwards, energy deposition at the start higher than at the end!" << std::endl;
    double TmpX, TmpY, TmpZ;
    TmpX = trkStart.X(); TmpY = trkStart.Y(); TmpZ = trkStart.Z();
    trkStart.SetX(trkEnd.X()); trkStart.SetY(trkEnd.Y()); trkStart.SetZ(trkEnd.Z());
    trkEnd.SetX(TmpX); trkEnd.SetY(TmpY); trkEnd.SetZ(TmpZ);
  }

  InvertTrack = !LowdEdxSt;
  std::cout << "InvertTrack = " << InvertTrack << std::endl;
}

double atm::Atmospheric::CalcPIDA(  std::vector<art::Ptr<anab::Calorimetry>> calos, bool IsInverted ){

  double PIDA  = 0;
  int UsedHits = 0, TotHits = 0;
  double dEdxSum = 0;

  // *********** How do I deal with backwards tracks....What do I do if only one is backwards?.....
  //  If the dEdx is the wrong way around then without truth I would assume that the track is backwards.
  //  This means that I should use whether the MC is correct as a later cut.
  //  So if MCCorOrient == false then get PIDA using start of 'track'
  for ( int Plane=0; Plane < (int)calos.size(); ++Plane ) { // Loop through planes
    double PlanePIDA=0; int PlaneHits=0;
    for ( int PlaneHit=0; PlaneHit < (int)calos[Plane]->dEdx().size(); ++PlaneHit ) { // loop through hits on each plane
      double ThisdEdx = 0;
      double ThisResR = 0;
      if(IsInverted){
        ThisdEdx = calos[Plane]->dEdx()[PlaneHit];
        ThisResR = calos[Plane]->ResidualRange()[(int)calos[Plane]->dEdx().size()-1-PlaneHit];
      } else {
        ThisdEdx = calos[Plane]->dEdx()[PlaneHit];
        ThisResR = calos[Plane]->ResidualRange()[PlaneHit];
      }
      // Increment TotHits and dEdx sum
      dEdxSum += ThisdEdx;
      ++TotHits;
      // Fill the calorimetry objects
      if (Plane==0 && CaloPlane0<kMaxHits) {
	dEdxPlane0[CaloPlane0] = ThisdEdx;
	ResRPlane0[CaloPlane0] = ThisResR;
	++CaloPlane0;
      } else if (Plane==1 && CaloPlane1<kMaxHits) {
	dEdxPlane1[CaloPlane1] = ThisdEdx;
	ResRPlane1[CaloPlane1] = ThisResR;
	++CaloPlane1;
      } else if (Plane==2 && CaloPlane2<kMaxHits) {
	dEdxPlane2[CaloPlane2] = ThisdEdx;
	ResRPlane2[CaloPlane2] = ThisResR;
	++CaloPlane2;
      }
      // Output some calorimetry information
      /*
      double ThisRang = calos[Plane]->Range();
      double ThisKinE = calos[Plane]->KineticEnergy();
      TVector3 ThisPos = calos[Plane]->XYZ()[PlaneHit];
      std::cout << "Looking at hit " << PlaneHit << " of " << (int)calos[Plane]->dEdx().size() << "...dEdx is " << ThisdEdx << ", ResR is " 
		<< ThisResR << ", this Range is " << ThisRang << ", TrackLength is " << TrackLength << ", KinE is " << ThisKinE 
		<< ". This is pos is (" << ThisPos[0] << ", " << ThisPos[1] << ", " << ThisPos[2] << ")." 
		<< std::endl;
      //*/
      
      // ==== If MCCorOrient == true
      // Work out PIDA if ResRange < 30 cm
      if ( ThisResR < 30 ) { // Only want PIDA for last 30 cm
	PlanePIDA += PIDAFunc( ThisdEdx, ThisResR );
	++PlaneHits;
      } // If ResRange < 30 cm
      // ===== This is where I need to do things...
    } // Loop over hits on each plane
      // Increment whole track PIDA.
    PIDA     += PlanePIDA;
    UsedHits += PlaneHits;
    // Work out PIDA for this plane
    PlanePIDA = PlanePIDA/PlaneHits;
    // if (Plane == 0 ) PIDA_Plane0 = PlanePIDA;
    // else if (Plane == 1 ) PIDA_Plane1 = PlanePIDA;
    // else if (Plane == 2 ) PIDA_Plane2 = PlanePIDA;
  } // Loop over planes
  
  if ( UsedHits ) // If had any hits, work out PIDA and calculate
    PIDA = PIDA / UsedHits;
  // AvdEdx = dEdxSum / TotHits;
  return PIDA;
} // CalcPIDA

double atm::Atmospheric::PIDAFunc( double dedx, double resrng ) {
  double Va = dedx * pow( resrng, 0.42 );
  return Va;
}

void atm::Atmospheric::beginJob()
{

  std::cout << " atm::Atmospheric::beginJob() - initializing..." << std::endl;

/*

  SunPositions.open("background_sun_pos.dat");
  int nPosSun = 100000;
  int nline_SunPos = rand() % nPosSun;
  int nline = 1;
  double x = 0, y = 0, z = 0;

  while (1)
  {
     if (SunPositions.eof())
        break;
     SunPositions >> x >> y >> z;
     // if (SunFile.fail())
     // throw std::runtime_error("Can't read x, y or z at line " + std::to_string(iLine));
     //std::cout << "x: " << x << " y: " << y << " z: " << z << std::endl;
     if(nline==nline_SunPos){
      // std::cout << "nline==nline_SunPos\n";
       RandomSunPosition.SetXYZ(x, y, z);
      // std::cout << "RandomSunPosition.Mag() = " << RandomSunPosition.Mag() << std::endl;
     //  std::cout << "x,y,z = " <<x << "\t" << y <<"\t" << z << std::endl;
     }
     nline++;
  }*/

  art::ServiceHandle<art::TFileService const> tfs;
  m_AtmTree = tfs->make<TTree>("Atm", "AtmosphericAnalysis");
  m_AllEvents = tfs->make<TTree>("AllEvents", "AllEvents");

  m_AllEvents->Branch("event", &m_event, "event/I");
  m_AllEvents->Branch("CCNC", &fCCNC);

  m_AtmTree->Branch("CCNC", &fCCNC);
  m_AtmTree->Branch("run", &m_run, "run/I");
  m_AtmTree->Branch("event", &m_event, "event/I");
  m_AtmTree->Branch("LongestTrack", &fLongestTrack);

  m_AtmTree->Branch("NPrimaryDaughters", &fNPrimaryDaughters);
  m_AtmTree->Branch("NPrimaries", &fNPrimaries);
  m_AtmTree->Branch("DaughterTrackID", &fDaughterTrackID);
  m_AtmTree->Branch("CVN_NCProbability", &fCVN_NCProbability);
  

  m_AtmTree->Branch("TrackStartX", &fTrackStartX);
  m_AtmTree->Branch("TrackStartY", &fTrackStartY);
  m_AtmTree->Branch("TrackStartZ", &fTrackStartZ);
  m_AtmTree->Branch("TrackEndX", &fTrackEndX);
  m_AtmTree->Branch("TrackEndY", &fTrackEndY);
  m_AtmTree->Branch("TrackEndZ", &fTrackEndZ);

  m_AtmTree->Branch("PrimaryRecoVertex", &fPrimaryRecoVertex);
  //m_AtmTree->Branch("DaughterTrackdEdx_Plane0", &fDaughterTrackdEdx_Plane0);
  //m_AtmTree->Branch("DaughterTrackResidualRange_Plane0", &fDaughterTrackResidualRange_Plane0);
  //m_AtmTree->Branch("DaughterKE_Plane0", &fDaughterKE_Plane0);
  //m_AtmTree->Branch("DaughterTrackdEdx_Plane1", &fDaughterTrackdEdx_Plane1);
  //m_AtmTree->Branch("DaughterTrackResidualRange_Plane1", &fDaughterTrackResidualRange_Plane1);
  //m_AtmTree->Branch("DaughterKE_Plane1", &fDaughterKE_Plane1);
  //m_AtmTree->Branch("DaughterTrackdEdx_Plane2", &fDaughterTrackdEdx_Plane2);
  //m_AtmTree->Branch("DaughterTrackResidualRange_Plane2", &fDaughterTrackResidualRange_Plane2);
  //m_AtmTree->Branch("DaughterKE_Plane2", &fDaughterKE_Plane2);
  m_AtmTree->Branch("PIDA_NoFlip", &fPIDANoFlip);
  m_AtmTree->Branch("PIDAwithFlipMaybe", &fPIDAwithFlipMaybe);
  m_AtmTree->Branch("HighestTrackSummedADC", &fHighestTrackSummedADC);
  m_AtmTree->Branch("HighestShowerSummedADC", &fHighestShowerSummedADC);
  m_AtmTree->Branch("PIDALongestTrack", &fPIDALongestTrack);
  m_AtmTree->Branch("DaughterTrackTruePDG", &fDaughterTrackTruePDG);
  m_AtmTree->Branch("IsTrackFlipped", &fIsTrackFlipped);
  m_AtmTree->Branch("PrimaryPDGReco", &fPrimaryPDGReco);
  m_AtmTree->Branch("LargeShowerOpenAngle", &fLargeShowerOpenAngle);
  m_AtmTree->Branch("LongestShower", &fLongestShower);

  m_AtmTree->Branch("FracTotalChargeLongTrack", &fFracTotalChargeLongTrack);
  m_AtmTree->Branch("AvarageTrackLength", &fAvarageTrackLength);

  m_AtmTree->Branch("ShowerID", &fShowerID);
  m_AtmTree->Branch("ShowerDirectionX", &fShowerDirectionX);
  m_AtmTree->Branch("ShowerDirectionY", &fShowerDirectionY);
  m_AtmTree->Branch("ShowerDirectionZ", &fShowerDirectionZ);
  m_AtmTree->Branch("ShowerDirectionErr", &fShowerDirectionErr);
  m_AtmTree->Branch("ShowerShowerStart", &fShowerShowerStart);
  m_AtmTree->Branch("ShowerShowerStartErr", &fShowerShowerStartErr);
  m_AtmTree->Branch("Showerbest_plane", &fShowerbest_plane);
  m_AtmTree->Branch("ShowerEnergy", &fShowerEnergy);
  m_AtmTree->Branch("ShowerLength", &fShowerLength);
  m_AtmTree->Branch("ShowerOpenAngle", &fShowerOpenAngle);

  m_AtmTree->Branch("TotalMomRecoRangeUnitVect",  &fTotalMomRecoRangeUnitVect);
  m_AtmTree->Branch("TotalMomRecoCalVectUnit",  &fTotalMomRecoCalVectUnit);
  m_AtmTree->Branch("CosThetaDetTotalMom", &fCosThetaDetTotalMom);  
  m_AtmTree->Branch("CosPhiDetTotalMom", &fCosPhiDetTotalMom);
  m_AtmTree->Branch("nTracks", &fnTracks);
  m_AtmTree->Branch("nShowers", &fnShowers);
  m_AtmTree->Branch("TotalMomentumP", &fTotalMomentumP);
  m_AtmTree->Branch("nSpacePoints", &fnSpacePoints);
  m_AtmTree->Branch("minChi2value", &fminChi2value);
  m_AtmTree->Branch("minChi2PDG", &fminChi2PDG);
  m_AtmTree->Branch("DistVertex", &fDistVertex);

 // m_AtmTree->Branch("SunDirectionFromTrueBDM", &fSunDirectionFromTrueBDM);
 // m_AtmTree->Branch("TruePrimaryBDMVertex", &fPrimaryBDMVertex);
  m_AtmTree->Branch("CCNC", &fCCNC);
  m_AtmTree->Branch("ThetaNuLepton", &fThetaNuLepton);
  m_AtmTree->Branch("MCPrimaryNuPDG", &fMCPrimaryNuPDG);
  m_AtmTree->Branch("MCNuMomentum", &fMCNuMomentum);
  m_AtmTree->Branch("MCInitialPositionNu", &fMCInitialPositionNu);
  m_AtmTree->Branch("MCCosAzimuthNu", &fMCCosAzimuthNu);
  m_AtmTree->Branch("nGeantParticles", &fnGeantParticles);
  m_AtmTree->Branch("MCTrackId", &fMCTrackId);
  m_AtmTree->Branch("MCPdgCode", &fMCPdgCode);
  m_AtmTree->Branch("MCMother", &fMCMother);
  m_AtmTree->Branch("MCProcess", &fMCProcess);
  m_AtmTree->Branch("MCEndProcess", &fMCEndProcess);
  m_AtmTree->Branch("MCNumberDaughters", &fMCNumberDaughters);
  m_AtmTree->Branch("MCPosition", &fMCPosition);
  m_AtmTree->Branch("MCEndPosition", &fMCEndPosition);
  m_AtmTree->Branch("MCMomentum", &fMCMomentum);
  m_AtmTree->Branch("MCStartEnergy", &fMCStartEnergy);
  m_AtmTree->Branch("MCEndMomentum", &fMCEndMomentum);
  m_AtmTree->Branch("MCEndEnergy", &fMCEndEnergy);
  m_AtmTree->Branch("MCStatusCode", &fMCStatusCode);
  m_AtmTree->Branch("isMCinside", &fisMCinside);
  m_AtmTree->Branch("TotalMomentumUnitVect", &fTotalMomentumUnitVect);

  //  std::cout << " atm::Atmospheric::beginJob() - End" << std::endl;
}

void atm::Atmospheric::endJob()
{
  // Implementation of optional member function here.
}

bool atm::Atmospheric::IsVisibleParticle(int Pdg, std::string process)
{

  Pdg = abs(Pdg);
  bool condition = false;

  if ((Pdg == 130 || Pdg == 130 || Pdg == 211 || (Pdg > 300 && Pdg < 400) || Pdg == 2212 || (Pdg > 3000 && Pdg < 4000) || Pdg == 22 || Pdg == 13 || Pdg == 11) && process == "primary")
  {
    condition = true;
  }

  return condition;
}

void atm::Atmospheric::GeoLimits(art::ServiceHandle<geo::Geometry const> &geom, float fFidVolCutX, float fFidVolCutY, float fFidVolCutZ)
{

  // Define histogram boundaries (cm).
  // For now only draw cryostat=0.
  double minx = 1e9;
  double maxx = -1e9;
  double miny = 1e9;
  double maxy = -1e9;
  double minz = 1e9;
  double maxz = -1e9;
  for (size_t i = 0; i < geom->NTPC(); ++i)
  {
    double local[3] = {0., 0., 0.};
    double world[3] = {0., 0., 0.};
    const geo::TPCGeo &tpc = geom->TPC(i);
    tpc.LocalToWorld(local, world);
    if (minx > world[0] - geom->DetHalfWidth(i))
      minx = world[0] - geom->DetHalfWidth(i);
    if (maxx < world[0] + geom->DetHalfWidth(i))
      maxx = world[0] + geom->DetHalfWidth(i);
    if (miny > world[1] - geom->DetHalfHeight(i))
      miny = world[1] - geom->DetHalfHeight(i);
    if (maxy < world[1] + geom->DetHalfHeight(i))
      maxy = world[1] + geom->DetHalfHeight(i);
    if (minz > world[2] - geom->DetLength(i) / 2.)
      minz = world[2] - geom->DetLength(i) / 2.;
    if (maxz < world[2] + geom->DetLength(i) / 2.)
      maxz = world[2] + geom->DetLength(i) / 2.;
  }

  fFidVolXmin = minx + fFidVolCutX;
  fFidVolXmax = maxx - fFidVolCutX;
  fFidVolYmin = miny + fFidVolCutY;
  fFidVolYmax = maxy - fFidVolCutY;
  fFidVolZmin = minz + fFidVolCutZ;
  fFidVolZmax = maxz - fFidVolCutZ;
} // GeoLimits()

bool atm::Atmospheric::insideFV(geo::Point_t const &vertex)
{

  double const x = vertex.X();
  double const y = vertex.Y();
  double const z = vertex.Z();

  return x > fFidVolXmin && x < fFidVolXmax &&
         y > fFidVolYmin && y < fFidVolYmax &&
         z > fFidVolZmin && z < fFidVolZmax;

} // insideFV()

DEFINE_ART_MODULE(atm::Atmospheric)
