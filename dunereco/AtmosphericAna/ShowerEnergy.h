
#ifndef DUNE_NEUTRINO_SHOWER_RECO_ALG_H
#define DUNE_NEUTRINO_SHOWER_RECO_ALG_H

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
#include "dunereco/FDSensOpt/NeutrinoEnergyRecoAlg/NeutrinoEnergyRecoAlg.h"
#include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"


namespace atm{
    
    class ShowerEnergy
    {
        public:
            ShowerEnergy(const fhicl::ParameterSet &pset,std::string showerlabel, std::string tracklabel, std::string pfplabel, std::string wirelabel)
            double ShowerEnergyAtm(const art::Ptr<recob::Shower> &pElectronShower, const art::Event &event)

        private:

            double ElectronShowerEnergy;
            dune::NeutrinoEnergyRecoAlg fNeutrinoEnergyRecoAlg;
    };
}
#endif