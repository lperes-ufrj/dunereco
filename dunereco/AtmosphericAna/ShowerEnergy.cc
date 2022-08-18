#include "dunereco/AtmosphericAna/ShowerEnergy.h"

namespace atm
{


ShowerEnergy::ShowerEnergy(fhicl::ParameterSet const &p,
      fShowerModuleLabel(p.get<std::string>("ShowerModuleLabel")),
      fTrackModuleLabel(p.get<std::string>("TrackModuleLabel")),
      fPFParticleModuleLabel(p.get<std::string>("PFParticleModuleLabel")),
      fWireModuleLabel(p.get<std::string>("WireModuleLabel")),
      fNeutrinoEnergyRecoAlg(p.get<fhicl::ParameterSet>("NeutrinoEnergyRecoAlg"),fTrackModuleLabel,fShowerModuleLabel,
        fHitModuleLabel,fWireModuleLabel,fTrackModuleLabel,fShowerModuleLabel,fPFParticleModuleLabel)
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

double atm::ShowerEnergy(const art::Ptr<recob::Shower> &pElectronShower, const art::Event &event){

    ElectronShowerEnergy = fNeutrinoEnergyRecoAlg.CalculateElectronEnergy(pElectronShower,event);
  fEventRecoEnergy = energyRecoHandle->fNuLorentzVector.E();

    return ElectronShowerEnergy;
}



}