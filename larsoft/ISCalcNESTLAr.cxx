////////////////////////////////////////////////////////////////////////
// Class:       ISCalcNESTLAr
// Plugin Type: Algorithm
// File:        ISCalcNESTLAr.cxx
// Description:
// Aug. 30 by Mu Wei
//  Edited by N. Carrara to use LArNEST (05/29/2023)
////////////////////////////////////////////////////////////////////////

#include "larsim/IonizationScintillation/ISCalcNESTLAr.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include <algorithm>

namespace {
  constexpr double LAr_Z{18};
  constexpr double Density_LAr{1.393};

  constexpr double scint_yield{1.0 / (19.5 * CLHEP::eV)};
  constexpr double resolution_scale{0.107}; // Doke 1976
}

namespace larg4 {

  //----------------------------------------------------------------------------
  ISCalcNESTLAr::ISCalcNESTLAr(CLHEP::HepRandomEngine& Engine)
    : fEngine(Engine), fSCE{lar::providerFrom<spacecharge::SpaceChargeService>()}
  {
    LArDetector* detector = new LArDetector();
    mLArNEST = new NEST::LArNEST(detector);

    std::cout << "ISCalcNESTLAr Initialize." << std::endl;
  }

  //----------------------------------------------------------------------------
  ISCalcData ISCalcNESTLAr::CalcIonAndScint(detinfo::DetectorPropertiesData const& detProp,
                                            sim::SimEnergyDeposit const& edep)
  {
    
    std::cout << "HERE IN NEST!" << std::endl;
    return {0.0, 0.0, 0.0, 0.0};

    // return {energyDeposit,
    //         static_cast<double>(NumElectrons),
    //         static_cast<double>(NumPhotons),
    //         GetScintYieldRatio(edep)};
  }

  

  //----------------------------------------------------------------------------
  double ISCalcNESTLAr::EFieldAtStep(double efield, sim::SimEnergyDeposit const& edep)
  {
    geo::Point_t pos = edep.MidPoint();
    double EField = efield;
    geo::Vector_t eFieldOffsets;
    if (fSCE->EnableSimEfieldSCE()) {
      eFieldOffsets = fSCE->GetEfieldOffsets(pos);
      EField =
        std::sqrt((efield + efield * eFieldOffsets.X()) * (efield + efield * eFieldOffsets.X()) +
                  (efield * eFieldOffsets.Y() * efield * eFieldOffsets.Y()) +
                  (efield * eFieldOffsets.Z() * efield * eFieldOffsets.Z()));
    }
    return EField;
  }

}
