////////////////////////////////////////////////////////////////////////
// Class:       ISCalcNESTLAr
// Plugin Type: Algorithm
// File:        ISCalcNESTLAr.cxx
// Description:
// Aug. 30 by Mu Wei
//  Edited by N. Carrara to use LArNEST (05/29/2023)
////////////////////////////////////////////////////////////////////////

#ifndef LARG4_ISCALCNESTLAr_H
#define LARG4_ISCALCNESTLAr_H

#include "larsim/IonizationScintillation/ISCalc.h"
#include "larsim/larnest/LArDetector.hh"
#include "larsim/larnest/LArNEST.hh"

namespace spacecharge {
  class SpaceCharge;
}

namespace CLHEP {
  class HepRandomEngine;
}

namespace larg4 {
  class ISCalcNESTLAr : public ISCalc {
  public:
    explicit ISCalcNESTLAr(CLHEP::HepRandomEngine& fEngine);

    double EFieldAtStep(double efield,
                        sim::SimEnergyDeposit const& edep)
      override; //value of field with any corrections for this step
    ISCalcData CalcIonAndScint(detinfo::DetectorPropertiesData const& detProp,
                               sim::SimEnergyDeposit const& edep) override;

  private:
    CLHEP::HepRandomEngine& fEngine; // random engine
    const spacecharge::SpaceCharge* fSCE;
    
    NEST::LArNEST* mLArNEST;
  };
}
#endif // LARG4_ISCALCNESTLAr_H