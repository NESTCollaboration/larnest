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
#include "larsim/IonizationScintillation/ISTPC.h"
#include "larsim/IonizationScintillation/LArNEST.h"

#include "CLHEP/Random/RandBinomial.h"

namespace spacecharge {
  class SpaceCharge;
}

namespace CLHEP {
  class HepRandomEngine;
}

namespace larg4 {
  class ISCalcNESTLAr : public ISCalc {
  public:
    explicit ISCalcNESTLAr(detinfo::DetectorPropertiesData const& detProp, CLHEP::HepRandomEngine& fEngine);

    double EFieldAtStep(double efield,
                        sim::SimEnergyDeposit const& edep)
      override; //value of field with any corrections for this step
    double AngleToEFieldAtStep(double efield, sim::SimEnergyDeposit const& edep);
    ISCalcData CalcIonAndScint(detinfo::DetectorPropertiesData const& detProp,
                               sim::SimEnergyDeposit const& edep) override;

  private:
    CLHEP::HepRandomEngine& fEngine; // random engine
    ISTPC fISTPC;
    CLHEP::RandBinomial fBinomialGen;
    const spacecharge::SpaceCharge* fSCE;
    
    larnest::LArNEST mLArNEST;

    double fGeVToElectrons;      ///< from LArG4Parameters service
    double fWion;                ///< W_ion (23.6 eV) == 1/fGeVToElectrons
    double fWph;                 ///< from LArG4Parameters service
    double fScintPreScale;       ///< scintillation pre-scaling factor from LArProperties service
    double fRecombA;             ///< from LArG4Parameters service
    double fRecombk;             ///< from LArG4Parameters service
    double fModBoxA;             ///< from LArG4Parameters service
    double fModBoxB;             ///< from LArG4Parameters service
    double fEllipsModBoxA;       ///< from LArG4Parameters service
    double fEllipsModBoxB;       ///< from LArG4Parameters service
    double fEllipsModBoxR;       ///< from LArG4Parameters service
    double fLarqlChi0A;          ///< from LArG4Parameters service
    double fLarqlChi0B;          ///< from LArG4Parameters service
    double fLarqlChi0C;          ///< from LArG4Parameters service
    double fLarqlChi0D;          ///< from LArG4Parameters service
    double fLarqlAlpha;          ///< from LArG4Parameters service
    double fLarqlBeta;           ///< from LArG4Parameters service
    double fQAlpha;              ///< from LArG4Parameters service
    bool fUseModBoxRecomb;       ///< from LArG4Parameters service
    bool fUseEllipsModBoxRecomb; ///< from LArG4Parameters service
    bool fUseModLarqlRecomb;     ///< from LArG4Parameters service
    bool fUseBinomialFlucts;     ///< from LArG4Parameters service

    double EscapingEFraction(
      double const dEdx) const; //LArQL chi0 function = fraction of escaping electrons
    double FieldCorrection(double const EF, double const dEdx)
      const; //LArQL f_corr function = correction factor for electric field dependence
  };
  };
}
#endif // LARG4_ISCALCNESTLAr_H
