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

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

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

namespace larg4 
{
    class ISCalcNESTLAr : public ISCalc 
    {
    public:
        explicit ISCalcNESTLAr(
            fhicl::ParameterSet const& pset, 
            detinfo::DetectorPropertiesData const& detProp, 
            CLHEP::HepRandomEngine& fEngine
        );

        double EFieldAtStep(
            double efield,
            sim::SimEnergyDeposit const& edep
        ) override; //value of field with any corrections for this step

        double AngleToEFieldAtStep(double efield, sim::SimEnergyDeposit const& edep);

        ISCalcData CalcIonAndScint(
            detinfo::DetectorPropertiesData const& detProp,
            sim::SimEnergyDeposit const& edep
        ) override;

    private:
        fhicl::ParameterSet fParameterSet;
        CLHEP::HepRandomEngine& fEngine; // random engine
        ISTPC fISTPC;
        CLHEP::RandBinomial fBinomialGen;
        const spacecharge::SpaceCharge* fSCE;
        
        larnest::LArNEST mLArNEST;

        art::InputTag fNESTMode;
    };
}
#endif // LARG4_ISCALCNESTLAr_H
