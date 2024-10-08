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
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larsim/Simulation/LArG4Parameters.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <algorithm>

namespace {
    constexpr double LAr_Z{18};
    constexpr double Density_LAr{1.393};

    constexpr double scint_yield{1.0 / (19.5 * CLHEP::eV)};
    constexpr double resolution_scale{0.107}; // Doke 1976
}

namespace larg4 
{

    //----------------------------------------------------------------------------
    ISCalcNESTLAr::ISCalcNESTLAr(
        fhicl::ParameterSet const& pset,
        detinfo::DetectorPropertiesData const& detProp, 
        CLHEP::HepRandomEngine& Engine
    )
    : fParameterSet{pset}
    , fEngine(Engine)
    , fISTPC{*(lar::providerFrom<geo::Geometry>())}
    , fBinomialGen{CLHEP::RandBinomial(Engine)}
    , fSCE{lar::providerFrom<spacecharge::SpaceChargeService>()}
    {
        std::cout << "ISCalcNESTLAr Initialize." << std::endl;

        art::ServiceHandle<sim::LArG4Parameters const> LArG4PropHandle;

        // Adjust NEST parameters
        fNESTMode = pset.get<art::InputTag>("NESTMode");
        std::cout << "NEST Mode: " << fNESTMode << std::endl;
    }

    //----------------------------------------------------------------------------
    ISCalcData ISCalcNESTLAr::CalcIonAndScint(
        detinfo::DetectorPropertiesData const& detProp,
        sim::SimEnergyDeposit const& edep
    )
    {
        // Get energy deposit and other parameters
        double const energy_deposit = edep.Energy();
        double dx = edep.StepLength();
        double EFieldStep = EFieldAtStep(detProp.Efield(), edep);

        // Determine species based on calculation mode and particle type
        larnest::LArInteraction species = larnest::LArInteraction::BOX;

        larnest::LArYieldResult yields = mLArNEST.GetYields(
            species, 
            energy_deposit * 1000,  // NEST expects keV 
            dx, 
            EFieldStep, 
            1.39
        );
        double NumElectrons = yields.Ne;
        double NumPhotons = yields.Nph;

        // std::cout << "HERE IN NEST! num_e: " << NumElectrons << ", num_ph: " << NumPhotons << std::endl;

        return {
            energy_deposit,
            static_cast<double>(NumElectrons),
            static_cast<double>(NumPhotons),
            0.0
        };
    }

    

    //----------------------------------------------------------------------------
    double ISCalcNESTLAr::EFieldAtStep(
        double efield, 
        sim::SimEnergyDeposit const& edep
    )
    {
        // electric field outside active volume set to zero
        if (!fISTPC.isScintInActiveVolume(edep.MidPoint())) return 0.;

        TVector3 elecvec;

        art::ServiceHandle<geo::Geometry const> fGeometry;
        geo::TPCID tpcid = fGeometry->PositionToTPCID(edep.MidPoint());
        if (!bool(tpcid)) return 0.;
        const geo::TPCGeo& tpcGeo = fGeometry->TPC(tpcid);

        if (tpcGeo.DetectDriftDirection() == 1) elecvec.SetXYZ(1, 0, 0);
        if (tpcGeo.DetectDriftDirection() == -1) elecvec.SetXYZ(-1, 0, 0);
        if (tpcGeo.DetectDriftDirection() == 2) elecvec.SetXYZ(0, 1, 0);
        if (tpcGeo.DetectDriftDirection() == -2) elecvec.SetXYZ(0, -1, 0);
        if (tpcGeo.DetectDriftDirection() == 3) elecvec.SetXYZ(0, 0, 1);
        if (tpcGeo.DetectDriftDirection() == -3) elecvec.SetXYZ(0, 0, -1);

        elecvec *= efield;

        if (fSCE->EnableSimEfieldSCE()) {
        auto const eFieldOffsets = fSCE->GetEfieldOffsets(edep.MidPoint());
        TVector3 scevec;

        if (tpcGeo.DetectDriftDirection() == 1)
            scevec.SetXYZ(
            efield * eFieldOffsets.X(), efield * eFieldOffsets.Y(), efield * eFieldOffsets.Z());
        if (tpcGeo.DetectDriftDirection() == -1)
            scevec.SetXYZ(
            -1 * efield * eFieldOffsets.X(), efield * eFieldOffsets.Y(), efield * eFieldOffsets.Z());
        if (tpcGeo.DetectDriftDirection() == 2)
            scevec.SetXYZ(
            efield * eFieldOffsets.X(), efield * eFieldOffsets.Y(), efield * eFieldOffsets.Z());
        if (tpcGeo.DetectDriftDirection() == -2)
            scevec.SetXYZ(
            efield * eFieldOffsets.X(), -1 * efield * eFieldOffsets.Y(), efield * eFieldOffsets.Z());
        if (tpcGeo.DetectDriftDirection() == 3)
            scevec.SetXYZ(
            efield * eFieldOffsets.X(), efield * eFieldOffsets.Y(), efield * eFieldOffsets.Z());
        if (tpcGeo.DetectDriftDirection() == -3)
            scevec.SetXYZ(
            efield * eFieldOffsets.X(), efield * eFieldOffsets.Y(), -1 * efield * eFieldOffsets.Z());

        elecvec += scevec;
        }

        return elecvec.Mag();
    }

    double ISCalcNESTLAr::AngleToEFieldAtStep(
        double efield, 
        sim::SimEnergyDeposit const& edep
    )
    {

        // electric field outside active volume set to zero
        if (!fISTPC.isScintInActiveVolume(edep.MidPoint())) return 0.;

        TVector3 stepvec(
        edep.StartX() - edep.EndX(), edep.StartY() - edep.EndY(), edep.StartZ() - edep.EndZ());

        TVector3 elecvec;

        art::ServiceHandle<geo::Geometry const> fGeometry;
        geo::TPCID tpcid = fGeometry->PositionToTPCID(edep.MidPoint());
        if (!bool(tpcid)) return 0.;
        const geo::TPCGeo& tpcGeo = fGeometry->TPC(tpcid);

        if (tpcGeo.DetectDriftDirection() == 1) elecvec.SetXYZ(1, 0, 0);
        if (tpcGeo.DetectDriftDirection() == -1) elecvec.SetXYZ(-1, 0, 0);
        if (tpcGeo.DetectDriftDirection() == 2) elecvec.SetXYZ(0, 1, 0);
        if (tpcGeo.DetectDriftDirection() == -2) elecvec.SetXYZ(0, -1, 0);
        if (tpcGeo.DetectDriftDirection() == 3) elecvec.SetXYZ(0, 0, 1);
        if (tpcGeo.DetectDriftDirection() == -3) elecvec.SetXYZ(0, 0, -1);

        elecvec *= efield;

        // electric field inside active volume
        if (fSCE->EnableSimEfieldSCE()) {
        auto const eFieldOffsets = fSCE->GetEfieldOffsets(edep.MidPoint());

        TVector3 scevec;

        if (tpcGeo.DetectDriftDirection() == 1)
            scevec.SetXYZ(
            efield * eFieldOffsets.X(), efield * eFieldOffsets.Y(), efield * eFieldOffsets.Z());
        if (tpcGeo.DetectDriftDirection() == -1)
            scevec.SetXYZ(
            -1 * efield * eFieldOffsets.X(), efield * eFieldOffsets.Y(), efield * eFieldOffsets.Z());
        if (tpcGeo.DetectDriftDirection() == 2)
            scevec.SetXYZ(
            efield * eFieldOffsets.X(), efield * eFieldOffsets.Y(), efield * eFieldOffsets.Z());
        if (tpcGeo.DetectDriftDirection() == -2)
            scevec.SetXYZ(
            efield * eFieldOffsets.X(), -1 * efield * eFieldOffsets.Y(), efield * eFieldOffsets.Z());
        if (tpcGeo.DetectDriftDirection() == 3)
            scevec.SetXYZ(
            efield * eFieldOffsets.X(), efield * eFieldOffsets.Y(), efield * eFieldOffsets.Z());
        if (tpcGeo.DetectDriftDirection() == -3)
            scevec.SetXYZ(
            efield * eFieldOffsets.X(), efield * eFieldOffsets.Y(), -1 * efield * eFieldOffsets.Z());

        elecvec += scevec;
        }

        double angle = std::acos(stepvec.Dot(elecvec) / (stepvec.Mag() * elecvec.Mag()));

        if (angle > TMath::PiOver2()) { angle = abs(TMath::Pi() - angle); }

        return angle;
    }

} // namespace larg4