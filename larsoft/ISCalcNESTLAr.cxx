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
        fNESTMode = pset.get<art::InputTag>("nest_mode");
        larnest::LArNRYieldsParameters NRYieldsParameters;
        larnest::LArERElectronYieldsAlphaParameters ERElectronYieldsAlphaParameters;
        larnest::LArERElectronYieldsBetaParameters ERElectronYieldsBetaParameters;
        larnest::LArERElectronYieldsGammaParameters ERElectronYieldsGammaParameters;
        larnest::LArERElectronYieldsDokeBirksParameters ERElectronYieldsDokeBirksParameters;
        larnest::LArERYieldsParameters ERYieldsParameters;
        larnest::LArAlphaElectronYieldsParameters AlphaElectronYieldsParameters;
        larnest::LArAlphaPhotonYieldsParameters AlphaPhotonYieldsParameters;
        larnest::LArAlphaYieldsParameters AlphaYieldsParameters;
        larnest::LArdEdxParameters dEdxParameters;
        larnest::BOXParameters BOXParameters;
        larnest::BIRKSParameters BIRKSParameters;

        NRYieldsParameters.alpha = pset.get<Double_t>("nest_nr_yields_parameters_alpha");
        NRYieldsParameters.beta = pset.get<Double_t>("nest_nr_yields_parameters_beta");
        NRYieldsParameters.gamma = pset.get<Double_t>("nest_nr_yields_parameters_gamma");
        NRYieldsParameters.delta = pset.get<Double_t>("nest_nr_yields_parameters_delta");
        NRYieldsParameters.epsilon = pset.get<Double_t>("nest_nr_yields_parameters_epsilon");
        NRYieldsParameters.zeta = pset.get<Double_t>("nest_nr_yields_parameters_zeta");
        NRYieldsParameters.eta = pset.get<Double_t>("nest_nr_yields_parameters_eta");

        ERElectronYieldsAlphaParameters.A = pset.get<Double_t>("nest_er_electron_yields_alpha_parameters_A");
        ERElectronYieldsAlphaParameters.B = pset.get<Double_t>("nest_er_electron_yields_alpha_parameters_B");
        ERElectronYieldsAlphaParameters.C = pset.get<Double_t>("nest_er_electron_yields_alpha_parameters_C");
        ERElectronYieldsAlphaParameters.D = pset.get<Double_t>("nest_er_electron_yields_alpha_parameters_D");
        ERElectronYieldsAlphaParameters.E = pset.get<Double_t>("nest_er_electron_yields_alpha_parameters_E");
        ERElectronYieldsAlphaParameters.F = pset.get<Double_t>("nest_er_electron_yields_alpha_parameters_F");
        ERElectronYieldsAlphaParameters.G = pset.get<Double_t>("nest_er_electron_yields_alpha_parameters_G");

        ERElectronYieldsBetaParameters.A = pset.get<Double_t>("nest_er_electron_yields_beta_parameters_A");
        ERElectronYieldsBetaParameters.B = pset.get<Double_t>("nest_er_electron_yields_beta_parameters_B");
        ERElectronYieldsBetaParameters.C = pset.get<Double_t>("nest_er_electron_yields_beta_parameters_C");
        ERElectronYieldsBetaParameters.D = pset.get<Double_t>("nest_er_electron_yields_beta_parameters_D");
        ERElectronYieldsBetaParameters.E = pset.get<Double_t>("nest_er_electron_yields_beta_parameters_E");
        ERElectronYieldsBetaParameters.F = pset.get<Double_t>("nest_er_electron_yields_beta_parameters_F");

        ERElectronYieldsGammaParameters.A = pset.get<Double_t>("nest_er_electron_yields_gamma_parameters_A");
        ERElectronYieldsGammaParameters.B = pset.get<Double_t>("nest_er_electron_yields_gamma_parameters_B");
        ERElectronYieldsGammaParameters.C = pset.get<Double_t>("nest_er_electron_yields_gamma_parameters_C");
        ERElectronYieldsGammaParameters.D = pset.get<Double_t>("nest_er_electron_yields_gamma_parameters_D");
        ERElectronYieldsGammaParameters.E = pset.get<Double_t>("nest_er_electron_yields_gamma_parameters_E");
        ERElectronYieldsGammaParameters.F = pset.get<Double_t>("nest_er_electron_yields_gamma_parameters_F");
        ERElectronYieldsGammaParameters.G = pset.get<Double_t>("nest_er_electron_yields_gamma_parameters_G");

        ERElectronYieldsDokeBirksParameters.A = pset.get<Double_t>("nest_er_electron_yields_doke_birks_parameters_A");
        ERElectronYieldsDokeBirksParameters.B = pset.get<Double_t>("nest_er_electron_yields_doke_birks_parameters_B");
        ERElectronYieldsDokeBirksParameters.C = pset.get<Double_t>("nest_er_electron_yields_doke_birks_parameters_C");
        ERElectronYieldsDokeBirksParameters.D = pset.get<Double_t>("nest_er_electron_yields_doke_birks_parameters_D");
        ERElectronYieldsDokeBirksParameters.E = pset.get<Double_t>("nest_er_electron_yields_doke_birks_parameters_E");

        ERYieldsParameters.alpha = ERElectronYieldsAlphaParameters;
        ERYieldsParameters.beta = ERElectronYieldsBetaParameters;
        ERYieldsParameters.gamma = ERElectronYieldsGammaParameters;
        ERYieldsParameters.doke_birks = ERElectronYieldsDokeBirksParameters;
        ERYieldsParameters.p1 = pset.get<Double_t>("nest_er_yields_parameters_p1");
        ERYieldsParameters.p2 = pset.get<Double_t>("nest_er_yields_parameters_p2");
        ERYieldsParameters.p3 = pset.get<Double_t>("nest_er_yields_parameters_p3");
        ERYieldsParameters.p4 = pset.get<Double_t>("nest_er_yields_parameters_p4");
        ERYieldsParameters.p5 = pset.get<Double_t>("nest_er_yields_parameters_p5");
        ERYieldsParameters.delta = pset.get<Double_t>("nest_er_yields_parameters_delta");
        ERYieldsParameters.let = pset.get<Double_t>("nest_er_yields_parameters_let");

        AlphaElectronYieldsParameters.A = pset.get<Double_t>("nest_alpha_electron_yields_parameters_A");
        AlphaElectronYieldsParameters.B = pset.get<Double_t>("nest_alpha_electron_yields_parameters_B");
        AlphaElectronYieldsParameters.C = pset.get<Double_t>("nest_alpha_electron_yields_parameters_C");
        AlphaElectronYieldsParameters.D = pset.get<Double_t>("nest_alpha_electron_yields_parameters_D");
        AlphaElectronYieldsParameters.E = pset.get<Double_t>("nest_alpha_electron_yields_parameters_E");
        AlphaElectronYieldsParameters.F = pset.get<Double_t>("nest_alpha_electron_yields_parameters_F");
        AlphaElectronYieldsParameters.G = pset.get<Double_t>("nest_alpha_electron_yields_parameters_G");
        AlphaElectronYieldsParameters.H = pset.get<Double_t>("nest_alpha_electron_yields_parameters_H");
        AlphaElectronYieldsParameters.I = pset.get<Double_t>("nest_alpha_electron_yields_parameters_I");
        AlphaElectronYieldsParameters.J = pset.get<Double_t>("nest_alpha_electron_yields_parameters_J");

        AlphaPhotonYieldsParameters.A = pset.get<Double_t>("nest_alpha_photon_yields_parameters_A");
        AlphaPhotonYieldsParameters.B = pset.get<Double_t>("nest_alpha_photon_yields_parameters_B");
        AlphaPhotonYieldsParameters.C = pset.get<Double_t>("nest_alpha_photon_yields_parameters_C");
        AlphaPhotonYieldsParameters.D = pset.get<Double_t>("nest_alpha_photon_yields_parameters_D");
        AlphaPhotonYieldsParameters.E = pset.get<Double_t>("nest_alpha_photon_yields_parameters_E");
        AlphaPhotonYieldsParameters.F = pset.get<Double_t>("nest_alpha_photon_yields_parameters_F");
        AlphaPhotonYieldsParameters.G = pset.get<Double_t>("nest_alpha_photon_yields_parameters_G");
        AlphaPhotonYieldsParameters.H = pset.get<Double_t>("nest_alpha_photon_yields_parameters_H");
        AlphaPhotonYieldsParameters.I = pset.get<Double_t>("nest_alpha_photon_yields_parameters_I");
        AlphaPhotonYieldsParameters.J = pset.get<Double_t>("nest_alpha_photon_yields_parameters_J");
        AlphaPhotonYieldsParameters.K = pset.get<Double_t>("nest_alpha_photon_yields_parameters_K");
        AlphaPhotonYieldsParameters.L = pset.get<Double_t>("nest_alpha_photon_yields_parameters_L");
        AlphaPhotonYieldsParameters.M = pset.get<Double_t>("nest_alpha_photon_yields_parameters_M");
        AlphaYieldsParameters.Ye = AlphaElectronYieldsParameters;
        AlphaYieldsParameters.Yph = AlphaPhotonYieldsParameters;

        BOXParameters.alpha = pset.get<Double_t>("nest_box_parameters_alpha");
        BOXParameters.beta = pset.get<Double_t>("nest_box_parameters_beta");
        BIRKSParameters.Ab = pset.get<Double_t>("nest_birks_parameters_Ab");
        BIRKSParameters.kb = pset.get<Double_t>("nest_birks_parameters_kb");

        mLArNEST.SetNRYieldsParameters(NRYieldsParameters);
        mLArNEST.SetERYieldsParameters(ERYieldsParameters);
        mLArNEST.SetAlphaYieldsParameters(AlphaYieldsParameters);
        mLArNEST.SetBOXParameters(BOXParameters);
        mLArNEST.SetBIRKSParameters(BIRKSParameters);
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
        larnest::LArInteraction species;
        if (fNESTMode == "Box") {
            species = larnest::LArInteraction::BOX;
        }
        else {
            species = larnest::LArInteraction::BIRKS;
        }

        larnest::LArYieldResult yields = mLArNEST.GetYields(
            species, 
            energy_deposit * 1000,  // NEST expects keV 
            dx, 
            EFieldStep, 
            1.39
        );
        double NumElectrons = yields.Ne;
        double NumPhotons = yields.Nph;

        std::cout << "HERE IN NEST! num_e: " << NumElectrons << ", num_ph: " << NumPhotons << std::endl;

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