/**
 * @file    Configuration.h
 * @author  Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief   A struct for holding LArSoft configuration parameters
 *          for the IonAndScint module.
 * @version 0.1
 * @date 2022-02-15
 */
#pragma once

// art framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace larg4
{
    struct IonAndScintConfiguration
    {
        fhicl::Atom<std::string> NESTMode
        {
            fhicl::Name("NESTMode"),
            fhicl::Comment("Mode for the NEST algorithm, can be Birks/Box/Legacy/NEST/Combined.")
        };
        // fhicl::Atom<Double_t> ShowerEnergyThreshold
        // {
        //     fhicl::Name("ShowerEnergyThreshold"),
        //     fhicl::Comment("The amount of energy required for an ionized electron to be considered a shower, rather than a blip.")
        // };
        // fhicl::Atom<bool> VoxelizeEnergyDeposits
        // {
        //     fhicl::Name("VoxelizeEnergyDeposits"),
        //     fhicl::Comment("Whether to save edep point clouds with true (x,y,z) or to voxelize them.")
        // };
    };

    using Parameters = art::EDAnalyzer::Table<IonAndScintConfiguration>;
}