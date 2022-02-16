/**
 * @file    Configuration.h
 * @author  Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief   A struct for holding LArSoft configuration parameters
 *          for the ParticleExtractor module.
 * @version 0.1
 * @date 2022-02-15
 * 
 * @copyright Copyright (c) 2022
 * 
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

namespace extractor
{
    /**
     * @brief A collection of fhicl parameters for the ParticleExtractor
     * module.  Each of these must be specified, which have default values
     * in the accompanying fhicl file.
     */
    struct Configuration
    {
        /**
         * @brief This option generates a TTree called "mc_neutron_captures",
         * which stores various neutron capture statistics from each event.
         * 
         */
        fhicl::Atom<bool> FillMCNeutronCaptures
        {
            fhicl::Name("FillMCNeutronCaptures"),
            fhicl::Comment("Whether to save neutron capture locations.")
        };

        /**
         * @brief LAr Geant4 configuration parameters.
         * The user must specify the name of the larg4 module
         * that was used to run the geant simulations.
         * This is used to get simb::MCParticle products.
         */
        fhicl::Atom<art::InputTag> LArGeantProducerLabel
        {
            fhicl::Name("LArGeantProducerLabel"),
            fhicl::Comment("Tag of the input data product for the largeant side of the simulation.")
        };

        /**
         * @brief Set of pdg codes to extract from each event
         * 
         */
        fhicl::Sequence<int> PDGCodes
        {
            fhicl::Name("PDGCodes"),
            fhicl::Comment("PDG IDs of the particles to extracted.")
        };

    };

    using Parameters = art::EDAnalyzer::Table<Configuration>;
}