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
         * LAr Geant4 and reco configuration parameters.
         * The user must specify the name of the larg4 module
         * that was used to run the geant simulations.
         * This is used to get simb::MCParticle products.
         */
        fhicl::Atom<art::InputTag> LArGeantProducerLabel
        {
            fhicl::Name("LArGeantProducerLabel"),
            fhicl::Comment("Tag of the input data product for the largeant side of the simulation.")
        };
        fhicl::Atom<art::InputTag> IonAndScintProducerLabel
        {
            fhicl::Name("IonAndScintProducerLabel"),
            fhicl::Comment("Tag of the input data product for the IonAndScint side of the simulation.")
        };
        fhicl::Atom<art::InputTag> SimChannelProducerLabel
        {
            fhicl::Name("SimChannelProducerLabel"),
            fhicl::Comment("Tag of the input data product for the SimChannel side of the simulation.")
        };
        fhicl::Atom<art::InputTag> SimChannelInstanceProducerLabel
        {
            fhicl::Name("SimChannelInstanceProducerLabel"),
            fhicl::Comment("Tag of the input data product for the SimChannel instance side of the simulation.")
        };
        fhicl::Atom<art::InputTag> HitProducerLabel
        {
            fhicl::Name("HitProducerLabel"),
            fhicl::Comment("Tag of the input data product for the Hit side of the simulation.")
        };
        fhicl::Atom<art::InputTag> SpacePointProducerLabel
        {
            fhicl::Name("SpacePointProducerLabel"),
            fhicl::Comment("Tag of the input data product for the SpacePoint side of the simulation.")
        };

        fhicl::Atom<art::InputTag> PandoraLabel
        {
            fhicl::Name("PandoraLabel"),
            fhicl::Comment("Tag of the input data product for the SpacePoints from Pandora.")
        };

        fhicl::Atom<art::InputTag> PandoraTrackLabel
        {
            fhicl::Name("PandoraTrackLabel"),
            fhicl::Comment("Tag of the input data product for the track hits from Pandora.")
        };

        /**
         * These options generate TTrees called "mc_neutron_captures",
         * "mc_energy_deposits", "reco_energy_deposits",
         * which stores various neutron capture statistics,
         * and energy deposition information from each event.
         */
        fhicl::Atom<bool> FillMCNeutronCaptures
        {
            fhicl::Name("FillMCNeutronCaptures"),
            fhicl::Comment("Whether to save neutron capture locations.")
        };
        fhicl::Atom<bool> FillMCEnergyDeposits
        {
            fhicl::Name("FillMCEnergyDeposits"),
            fhicl::Comment("Whether to save MC edep information.")
        };
        fhicl::Atom<bool> FillRecoEnergyDeposits
        {
            fhicl::Name("FillRecoEnergyDeposits"),
            fhicl::Comment("Whether to save Reco edep information.")
        };
        fhicl::Atom<std::string> RecoEdepBoundingBox
        {
            fhicl::Name("RecoEdepBoundingBox"),
            fhicl::Comment("Which bounding box to use for the Reco edeps.")
        };

        /**
         * Below are a set of parameters for extracting energy deposition
         * information from MC truth.
         * Set of pdg codes to extract from each event
         */
        fhicl::Atom<std::string> MCEdepBoundingBox
        {
            fhicl::Name("MCEdepBoundingBox"),
            fhicl::Comment("Which bounding box to use for the MC edeps.")
        };
        fhicl::Sequence<int> MCEdepPDGCodes
        {
            fhicl::Name("MCEdepPDGCodes"),
            fhicl::Comment("PDG IDs of the particles to extracted.")
        };
        fhicl::Sequence<std::string> MCEdepPDGLevels
        {
            fhicl::Name("MCEdepPDGLevels"),
            fhicl::Comment("Heirarchy of particles to keep.")
        };
        fhicl::Atom<double> MCEdepEnergyCutoff
        {
            fhicl::Name("MCEdepEnergyCutoff"),
            fhicl::Comment("Cutoff for storing Edeps.")
        };
        /**
         * Below are a set of voxelization parameters for MC energy depositions.
         * One must specify a voxel size, as well as a bounding box to use
         * for each event (typically just the active tpc volume).
         * 
         */
        fhicl::Atom<bool> FillMCVoxels
        {
            fhicl::Name("FillMCVoxels"),
            fhicl::Comment("Whether to save MC voxel information.")
        };
        fhicl::Sequence<int> MCEdepPDGLabels
        {
            fhicl::Name("MCEdepPDGLabels"),
            fhicl::Comment("Labeling scheme for the PDGCodes.")
        };
        fhicl::Atom<double> MCVoxelSize
        {
            fhicl::Name("MCVoxelSize"),
            fhicl::Comment("Size of the voxels in mm.")
        };
        fhicl::Atom<std::string> MCVoxelBoundingBox
        {
            fhicl::Name("MCVoxelBoundingBox"),
            fhicl::Comment("Which bounding box to use for the voxelization.")
        };
        fhicl::Atom<std::string> MCVoxelLabeling
        {
            fhicl::Name("MCVoxelLabeling"),
            fhicl::Comment("Labeling scheme for the voxels.")
        };

        /**
         * 
         */
        fhicl::Atom<bool> FillRawDecoder
        {
            fhicl::Name("FillRawDecoder"),
            fhicl::Comment("Whether to save MC voxel information.")
        };
        fhicl::Sequence<int> RawDecoderPDGCodes
        {
            fhicl::Name("RawDecoderPDGCodes"),
            fhicl::Comment("PDG IDs of the particles to extracted.")
        };
        fhicl::Sequence<std::string> RawDecoderPDGLevels
        {
            fhicl::Name("RawDecoderPDGLevels"),
            fhicl::Comment("Heirarchy of particles to keep.")
        };
        fhicl::Atom<double> RawDecoderEnergyCutoff
        {
            fhicl::Name("RawDecoderEnergyCutoff"),
            fhicl::Comment("Cutoff for storing Edeps.")
        };

        /**
         * Below are a set of voxelization parameters for Reco energy depositions.
         */
        fhicl::Atom<bool> FillRecoVoxels
        {
            fhicl::Name("FillRecoVoxels"),
            fhicl::Comment("Whether to save Reco voxel information.")
        };
        fhicl::Sequence<int> RecoEdepPDGLabels
        {
            fhicl::Name("RecoEdepPDGLabels"),
            fhicl::Comment("Labeling scheme for the PDGCodes.")
        };
        fhicl::Atom<double> RecoVoxelSize
        {
            fhicl::Name("RecoVoxelSize"),
            fhicl::Comment("Size of the voxels in mm.")
        };
        fhicl::Atom<std::string> RecoVoxelBoundingBox
        {
            fhicl::Name("RecoVoxelBoundingBox"),
            fhicl::Comment("Which bounding box to use for the voxelization.")
        };
        fhicl::Atom<std::string> RecoVoxelLabeling
        {
            fhicl::Name("RecoVoxelLabeling"),
            fhicl::Comment("Labeling scheme for the voxels.")
        };

        // RecoTracks

        fhicl::Atom<bool> FillRecoTracks
        {
            fhicl::Name("FillRecoTracks"),
            fhicl::Comment("Whether to save Reco tracks information.")
        };
    };

    using Parameters = art::EDAnalyzer::Table<Configuration>;

    // allowed values of MCEdepPDGType
    std::vector<std::string> allowed_mc_edep_levels = 
    {
        "parent",
        "daughters",
        "electrons",
        "parent_electrons",
        "all"
    };

    std::vector<std::string> allowed_mc_voxel_labeling = 
    {
        "largest",
        "mixed"  
    };

    std::vector<std::string> allowed_bounding_boxes = 
    {
        "TPC",
        "tpc",
        "cryo",
        "Cryo",
        "World",
        "world"
    };
}