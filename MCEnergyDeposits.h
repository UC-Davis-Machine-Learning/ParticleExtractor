/**
 * @file MCEnergyDeposits.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-02-16
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "DetectorGeometry.h"

namespace extractor
{
    /**
     * @brief 
     * 
     */
    struct MCEdep
    {

    };

    /**
     * @brief 
     * 
     */
    class MCEnergyDeposits
    {
    public:
        MCEnergyDeposits();
        ~MCEnergyDeposits();

        void processEvent(
            art::ValidHandle<std::vector<simb::MCParticle>> mcParticles,
            art::ValidHandle<std::vector<sim::SimEnergyDeposit>> mcEnergyDeposits
        );

    private:
        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.
         */ 
        art::ServiceHandle<art::TFileService> fTFileService;
        TTree *fMCNeutronCapturesTree;
        // geometry information
        DetectorGeometry* fGeometry = DetectorGeometry::getInstance("MCEnergyDeposits");
        
        // struct for holding event information
        MCEdep fMCEdep;
    };
}