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

#include <vector>
#include <algorithm>

#include "DetectorGeometry.h"

namespace extractor
{
    /**
     * @brief 
     * 
     */
    struct MCEdep
    {
        std::vector<Int_t> pdg;
        std::vector<Int_t> track_id;
        std::vector<Int_t> ancestor_id;
        std::vector<Int_t> level;
        std::vector<Double_t> edep_x;
        std::vector<Double_t> edep_y;
        std::vector<Double_t> edep_z;
        std::vector<Double_t> energy;
        std::vector<Int_t> num_electrons;
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

        void setBoundingBoxType(VolumeType volumeType) { fBoundingBoxType = volumeType; }
        void setBoundingBoxType(std::string volumeType);
        void setPDGCodes(std::vector<Int_t> PDGCodes) { fPDGCodes = PDGCodes; }
        void setPDGLevels(std::vector<Int_t> PDGLevels) { fPDGLevels = PDGLevels; }
        void setPDGLevels(std::vector<std::string> PDGLevels);
        void setEnergyCutoff(Double_t energyCutoff) { fEnergyCutoff = energyCutoff; }

        void processEvent(
            const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
            const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>& mcEnergyDeposits
        );

        MCEdep getMCEdep() const { return fMCEdep; }

    private:
        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.
         */ 
        art::ServiceHandle<art::TFileService> fTFileService;
        TTree *fMCEnergyDepositsTree;
        // geometry information
        DetectorGeometry* fGeometry = DetectorGeometry::getInstance("MCEnergyDeposits");

        // pdg codes to construct
        VolumeType fBoundingBoxType;
        std::vector<Int_t> fPDGCodes;
        std::vector<Int_t> fPDGLevels;
        Double_t fEnergyCutoff;

        // struct for holding event information
        MCEdep fMCEdep;
    };
}