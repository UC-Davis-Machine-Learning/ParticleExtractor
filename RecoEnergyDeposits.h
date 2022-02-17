/**
 * @file RecoEnergyDeposits.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-02-17
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
    struct RecoEdep
    {
        std::vector<Int_t> pdg;
        std::vector<Int_t> track_id;
        std::vector<Int_t> channel_id;
        std::vector<Double_t> sp_x;
        std::vector<Double_t> sp_y;
        std::vector<Double_t> sp_z;
        std::vector<Double_t> energy;
        std::vector<Int_t> num_electrons;
    };

    /**
     * @brief 
     * 
     */
    class RecoEnergyDeposits
    {
    public:
        RecoEnergyDeposits();
        ~RecoEnergyDeposits();

        void setPDGCodes(std::vector<Int_t> PDGCodes) { fPDGCodes = PDGCodes; }
        void setPDGLevels(std::vector<Int_t> PDGLevels) { fPDGLevels = PDGLevels; }
        void setPDGLevels(std::vector<std::string> PDGLevels);

        void processEvent(
            const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
            const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
            const art::ValidHandle<std::vector<recob::Hit>>& recoHits,
            const art::ValidHandle<std::vector<recob::SpacePoint>>& recoSpacePoints,
            const art::FindManyP<recob::Hit>& hitSpacePointAssn
        );

        RecoEdep getRecoEdep() const { return fRecoEdep; }

    private:
        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.
         */ 
        art::ServiceHandle<art::TFileService> fTFileService;
        TTree *fRecoEnergyDepositsTree;
        // geometry information
        DetectorGeometry* fGeometry = DetectorGeometry::getInstance("RecoEnergyDeposits");

        // pdg codes to construct
        std::vector<Int_t> fPDGCodes;
        std::vector<Int_t> fPDGLevels;

        // struct for holding event information
        RecoEdep fRecoEdep;
    };
}