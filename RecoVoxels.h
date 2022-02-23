/**
 * @file RecoVoxels.h
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
#include <numeric>
#include <map>

#include "DetectorGeometry.h"
#include "RecoEnergyDeposits.h"
#include "MCVoxels.h"

namespace extractor
{
    class RecoVoxels
    {
    public:
        RecoVoxels();
        ~RecoVoxels();

        void setPDGCodes(std::vector<Int_t> pdgCodes) { fPDGCodes = pdgCodes; }
        void setVoxelLabels(std::vector<Int_t> voxelLabels);
        void setVoxelSize(Double_t voxelSize) { fVoxelSize = voxelSize; }
        void setBoundingBox(BoundingBox boundingBox) { fBoundingBox = boundingBox; }
        void setBoundingBox(std::string boundingBox);
        void setVoxelLabeling(Int_t voxelLabeling) { fVoxelLabeling = voxelLabeling; }
        void setVoxelLabeling(std::string voxelLabeling);

        void processEvent(const RecoEnergyDeposits& energyDeposits);
    private:
        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.
         */ 
        art::ServiceHandle<art::TFileService> fTFileService;
        TTree *fRecoVoxelsTree;
        // geometry information
        DetectorGeometry* fGeometry = DetectorGeometry::getInstance("MCVoxels");

        // parameters
        BoundingBox fBoundingBox;
        std::vector<Int_t> fPDGCodes;
        std::vector<Int_t> fVoxelLabels;
        Double_t fVoxelSize;
        Int_t fVoxelLabeling;
        Int_t fMixedLabel;
        std::map<Int_t, Int_t> fPDGLabelMap;

        // temp struct
        Voxels fVoxels;
    };
}