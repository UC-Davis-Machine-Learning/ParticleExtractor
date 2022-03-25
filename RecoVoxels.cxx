/**
 * @file RecoVoxels.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-02-17
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "RecoVoxels.h"

namespace extractor
{
    RecoVoxels::RecoVoxels()
    {
        fRecoVoxelsTree = fTFileService->make<TTree>("reco_voxels", "reco_voxels");
        fRecoVoxelsTree->Branch("voxels", &fVoxels.voxels);
        fRecoVoxelsTree->Branch("labels", &fVoxels.labels);
        fRecoVoxelsTree->Branch("energy", &fVoxels.energy);
        fRecoVoxelsTree->Branch("edep_idxs", &fVoxels.edep_idxs);
        fRecoVoxelsTree->Branch("levels", &fVoxels.levels);
    }

    RecoVoxels::~RecoVoxels()
    {}

    void RecoVoxels::setBoundingBox(std::string boundingBox)
    {
        if (boundingBox == "TPC" or boundingBox == "tpc") { 
            fBoundingBox = fGeometry->GetTotalActiveTPCBox();
        }
        else if (boundingBox == "Cryo" or boundingBox == "cryo") {
            fBoundingBox = fGeometry->GetCryostatBox();
        }
        else {
            fBoundingBox = fGeometry->GetWorldBox();
        }
    }

    void RecoVoxels::setVoxelLabels(std::vector<Int_t> voxelLabels)
    {
        fVoxelLabels = voxelLabels;
        for (size_t i = 0; i < fVoxelLabels.size(); i++)
        {
            fPDGLabelMap[fPDGCodes[i]] = fVoxelLabels[i];
        }
    }

    void RecoVoxels::setVoxelLabeling(std::string voxelLabeling)
    {
        if (voxelLabeling == "largest"){
            fVoxelLabeling = 0;
        }
        else {
            fVoxelLabeling = 1;
            fMixedLabel = *std::max_element(fVoxelLabels.begin(), fVoxelLabels.end()) + 1;
        }
    }

    void RecoVoxels::processEvent(const RecoEnergyDeposits& energyDeposits)
    {
        Voxels voxels;
        std::vector<std::vector<Int_t>> temp_voxels;
        std::vector<std::vector<Int_t>> temp_labels;
        std::vector<std::vector<Double_t>> temp_energy;
        std::vector<std::vector<Int_t>> temp_edep_idxs;
        std::vector<std::vector<Int_t>> temp_levels;

        std::vector<Int_t> xyz(3);
        Int_t label;
        Double_t energy;
        Int_t edep_idx;
        Int_t level;

        RecoEdep recoEdep = energyDeposits.getRecoEdep();
        /**
         * Now we loop through all of the energy depositions in energyDeposits
         * and construct unique voxels.  For now, we store the labels, energy
         * and edep_idxs for each edep that lands in each voxel, and later decide
         * how to deal with them.
         */
        for (size_t i = 0; i < recoEdep.pdg.size(); i++)
        {
            for (size_t j = 0; j < recoEdep.pdg[i].size(); j++)
            {
                // gather variables
                xyz[0] = int((recoEdep.sp_x[i] - fBoundingBox.x_min) / fVoxelSize);
                xyz[1] = int((recoEdep.sp_y[i] - fBoundingBox.y_min) / fVoxelSize);
                xyz[2] = int((recoEdep.sp_z[i] - fBoundingBox.z_min) / fVoxelSize);
                // find if pdg is in the map, otherwise use -1
                auto it = fPDGLabelMap.find(recoEdep.pdg[i][j]);
                if (it != fPDGLabelMap.end()) {
                    label = fPDGLabelMap[recoEdep.pdg[i][j]];
                }
                else {
                    label = -1;
                }
                energy = recoEdep.summed_adc[i][j];
                edep_idx = static_cast<int>(i);
                level = recoEdep.level[i][j];

                // see if xyz is in temp_voxels
                auto voxel_exists = std::find(
                    temp_voxels.begin(), 
                    temp_voxels.end(), 
                    xyz
                );
                if (voxel_exists == temp_voxels.end())
                {
                    temp_voxels.emplace_back(xyz);
                    temp_labels.emplace_back(std::vector<Int_t>({label}));
                    temp_energy.emplace_back(std::vector<Double_t>({energy}));
                    temp_edep_idxs.emplace_back(std::vector<Int_t>({edep_idx}));
                    temp_levels.emplace_back(std::vector<Int_t>({level}));
                }
                else
                {
                    Int_t voxel_index = std::distance(temp_voxels.begin(), voxel_exists);
                    temp_labels[voxel_index].emplace_back(label);
                    temp_energy[voxel_index].emplace_back(energy);
                    temp_edep_idxs[voxel_index].emplace_back(edep_idx);
                    temp_levels[voxel_index].emplace_back(level);
                }
            }
        }
        /**
         * Now we itrerate through the unique voxels once more to decide
         * how to label them.  If the labeling is "largest", then we
         * simply assign the label corresponding to the largest amount
         * of energy deposited in the voxel, otherwise we assign a "mixed"
         * labeling.
         */
        for (size_t i = 0; i < temp_voxels.size(); i++)
        {
            std::vector<Int_t> unique_labels;
            std::vector<Double_t> unique_energy;
            std::vector<Int_t> unique_levels;
            // collect unique labels and the total energy
            for (size_t j = 0; j < temp_labels[i].size(); j++)
            {
                auto label_exists = std::find(
                    unique_labels.begin(), 
                    unique_labels.end(), 
                    temp_labels[i][j]
                );
                if (label_exists == unique_labels.end())
                {
                    unique_labels.emplace_back(temp_labels[i][j]);
                    unique_energy.emplace_back(temp_energy[i][j]);
                    unique_levels.emplace_back(temp_levels[i][j]);
                }
                else
                {
                    Int_t label_index = std::distance(unique_labels.begin(), label_exists);
                    unique_energy[label_index] += temp_energy[i][j];      
                }
            }
            voxels.voxels.emplace_back(temp_voxels[i]);
            voxels.energy.emplace_back(
                std::accumulate(unique_energy.begin(), unique_energy.end(), 0.)
            );
            voxels.edep_idxs.emplace_back(temp_edep_idxs[i]);

            if (unique_labels.size() > 1) 
            {   
                // if using a mixed labeling
                if (fVoxelLabeling == 1) {
                    voxels.labels.emplace_back(fMixedLabel);  
                }
                else
                {
                    Int_t arg_max_energy = std::distance(
                        unique_energy.begin(), 
                        std::max_element(unique_energy.begin(), unique_energy.end())
                    );
                    voxels.labels.emplace_back(unique_labels[arg_max_energy]);
                    voxels.levels.emplace_back(unique_levels[arg_max_energy]);
                }
            }
            else {
                voxels.labels.emplace_back(unique_labels[0]);
            }
        }
        fVoxels = voxels;
        fRecoVoxelsTree->Fill();
    }
}