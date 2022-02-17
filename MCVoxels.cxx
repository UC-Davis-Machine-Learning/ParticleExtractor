/**
 * @file MCVoxels.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-02-17
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "MCVoxels.h"

namespace extractor
{
    MCVoxels::MCVoxels()
    {
        fMCVoxelsTree = fTFileService->make<TTree>("mc_voxels", "mc_voxels");
        fMCVoxelsTree->Branch("voxels", &fVoxels.voxels);
        fMCVoxelsTree->Branch("labels", &fVoxels.labels);
        fMCVoxelsTree->Branch("energy", &fVoxels.energy);
        fMCVoxelsTree->Branch("edep_idxs", &fVoxels.edep_idxs);
    }

    MCVoxels::~MCVoxels()
    {}

    void MCVoxels::setBoundingBox(std::string boundingBox)
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

    void MCVoxels::setVoxelLabeling(std::string voxelLabeling)
    {
        if (voxelLabeling == "largest"){
            fVoxelLabeling = 0;
        }
        else {
            fVoxelLabeling = 1;
            fMixedLabel = *std::max_element(fVoxelLabels.begin(), fVoxelLabels.end()) + 1;
        }
    }

    void MCVoxels::processEvent(const MCEnergyDeposits& energyDeposits)
    {
        Voxels voxels;
        std::vector<std::vector<Int_t>> temp_voxels;
        std::vector<std::vector<Int_t>> temp_labels;
        std::vector<std::vector<Double_t>> temp_energy;
        std::vector<std::vector<Int_t>> temp_edep_idxs;

        std::vector<Int_t> xyz(3);
        Int_t label;
        Double_t energy;
        Int_t edep_idx;

        MCEdep mcEdep = energyDeposits.getMCEdep();
        /**
         * Now we loop through all of the energy depositions in energyDeposits
         * and construct unique voxels.  For now, we store the labels, energy
         * and edep_idxs for each edep that lands in each voxel, and later decide
         * how to deal with them.
         */
        for (size_t i = 0; i < mcEdep.pdg.size(); i++)
        {
            // gather variables
            xyz[0] = int((mcEdep.edep_x[i] - fBoundingBox.x_min) / fVoxelSize);
            xyz[1] = int((mcEdep.edep_y[i] - fBoundingBox.y_min) / fVoxelSize);
            xyz[2] = int((mcEdep.edep_z[i] - fBoundingBox.z_min) / fVoxelSize);
            label = mcEdep.pdg[i];
            energy = mcEdep.energy[i];
            edep_idx = static_cast<int>(i);

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
            }
            else
            {
                Int_t voxel_index = std::distance(temp_voxels.begin(), voxel_exists);
                temp_labels[voxel_index].emplace_back(label);
                temp_energy[voxel_index].emplace_back(energy);
                temp_edep_idxs[voxel_index].emplace_back(edep_idx);
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
            std::vector<Int_t> unique_energy;
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
                }
            }
            else {
                voxels.labels.emplace_back(unique_labels[0]);
            }
        }
        fVoxels = voxels;
        fMCVoxelsTree->Fill();
    }
}