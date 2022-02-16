/**
 * @file MCEnergyDeposits.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-02-16
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "MCEnergyDeposits.h"

namespace extractor
{
    MCEnergyDeposits::MCEnergyDeposits()
    {
        fMCEnergyDepositsTree->Branch("pdg", &fMCEdep.pdg);
        fMCEnergyDepositsTree->Branch("track_id", &fMCEdep.track_id);
        fMCEnergyDepositsTree->Branch("ancestor_id", &fMCEdep.ancestor_id);
        fMCEnergyDepositsTree->Branch("level", &fMCEdep.level);
        fMCEnergyDepositsTree->Branch("edep_x", &fMCEdep.edep_x);
        fMCEnergyDepositsTree->Branch("edep_y", &fMCEdep.edep_y);
        fMCEnergyDepositsTree->Branch("edep_z", &fMCEdep.edep_z);
        fMCEnergyDepositsTree->Branch("energy", &fMCEdep.energy);
        fMCEnergyDepositsTree->Branch("num_electrons", &fMCEdep.num_electrons);
    }

    MCEnergyDeposits::~MCEnergyDeposits()
    {}

    void MCEnergyDeposits::setPDGLevels(std::vector<std::string> PDGLevels)
    {
        std::vector<Int_t> levels;
        for (size_t i = 0; i < PDGLevels.size(); i++)
        {
            if (PDGLevels[i] == "parent") {
                levels.emplace_back(0);
            }
            else if (PDGLevels[i] == "daughters") { 
                levels.emplace_back(1);
            }
            else {
                levels.emplace_back(2);
            }
        }
        fPDGLevels = levels;
    }

    void MCEnergyDeposits::processEvent(
        art::ValidHandle<std::vector<simb::MCParticle>> mcParticles,
        art::ValidHandle<std::vector<sim::SimEnergyDeposit>> mcEnergyDeposits
    )
    {
        MCEdep mcEdep;
        if (mcParticles.isValid() and mcEnergyDeposits.isValid())
        {
            /**
             * We first iterate through all particles and create a map of 
             * parent-daughter pairs for track ids.  This way we can search
             * recursively for ancestors of each particle.
             */
            std::map<Int_t, Int_t> parentDaughterMap;
            std::map<Int_t, Int_t> particlePDGMap;
            for (auto particle : *mcParticles)
            {
                parentDaughterMap[particle.TrackId()] = particle.Mother();
                particlePDGMap[particle.TrackId()] = particle.PdgCode();
            }
            /**
             * Now we work through the energy deposit list, looking first
             * at each pdg code and seeing if it is in our desired list.
             * Then we determine how deep to search through each edep until
             * we find a match.
             */
            
            for (auto energyDeposit : *mcEnergyDeposits)
            {
                // find the track id of the primary ancestor
                Int_t level = 0;
                Int_t track_id = energyDeposit.TrackID();
                Int_t mother = parentDaughterMap[track_id];
                while(mother != 0)
                {
                    level += 1;
                    track_id = mother;
                    mother = parentDaughterMap[mother];
                }
                // see if the ancestors pdg code is in the list
                auto pdg_index = std::find(
                    fPDGCodes.begin(), 
                    fPDGCodes.end(), 
                    particlePDGMap[track_id]
                )
                if (pdg_index != fPDGCodes.end())
                {
                    if (
                        (fPDGLevels[pdg_index] == 0 and level == 0) or
                        (fPDGLevels[pdg_index] == 1 and level != 0) or
                        (fPDGLevels[pdg_index] == 2)
                    )
                    {
                        fMCEdep.pdg.emplace_back(particlePDGMap[track_id]);
                        fMCEdep.track_id.emplace_back(energyDeposit.TrackID());
                        fMCEdep.ancestor_id.emplace_back(track_id);
                        fMCEdep.level.emplace_back(level);
                        fMCEdep.edep_x.emplace_back(energyDeposit.StartX());
                        fMCEdep.edep_y.emplace_back(energyDeposit.StartY());
                        fMCEdep.edep_z.emplace_back(energyDeposit.StartZ());
                        fMCEdep.energy.emplace_back(energyDeposit.Energy());
                        fMCEdep.num_electrons.emplace_back(energyDeposit.NumElectrons());
                    }
                }
            }
        }
        fMCEdep = mcEdep;
        fMCEnergyDepositsTree->Fill();
    }
}