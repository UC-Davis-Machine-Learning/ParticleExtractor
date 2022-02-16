/**
 * @file MCNeutronCaptures.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-02-15
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "MCNeutronCaptures.h"

namespace extractor
{
    MCNeutronCaptures::MCNeutronCaptures()
    {
        fMCNeutronCapturesTree = fTFileService->make<TTree>("mc_neutron_captures", "mc_neutron_captures");
        fMCNeutronCapturesTree->Branch("neutron_ids", &fMCNeutronStatistics.neutron_ids);
        fMCNeutronCapturesTree->Branch("total_number_steps", &fMCNeutronStatistics.total_number_steps);
        fMCNeutronCapturesTree->Branch("neutron_capture_x", &fMCNeutronStatistics.neutron_capture_x);
        fMCNeutronCapturesTree->Branch("neutron_capture_y", &fMCNeutronStatistics.neutron_capture_y);
        fMCNeutronCapturesTree->Branch("neutron_capture_z", &fMCNeutronStatistics.neutron_capture_z);
    }

    MCNeutronCaptures::~MCNeutronCaptures()
    {}

    void MCNeutronCaptures::processEvent(art::ValidHandle<std::vector<simb::MCParticle>> mcParticles)
    {
        if (mcParticles.isValid())
        {
            NeutronStatistics neutronStatistics;
            for (auto particle : *mcParticles)
            {
                // check if the particle is a neutron
                if (particle.PdgCode() == 2112)
                {
                    if (particle.Mother() == 0)
                    {
                        neutronStatistics.neutron_ids.emplace_back(particle.TrackId());
                        neutronStatistics.total_number_steps.emplace_back(particle.NumberTrajectoryPoints());

                        if (particle.EndProcess() == "nCapture")
                        {
                            neutronStatistics.neutron_capture_x.emplace_back(particle.EndX());
                            neutronStatistics.neutron_capture_y.emplace_back(particle.EndY());
                            neutronStatistics.neutron_capture_z.emplace_back(particle.EndZ());
                        }
                    }
                }
            }
        }
        fMCNeutronStatistics = neutronStatistics;
        fMCNeutronCapturesTree->Fill();
    }
}