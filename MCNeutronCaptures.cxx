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
        fMCNeutronCapturesTree->Branch("primary", &fMCNeutronStatistics.primary);

        fMCNeutronCapturesTree->Branch("total_number_steps", &fMCNeutronStatistics.total_number_steps);
        fMCNeutronCapturesTree->Branch("tpc_number_steps", &fMCNeutronStatistics.tpc_number_steps);
        fMCNeutronCapturesTree->Branch("lar_number_steps", &fMCNeutronStatistics.lar_number_steps);

        fMCNeutronCapturesTree->Branch("total_distance", &fMCNeutronStatistics.total_distance);
        fMCNeutronCapturesTree->Branch("tpc_distance", &fMCNeutronStatistics.tpc_distance);

        fMCNeutronCapturesTree->Branch("neutron_capture_x", &fMCNeutronStatistics.neutron_capture_x);
        fMCNeutronCapturesTree->Branch("neutron_capture_y", &fMCNeutronStatistics.neutron_capture_y);
        fMCNeutronCapturesTree->Branch("neutron_capture_z", &fMCNeutronStatistics.neutron_capture_z);
    }

    MCNeutronCaptures::~MCNeutronCaptures()
    {}

    void MCNeutronCaptures::processEvent(art::ValidHandle<std::vector<simb::MCParticle>> mcParticles)
    {
        MCNeutronStatistics neutronStatistics;
        if (mcParticles.isValid())
        {
            for (auto particle : *mcParticles)
            {
                // check if the particle is a neutron
                if (particle.PdgCode() == 2112)
                {
                    neutronStatistics.neutron_ids.emplace_back(particle.TrackId());
                    neutronStatistics.total_number_steps.emplace_back(particle.NumberTrajectoryPoints());
                    /// wether neutron is a primary
                    if (particle.Mother() == 0) {
                        neutronStatistics.primary.emplace_back(true);
                    }
                    else {
                        neutronStatistics.primary.emplace_back(false);
                    }
                    /// get capture locations
                    if (particle.EndProcess() == "nCapture")
                    {
                        neutronStatistics.neutron_capture_x.emplace_back(particle.EndX());
                        neutronStatistics.neutron_capture_y.emplace_back(particle.EndY());
                        neutronStatistics.neutron_capture_z.emplace_back(particle.EndZ());
                    }
                    else
                    {
                        neutronStatistics.neutron_capture_x.emplace_back(-1e9);
                        neutronStatistics.neutron_capture_y.emplace_back(-1e9);
                        neutronStatistics.neutron_capture_z.emplace_back(-1e9);
                    }
                    /// iterate over trajectory
                    Int_t cryo_number_steps = 0;
                    Int_t tpc_number_steps = 0;
                    Int_t lar_number_steps = 0;
                    Double_t total_distance = 0.0;
                    Double_t cryo_distance = 0.0;
                    Double_t tpc_distance = 0.0;
                    for (size_t i = 0; i < particle.NumberTrajectoryPoints(); i++)
                    {
                        /// check what volume the step is in
                        DetectorVolume volume = fGeometry->getVolume(
                            particle.Vx(i), particle.Vy(i), particle.Vz(i)
                        );
                        std::cout << volume.volume_type << "," << volume.volume_name << "," << volume.material_name << "," << volume.material << std::endl;
                    }
                }
            }
        }
        fMCNeutronStatistics = neutronStatistics;
        fMCNeutronCapturesTree->Fill();
    }
}