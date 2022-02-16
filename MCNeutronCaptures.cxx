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
        fMCNeutronCapturesTree->Branch("capture", &fMCNeutronStatistics.capture);
        fMCNeutronCapturesTree->Branch("capture_tpc", &fMCNeutronStatistics.capture_tpc);
        fMCNeutronCapturesTree->Branch("capture_tpc_lar", &fMCNeutronStatistics.capture_tpc_lar);
        fMCNeutronCapturesTree->Branch("inelastic", &fMCNeutronStatistics.inelastic);

        fMCNeutronCapturesTree->Branch("total_number_steps", &fMCNeutronStatistics.total_number_steps);
        fMCNeutronCapturesTree->Branch("cryo_number_steps", &fMCNeutronStatistics.cryo_number_steps);
        fMCNeutronCapturesTree->Branch("tpc_number_steps", &fMCNeutronStatistics.tpc_number_steps);
        fMCNeutronCapturesTree->Branch("lar_number_steps", &fMCNeutronStatistics.lar_number_steps);

        fMCNeutronCapturesTree->Branch("tpc_avg_material", &fMCNeutronStatistics.tpc_avg_material);

        fMCNeutronCapturesTree->Branch("total_distance", &fMCNeutronStatistics.total_distance);
        fMCNeutronCapturesTree->Branch("cryo_distance", &fMCNeutronStatistics.cryo_distance);
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
                        neutronStatistics.capture.emplace_back(true);
                        // get the ending volume
                        DetectorVolume ending_volume = fGeometry->getVolume(
                            particle.EndX(), particle.EndY(), particle.EndZ()
                        );
                        if (ending_volume.volume_type == 2) {
                            neutronStatistics.capture_tpc.emplace_back(true);
                        }
                        else {
                            neutronStatistics.capture_tpc.emplace_back(false);
                        }
                        if (ending_volume.material_name == "LAr") {
                            neutronStatistics.capture_tpc_lar.emplace_back(true);
                        }
                        else {
                            neutronStatistics.capture_tpc_lar.emplace_back(false);
                        }
                        neutronStatistics.inelastic.emplace_back(false);
                        neutronStatistics.neutron_capture_x.emplace_back(particle.EndX());
                        neutronStatistics.neutron_capture_y.emplace_back(particle.EndY());
                        neutronStatistics.neutron_capture_z.emplace_back(particle.EndZ());
                    }
                    else
                    {
                        neutronStatistics.capture.emplace_back(false);
                        neutronStatistics.capture_tpc.emplace_back(false);
                        neutronStatistics.capture_tpc_lar.emplace_back(false);
                        if (particle.EndProcess() == "neutronInelastic") {
                            neutronStatistics.inelastic.emplace_back(true);
                        }
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

                    Double_t tpc_avg_material = 0.0;
                    for (size_t i = 0; i < particle.NumberTrajectoryPoints(); i++)
                    {
                        /// check what volume the step is in
                        DetectorVolume current_volume = fGeometry->getVolume(
                            particle.Vx(i), particle.Vy(i), particle.Vz(i)
                        );
                        if (current_volume.volume_type == 1) {
                            cryo_number_steps += 1;
                        }
                        else if (current_volume.volume_type == 2) 
                        {
                            tpc_number_steps += 1;
                            tpc_avg_material += current_volume.material;
                            if (current_volume.material_name == "LAr") 
                            {
                                lar_number_steps += 1;
                            }
                        }
                        /// compute distances and average material interactions
                        if (i > 0)
                        {
                            DetectorVolume previous_volume = fGeometry->getVolume(
                                particle.Vx(i-1), particle.Vy(i-1), particle.Vz(i-1)
                            );
                            Double_t step_distance = euclidean_distance(
                                particle.Vx(i), particle.Vy(i), particle.Vz(i),
                                particle.Vx(i-1), particle.Vy(i-1), particle.Vz(i-1)
                            );
                            total_distance += step_distance;
                            /**
                             * Check wether the last and current step are in the same
                             * volume type and if so, add that distance to that type,
                             * otherwise, divide the distance in half for the previous
                             * and current steps (basically assume an average distance).
                             */
                            // Cryostat (volume_type == 1)
                            if (current_volume.volume_type == 1)
                            {
                                if (previous_volume.volume_type == 1) {
                                    cryo_distance += step_distance;
                                }
                                else 
                                {
                                    cryo_distance += step_distance / 2.0;
                                    if (previous_volume.volume_type == 2) {
                                        tpc_distance += step_distance / 2.0;
                                    }
                                }
                            }
                            // TPC (volume_type == 2)
                            else if (current_volume.volume_type == 2)
                            {
                                if (previous_volume.volume_type == 2) {
                                    tpc_distance += step_distance;
                                }
                                else
                                {
                                    tpc_distance += step_distance / 2.0;
                                    if (previous_volume.volume_type == 1) {
                                        cryo_distance += step_distance / 2.0;
                                    }
                                }
                            }
                        }
                    }
                    // accumulate step results
                    tpc_avg_material /= tpc_number_steps;

                    neutronStatistics.cryo_number_steps.emplace_back(cryo_number_steps);
                    neutronStatistics.tpc_number_steps.emplace_back(tpc_number_steps);
                    neutronStatistics.lar_number_steps.emplace_back(lar_number_steps);

                    neutronStatistics.tpc_avg_material.emplace_back(tpc_avg_material);

                    neutronStatistics.total_distance.emplace_back(total_distance);
                    neutronStatistics.cryo_distance.emplace_back(cryo_distance);
                    neutronStatistics.tpc_distance.emplace_back(tpc_distance);
                }
            }
        }
        fMCNeutronStatistics = neutronStatistics;
        fMCNeutronCapturesTree->Fill();
    }
}