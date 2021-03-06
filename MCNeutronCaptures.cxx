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
        fMCNeutronCapturesTree->Branch("gamma_ids", &fMCNeutronStatistics.gamma_ids);              
        fMCNeutronCapturesTree->Branch("gamma_neutron_ids", &fMCNeutronStatistics.gamma_neutron_ids);       
        fMCNeutronCapturesTree->Branch("gamma_energy", &fMCNeutronStatistics.gamma_energy);         
        fMCNeutronCapturesTree->Branch("gamma_electron_energy", &fMCNeutronStatistics.gamma_electron_energy);
        fMCNeutronCapturesTree->Branch("gamma_edep_energy", &fMCNeutronStatistics.gamma_edep_energy);    

        fMCNeutronCapturesTree->Branch("electron_ids", &fMCNeutronStatistics.electron_ids);            
        fMCNeutronCapturesTree->Branch("electron_parent", &fMCNeutronStatistics.electron_parent);         
        fMCNeutronCapturesTree->Branch("electron_gamma_ids", &fMCNeutronStatistics.electron_gamma_ids);      
        fMCNeutronCapturesTree->Branch("electron_neutron_ids", &fMCNeutronStatistics.electron_neutron_ids);    
        fMCNeutronCapturesTree->Branch("electron_energy", &fMCNeutronStatistics.electron_energy);

        fMCNeutronCapturesTree->Branch("primary", &fMCNeutronStatistics.primary);
        fMCNeutronCapturesTree->Branch("capture", &fMCNeutronStatistics.capture);
        fMCNeutronCapturesTree->Branch("capture_tpc", &fMCNeutronStatistics.capture_tpc);
        fMCNeutronCapturesTree->Branch("capture_tpc_lar", &fMCNeutronStatistics.capture_tpc_lar);
        fMCNeutronCapturesTree->Branch("inelastic", &fMCNeutronStatistics.inelastic);

        fMCNeutronCapturesTree->Branch("total_number_steps", &fMCNeutronStatistics.total_number_steps);
        fMCNeutronCapturesTree->Branch("cryo_number_steps", &fMCNeutronStatistics.cryo_number_steps);
        fMCNeutronCapturesTree->Branch("tpc_number_steps", &fMCNeutronStatistics.tpc_number_steps);
        fMCNeutronCapturesTree->Branch("lar_number_steps", &fMCNeutronStatistics.lar_number_steps);

        fMCNeutronCapturesTree->Branch("entered_tpc", &fMCNeutronStatistics.entered_tpc);
        fMCNeutronCapturesTree->Branch("entered_tpc_step", &fMCNeutronStatistics.entered_tpc_step);
        fMCNeutronCapturesTree->Branch("entered_tpc_time", &fMCNeutronStatistics.entered_tpc_time);
        fMCNeutronCapturesTree->Branch("entered_tpc_energy", &fMCNeutronStatistics.entered_tpc_energy);

        fMCNeutronCapturesTree->Branch("exited_tpc", &fMCNeutronStatistics.exited_tpc);
        fMCNeutronCapturesTree->Branch("exited_tpc_step", &fMCNeutronStatistics.exited_tpc_step);
        fMCNeutronCapturesTree->Branch("exited_tpc_time", &fMCNeutronStatistics.exited_tpc_time);
        fMCNeutronCapturesTree->Branch("exited_tpc_energy", &fMCNeutronStatistics.exited_tpc_energy);

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

    void MCNeutronCaptures::processEvent(
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
        const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>& mcEnergyDeposits
    )
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
                    /**
                     * We want to accumulate statistics about neutron captures, specifically
                     * when the neutron captures on LAr and what the location in the detector
                     * is.
                     * 
                     * There are some subtleties here regarding the bookkeeping that Geant4 does.
                     * For instance, if a neutron undergoes an inelastic scatter (i.e. when
                     * particle.EndProcess() == "neutronInelastic"), then Geant4
                     * creates a new neutron and kills the previous one.  Thus, if we want a "true"
                     * calculation of neutron lifetime, we need to stick together the trajectories
                     * of the parent and daughter neutron, otherwise we will miss (albeit small)
                     * correction.
                     */
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
                    /// set up trajectory information
                    Int_t cryo_number_steps = 0;
                    Int_t tpc_number_steps = 0;
                    Int_t lar_number_steps = 0;

                    Double_t total_distance = 0.0;
                    Double_t cryo_distance = 0.0;
                    Double_t tpc_distance = 0.0;

                    Double_t tpc_avg_material = 0.0;

                    bool entered_tpc = false;
                    Int_t entered_tpc_step = -1;
                    Double_t entered_tpc_time = -1;
                    Double_t entered_tpc_energy = -1;

                    bool exited_tpc = false;
                    Int_t exited_tpc_step = -1;
                    Double_t exited_tpc_time = -1;
                    Double_t exited_tpc_energy = -1;
                    /// iterate over trajectory
                    for (size_t i = 0; i < particle.NumberTrajectoryPoints(); i++)
                    {
                        /// check what volume the step is in
                        DetectorVolume current_volume = fGeometry->getVolume(
                            particle.Vx(i), particle.Vy(i), particle.Vz(i)
                        );
                        /**
                         * We want to gather statistics about when a neutron enters the
                         * active TPC region, if and how long it stays there, and what
                         * sorts of interactions happen inside (i.e. what is the average
                         * density of materials the neutron is interacting with, is it all
                         * LAr?).
                         * 
                         * Knowing these quantities helps us calculate an effective lifetime
                         * for the neutron inside the TPC, which is affected by the materials
                         * present (e.g. if scattering off a denser material than LAr like iron,
                         * then the neutron will recoil with a different expected energy change).
                         */
                        if (current_volume.volume_type == 1) 
                        {
                            if (entered_tpc == true and exited_tpc == false) 
                            {
                                exited_tpc = true;
                                exited_tpc_step = i;
                                exited_tpc_time = particle.T(i);
                                exited_tpc_energy = particle.E(i);
                            }
                            cryo_number_steps += 1;
                        }
                        else if (current_volume.volume_type == 2) 
                        {
                            if (entered_tpc == false and exited_tpc == false) 
                            {
                                entered_tpc = true;
                                entered_tpc_step = i;
                                entered_tpc_time = particle.T(i);
                                entered_tpc_energy = particle.E(i);
                            }
                            tpc_number_steps += 1;
                            tpc_avg_material += current_volume.material;
                            if (current_volume.material_name == "LAr") 
                            {
                                lar_number_steps += 1;
                            }
                        }
                        else
                        {
                            if (entered_tpc == true and exited_tpc == false) 
                            {
                                exited_tpc = true;
                                exited_tpc_step = i;
                                exited_tpc_time = particle.T(i);
                                exited_tpc_energy = particle.E(i);
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
                    neutronStatistics.entered_tpc.emplace_back(entered_tpc);
                    neutronStatistics.entered_tpc_step.emplace_back(entered_tpc_step);
                    neutronStatistics.entered_tpc_time.emplace_back(entered_tpc_time);
                    neutronStatistics.entered_tpc_energy.emplace_back(entered_tpc_energy);

                    neutronStatistics.exited_tpc.emplace_back(exited_tpc);
                    neutronStatistics.exited_tpc_step.emplace_back(exited_tpc_step);
                    neutronStatistics.exited_tpc_time.emplace_back(exited_tpc_time);
                    neutronStatistics.exited_tpc_energy.emplace_back(exited_tpc_energy);

                    if (tpc_number_steps > 0) {
                        tpc_avg_material /= tpc_number_steps;
                    }

                    neutronStatistics.cryo_number_steps.emplace_back(cryo_number_steps);
                    neutronStatistics.tpc_number_steps.emplace_back(tpc_number_steps);
                    neutronStatistics.lar_number_steps.emplace_back(lar_number_steps);

                    neutronStatistics.tpc_avg_material.emplace_back(tpc_avg_material);

                    neutronStatistics.total_distance.emplace_back(total_distance);
                    neutronStatistics.cryo_distance.emplace_back(cryo_distance);
                    neutronStatistics.tpc_distance.emplace_back(tpc_distance);
                    /**
                     * TODO:
                     * get lifetimes
                     * accumulate trajectories
                     */
                }
                // check if the particle is a gamma
                if (particle.PdgCode() == 22)
                {
                    for(size_t i = 0; i < neutronStatistics.neutron_ids.size(); i++)
                    {
                        if (neutronStatistics.neutron_ids[i] == particle.Mother())
                        {
                            neutronStatistics.gamma_ids.emplace_back(particle.TrackId());
                            neutronStatistics.gamma_neutron_ids.emplace_back(particle.Mother());
                            neutronStatistics.gamma_energy.emplace_back(particle.E());
                            neutronStatistics.gamma_electron_energy.emplace_back(0);
                            neutronStatistics.gamma_edep_energy.emplace_back(0);
                        }
                    }
                }
                // check if the particle is an electron
                if (particle.PdgCode() == 11)
                {
                    // check for gammas first
                    for (size_t i = 0; i < neutronStatistics.gamma_ids.size(); i++)
                    {
                        if (neutronStatistics.gamma_ids[i] == particle.Mother())
                        {
                            neutronStatistics.electron_ids.emplace_back(particle.TrackId());
                            neutronStatistics.electron_parent.emplace_back(particle.Mother());
                            neutronStatistics.electron_gamma_ids.emplace_back(particle.Mother());
                            neutronStatistics.gamma_electron_energy[i] += particle.E();
                            // find the corresponding neutron id
                            for(size_t j = 0; j < neutronStatistics.neutron_ids.size(); j++)
                            {
                                if (neutronStatistics.neutron_ids[j] == neutronStatistics.gamma_neutron_ids[i])
                                {
                                    neutronStatistics.electron_neutron_ids.emplace_back(neutronStatistics.neutron_ids[j]);
                                }
                            }
                            neutronStatistics.electron_energy.emplace_back(particle.E());
                        }
                    }
                    // then check electrons
                    for (size_t i = 0; i < neutronStatistics.electron_ids.size(); i++)
                    {
                        if (neutronStatistics.electron_ids[i] == particle.Mother())
                        {
                            neutronStatistics.electron_ids.emplace_back(particle.TrackId());
                            neutronStatistics.electron_parent.emplace_back(particle.Mother());
                            // find the corresponding gamma
                            neutronStatistics.electron_gamma_ids.emplace_back(neutronStatistics.electron_gamma_ids[i]);
                            neutronStatistics.electron_neutron_ids.emplace_back(neutronStatistics.electron_neutron_ids[i]);
                            neutronStatistics.electron_energy.emplace_back(particle.E());
                        }
                    }
                }
            }
        }
        if (mcEnergyDeposits.isValid())
        {
            for (auto energyDeposit : *mcEnergyDeposits)
            {
                // check the list of electrons
                for (size_t i = 0; i < neutronStatistics.electron_ids.size(); i++)
                {
                    if (neutronStatistics.electron_ids[i] == energyDeposit.TrackID())
                    {
                        neutronStatistics.neutron_edep_ids.emplace_back(neutronStatistics.electron_neutron_ids[i]);
                        neutronStatistics.neutron_edep_parent.emplace_back(energyDeposit.TrackID());
                        neutronStatistics.neutron_edep_gamma_ids.emplace_back(neutronStatistics.electron_gamma_ids[i]);
                        for (size_t j = 0; j < neutronStatistics.gamma_ids.size(); j++)
                        {
                            if (neutronStatistics.gamma_ids[j] == neutronStatistics.electron_gamma_ids[i])
                            {
                                neutronStatistics.gamma_edep_energy[j] += energyDeposit.Energy();
                            }
                        }
                        neutronStatistics.neutron_edep_energy.emplace_back(energyDeposit.Energy());
                        neutronStatistics.neutron_edep_num_electrons.emplace_back(energyDeposit.NumElectrons());
                        neutronStatistics.neutron_edep_x.emplace_back(energyDeposit.StartX());
                        neutronStatistics.neutron_edep_y.emplace_back(energyDeposit.StartY());
                        neutronStatistics.neutron_edep_z.emplace_back(energyDeposit.StartZ());
                        break;
                    }
                }
            }
        }
        fMCNeutronStatistics = neutronStatistics;
        fMCNeutronCapturesTree->Fill();
    }
}