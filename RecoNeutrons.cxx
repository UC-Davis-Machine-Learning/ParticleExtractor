/**
 * @file RecoNeutrons.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-02-17
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "RecoNeutrons.h"

namespace extractor
{
    RecoNeutrons::RecoNeutrons()
    {
        fRecoNeutronsTree = fTFileService->make<TTree>("reco_neutrons", "reco_neutrons");
        fRecoNeutronsTree->Branch("sp_x", &fRecoNeutronsSet.sp_x);
        fRecoNeutronsTree->Branch("sp_y", &fRecoNeutronsSet.sp_y);
        fRecoNeutronsTree->Branch("sp_z", &fRecoNeutronsSet.sp_z);
        fRecoNeutronsTree->Branch("neutron_id", &fRecoNeutronsSet.neutron_id);
        fRecoNeutronsTree->Branch("gamma_id", &fRecoNeutronsSet.gamma_id);
        fRecoNeutronsTree->Branch("gamma_energy", &fRecoNeutronsSet.gamma_energy);
        fRecoNeutronsTree->Branch("summed_adc", &fRecoNeutronsSet.summed_adc);
        fRecoNeutronsTree->Branch("mean_adc", &fRecoNeutronsSet.mean_adc);
        fRecoNeutronsTree->Branch("peak_adc", &fRecoNeutronsSet.peak_adc);
        fRecoNeutronsTree->Branch("sigma_adc", &fRecoNeutronsSet.sigma_adc);
    }

    RecoNeutrons::~RecoNeutrons()
    {}

    void RecoNeutrons::setBoundingBoxType(std::string volumeType)
    {
        if (volumeType == "TPC" or volumeType == "tpc") { 
            fBoundingBoxType = VolumeType::TPC;
        }
        else if (volumeType == "Cryo" or volumeType == "cryo") {
            fBoundingBoxType = VolumeType::Cryostat;
        }
        else {
            fBoundingBoxType = VolumeType::World;
        }
    }

    void RecoNeutrons::processEvent(
        detinfo::DetectorClocksData const& clockData,
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
        const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
        const art::ValidHandle<std::vector<recob::SpacePoint>>& recoSpacePoints,
        const art::FindManyP<recob::Hit>& hitPandoraSPsAssn //to associate space points from pandora to hits
    )
    {
        RecoNeutronsSet RecoNeutronsSet;
        if (mcParticles.isValid() and mcChannels.isValid() and recoSpacePoints.isValid())
        {
            /**
             * We first iterate through all particles and create a map of 
             * parent-daughter pairs for track ids.  This way we can search
             * recursively for ancestors of each particle.
             */
            std::map<Int_t, Int_t> parentDaughterMap;
            std::map<Int_t, Int_t> particlePDGMap;

            std::vector<int> neutron_captures;
            std::vector<std::vector<int>> gamma_ids;
            std::vector<std::vector<double>> gamma_energy;

            std::map<Int_t, Int_t> neutronMap;
            std::map<Int_t, Int_t> gammaMap;
            std::map<Int_t, bool> neutronCapture;

            for (auto particle : *mcParticles)
            {
                parentDaughterMap[particle.TrackId()] = particle.Mother();
                particlePDGMap[particle.TrackId()] = particle.PdgCode();
            }
            for (auto particle : *mcParticles)
            {
                // check if the particle is a neutron
                if (particle.PdgCode() == 2112)
                {
                    if (particle.EndProcess() == "nCapture")
                    {
                        neutron_captures.emplace_back(particle.TrackId());
                        gamma_ids.emplace_back(std::vector<int>());
                        gamma_energy.emplace_back(std::vector<double>());
                    }
                }
                // check if the particle is a gamma
                else if (particle.PdgCode() == 22)
                {
                    for(size_t i = 0; i < neutron_captures.size(); i++)
                    {
                        if (neutron_captures[i] == particle.Mother())
                        {
                            gamma_ids[i].emplace_back(particle.TrackId());
                            gamma_energy[i].emplace_back(particle.E());
                        }
                    }
                }
                // otherwise see if the particle is from a neutron capture
                else
                {
                    Int_t mother = particle.Mother();
                    Int_t track_id = particle.TrackId();
                    while (mother != 0)
                    {
                        for(size_t i = 0; i < neutron_captures.size(); i++)
                        {
                            if (neutron_captures[i] == mother)
                            {
                                std::cout << "Found capture: " << neutron_captures[i] << " for track id: " << particle.TrackId() << std::endl;
                                neutronMap[particle.TrackId()] = i;
                                for (size_t j = 0; j < gamma_ids[i].size(); j++)
                                {
                                    if (gamma_ids[i][j] == track_id)
                                    {
                                        std::cout << "Found gamma: " << gamma_ids[i][j] << " for track id: " << particle.TrackId() << std::endl;
                                        gammaMap[particle.TrackId()] = j;
                                        neutronCapture[particle.TrackId()] = true;
                                        break;
                                    }
                                }
                                break;
                            }
                        }
                        track_id = mother;
                        mother = parentDaughterMap[track_id];
                    }
                }
            }
            
            std::vector<art::Ptr<recob::SpacePoint>> pointsList;
            art::fill_ptr_vector(pointsList, recoSpacePoints);            
            for (size_t i = 0; i < pointsList.size(); i++)
            {
                auto& spsHit = hitPandoraSPsAssn.at(i);
                for (auto hit : spsHit)
                {  
                    Int_t track_id = TruthMatchUtils::TrueParticleID(
                        clockData, hit, false
                    );
                    // check that track_id is present in parentDaughterMap
                    if (parentDaughterMap.find(track_id) == parentDaughterMap.end())
                    {
                        std::cout << "Track ID: " << track_id << " does not have an associated mother!" << std::endl;
                        continue;
                    }
                    if (particlePDGMap.find(track_id) == particlePDGMap.end())
                    {
                        std::cout << "Track ID: " << track_id << " does not have an associated pdg!" << std::endl;
                        continue;
                    }
                    if (neutronCapture[track_id])
                    {
                        std::cout << "Found hit: " << track_id << " with capture: ";
                        std::cout << neutron_captures[neutronMap[track_id]] << " and gamma: ";
                        std::cout << gamma_ids[neutronMap[track_id]][gammaMap[track_id]] << std::endl;
                        // collect results
                        auto xyz = pointsList[i]->XYZ();
                        // check if point is in active volume
                        // Determine if edep is within the desired volume
                        DetectorVolume edep_volume = fGeometry->getVolume(
                            xyz[0], xyz[1], xyz[2]
                        );
                        if (edep_volume.volume_type != fBoundingBoxType) {
                            continue;
                        }
                        
                        RecoNeutronsSet.sp_x.emplace_back(xyz[0]);
                        RecoNeutronsSet.sp_y.emplace_back(xyz[1]);
                        RecoNeutronsSet.sp_z.emplace_back(xyz[2]);
                        RecoNeutronsSet.neutron_id.emplace_back(neutron_captures[neutronMap[track_id]]);
                        RecoNeutronsSet.gamma_id.emplace_back(gamma_ids[neutronMap[track_id]][gammaMap[track_id]]);
                        RecoNeutronsSet.gamma_energy.emplace_back(gamma_energy[neutronMap[track_id]][gammaMap[track_id]]);
                        RecoNeutronsSet.summed_adc.emplace_back(hit->SummedADC());
                        RecoNeutronsSet.mean_adc.emplace_back(hit->PeakTime());
                        RecoNeutronsSet.peak_adc.emplace_back(hit->PeakAmplitude());
                        RecoNeutronsSet.sigma_adc.emplace_back(hit->RMS());
                    }
                }
            }
        }
        fRecoNeutronsSet = RecoNeutronsSet;
        fRecoNeutronsTree->Fill();
    }
}