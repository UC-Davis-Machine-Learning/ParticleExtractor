/**
 * @file RecoDBScan3D.cxx
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-04-29
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "RecoDBScan3D.h"

namespace extractor
{
    RecoDBScan3D::RecoDBScan3D()
    {
        fRecoDBScan3DTree = fTFileService->make<TTree>("reco_DBScan3D", "reco_DBScan3D");
        fRecoDBScan3DTree->Branch("sp_cluster_size", &fRecoDBScanSP.sp_cluster_size);
        fRecoDBScan3DTree->Branch("sp_cluster_num", &fRecoDBScanSP.sp_cluster_num);
        fRecoDBScan3DTree->Branch("sp_pdg", &fRecoDBScanSP.sp_pdg);
        fRecoDBScan3DTree->Branch("sp_x", &fRecoDBScanSP.sp_x);
        fRecoDBScan3DTree->Branch("sp_y", &fRecoDBScanSP.sp_y);
        fRecoDBScan3DTree->Branch("sp_z", &fRecoDBScanSP.sp_z);
        fRecoDBScan3DTree->Branch("summed_adc", &fRecoDBScanSP.summed_adc);
    }

    RecoDBScan3D::~RecoDBScan3D()
    {}

    void RecoDBScan3D::setBoundingBoxType(std::string volumeType)
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

    void RecoDBScan3D::processEvent(
        detinfo::DetectorClocksData const& clockData,
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
        const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
        const art::ValidHandle<std::vector<recob::SpacePoint>>& recoSpacePoints,
        const art::ValidHandle<std::vector<recob::Slice>>& recoSlices,
        const art::FindManyP<recob::Hit>& hitSpacePointAssn, //to associate space points to hits
        const art::FindManyP<recob::SpacePoint>& spacePointSliceAssn //to associate slices created by DBScan to space points
    )
    {
        RecoDBScanSP recoDBScanSP;
        if (mcParticles.isValid() and mcChannels.isValid() and recoSpacePoints.isValid() and recoSlices.isValid())
        {
            /**
             * We first iterate through all particles and create a map of 
             * parent-daughter pairs for track ids.  This way we can search
             * recursively for ancestors of each particle.
             */
            std::map<Int_t, Int_t> particlePDGMap;
            for (auto particle : *mcParticles)
            {
                particlePDGMap[particle.TrackId()] = particle.PdgCode();
            }

            // Spacepoints
            std::vector<art::Ptr<recob::SpacePoint>> pointsList;
            art::fill_ptr_vector(pointsList, recoSpacePoints);
            // Slices
            std::vector<art::Ptr<recob::Slice>> sliceList;
            art::fill_ptr_vector(sliceList, recoSlices);

            // Storing spacepoints and their parent slice info
            std::unordered_map<sptStruct, SlcNumSizeStruct> SptClusterMap;
            for (int i = 0; i < (int) sliceList.size(); ++i) {

                auto & sliceSps = spacePointSliceAssn.at(i);

                SlcNumSizeStruct SlcInfo;
                SlcInfo.SlcNum = i+1;
                SlcInfo.SlcSize = (int) sliceSps.size();

                for (int j = 0; j < (int) sliceSps.size(); ++j){
                    const double* xyz = sliceSps[j]->XYZ();

                    sptStruct spt;
                    spt.posX = xyz[0];
                    spt.posY = xyz[1];
                    spt.posZ = xyz[2];

                    SptClusterMap.emplace(spt, SlcInfo); //sliceSps.size()
                }
            }

            // Looping over spacepoints
            for (int i = 0; i < (int) pointsList.size(); ++i){

                const double* xyz = pointsList[i]->XYZ();

                // check if point is in active volume
                // Determine if edep is within the desired volume
                DetectorVolume edep_volume = fGeometry->getVolume(
                    xyz[0], xyz[1], xyz[2]
                );
                if (edep_volume.volume_type != fBoundingBoxType) {
                    continue;
                }

                sptStruct spt;
                spt.posX = xyz[0];
                spt.posY = xyz[1];
                spt.posZ = xyz[2];

                std::unordered_map<sptStruct, SlcNumSizeStruct>::iterator itr = SptClusterMap.find(spt);
                if (itr != SptClusterMap.end())
                {
                    recoDBScanSP.sp_cluster_num.emplace_back(itr->second.SlcNum); //Cluster space point
                    recoDBScanSP.sp_cluster_size.emplace_back(itr->second.SlcSize);
                } else {
                    recoDBScanSP.sp_cluster_num.emplace_back(0); //Non-Cluster space point
                    recoDBScanSP.sp_cluster_size.emplace_back(0);
                }

                recoDBScanSP.sp_x.emplace_back(xyz[0]);
                recoDBScanSP.sp_x.emplace_back(xyz[1]);
                recoDBScanSP.sp_x.emplace_back(xyz[2]);

                std::vector<Int_t> temp_pdg;
                std::vector<Double_t> temp_summed_adc;

                auto& spsHit = hitSpacePointAssn.at(i);
                for (auto hit : spsHit)
                {
                    Int_t track_id = TruthMatchUtils::TrueParticleID(
                        clockData, hit, false
                    );

                    Int_t pdg = particlePDGMap[track_id];
                    temp_pdg.emplace_back(pdg);
                    temp_summed_adc.emplace_back(hit->SummedADC());
                }

                if (temp_pdg.size() == 0) 
                {
                    recoDBScanSP.sp_pdg.emplace_back(-999);
                    recoDBScanSP.summed_adc.emplace_back(-999);
                }
                else 
                {
                    // now decide how to assign labels
                    // take the value with the largest summed adc
                    int max_index = std::distance(
                        temp_summed_adc.begin(), 
                        std::max_element(temp_summed_adc.begin(), temp_summed_adc.end())
                    );
                    
                    recoDBScanSP.sp_pdg.emplace_back(temp_pdg[max_index]);
                    recoDBScanSP.summed_adc.emplace_back(temp_summed_adc[max_index]);
                }
            }
        }
        fRecoDBScanSP = recoDBScanSP;
        fRecoDBScan3DTree->Fill();
    }
}