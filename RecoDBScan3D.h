/**
 * @file RecoDBScan3D.h
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-04-29
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLine.h"
#include "TAxis.h"
#include "TTimeStamp.h"

#include <vector>
#include <fstream>
#include "TPaveStats.h"
#include <iostream>
#include <string>
#include "math.h"
#include "stdio.h"
#include <iterator>
#include <map>
#include <list>

#include "larsim/Utils/TruthMatchUtils.h"

#include "DetectorGeometry.h"

namespace detinfo {
  class DetectorClocksData;
}

struct sptStruct {
    float posX;
    float posY;
    float posZ;

    bool operator==(const sptStruct &other) const{ 
        return (posX == other.posX
        && posY == other.posY
        && posZ == other.posZ);
    }
};

namespace std{

    template <>
    struct hash<sptStruct>
    {
        std::size_t operator()(const sptStruct& k) const
        {
            using std::size_t;
            using std::hash;
            //using std::float;

            // Compute individual hash values for first,
            // second and third and combine them using XOR
            // and bit shifting:

            return ((hash<float>()(k.posX)
            ^ (hash<float>()(k.posY) << 1)) >> 1)
            ^ (hash<float>()(k.posZ) << 1);
        }
    };
}

namespace extractor
{
    struct RecoDBScanSP
    {
        std::vector<Int_t> sp_cluster_size;
        std::vector<Int_t> sp_cluster_num; //0 = non-cluster spt; > 0 -> cluster spt
        std::vector<Int_t> sp_pdg;
        std::vector<Double_t> sp_x;
        std::vector<Double_t> sp_y;
        std::vector<Double_t> sp_z;
        std::vector<Double_t> summed_adc;
    };

    struct SlcNumSizeStruct {
        int SlcNum;
        int SlcSize;
    };

    class RecoDBScan3D
    {
    public:
        RecoDBScan3D();
        ~RecoDBScan3D();

        void setPDGCodes(std::vector<Int_t> PDGCodes) { fPDGCodes = PDGCodes; }
        void setPDGLevels(std::vector<Int_t> PDGLevels) { fPDGLevels = PDGLevels; }
        void setPDGLevels(std::vector<std::string> PDGLevels);

        void setBoundingBoxType(std::string volumeType);

        void processEvent(
            detinfo::DetectorClocksData const& clockData,
            const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
            const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
            const art::ValidHandle<std::vector<recob::SpacePoint>>& recoSpacePoints,
            const art::ValidHandle<std::vector<recob::Slice>>& recoSlices,
            const art::FindManyP<recob::Hit>& hitSpacePointAssn, //to associate space points to hits
            const art::FindManyP<recob::SpacePoint>& spacePointSliceAssn //to associate slices created by DBScan to space points
        );

        RecoDBScanSP getRecoDBScanSP() const { return fRecoDBScanSP; }

    private:
        art::ServiceHandle<art::TFileService> fTFileService;
        TTree *fRecoDBScan3DTree;
        // geometry information
        DetectorGeometry* fGeometry = DetectorGeometry::getInstance("RecoDBScan3D");

        // pdg codes to construct
        VolumeType fBoundingBoxType;
        std::vector<Int_t> fPDGCodes;
        std::vector<Int_t> fPDGLevels;

        // struct for holding event information
        RecoDBScanSP fRecoDBScanSP;

    };
}
