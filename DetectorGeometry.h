/**
 * @file    GeometryExtractor.h
 * @brief   A class for extracting geometry information from art root files.
 * @ingroup GeometryExtractor
 * @author  Nicholas Carrara (nmcarrara@ucdavis.edu),
**/
#pragma once

// art includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// special utility includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "art_root_io/TFileService.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// necessary ROOT libraries
#include <TTree.h>
#include <TH1.h>
#include "TH1F.h"
#include "TGeoMaterial.h"
#include "TGeoElement.h"

// std includes
#include <math.h>
#include <string>
#include <vector>
#include <memory>

namespace extractor 
{
    inline Double_t euclidean_distance(
        Double_t x1, Double_t y1, Double_t z1,
        Double_t x2, Double_t y2, Double_t z2
    )
    {
        return sqrt(
            (x1 - x2)*(x1 - x2) + 
            (y1 - y2)*(y1 - y2) + 
            (z1 - z2)*(z1 - z2)
        );
    }
    // list of materials in the detector
    enum MaterialList 
    {

    };
    // list of effective atomic numbers for the materials

    enum VolumeType 
    {
        World,
        Cryostat,
        TPC,
    };
    ///////////////////////////////////////////////////////////////////////////////////////
    // struct for detector volume information
    ///////////////////////////////////////////////////////////////////////////////////////
    struct DetectorVolume
    {
        VolumeType volume_type;
        std::string volume_name;
        std::string material_name;
        double material;
        DetectorVolume() {}
        DetectorVolume(VolumeType volumeType, std::string volumeName, 
            std::string materialName, double material)
        : volume_type(volumeType)
        , volume_name(volumeName)
        , material_name(materialName)
        , material(material)
        {}
    };
    ///////////////////////////////////////////////////////////////////////////////////////
    // struct for bounding boxes
    ///////////////////////////////////////////////////////////////////////////////////////
    struct BoundingBox
    {
        double x_min = 0; double x_max = 0;
        double y_min = 0; double y_max = 0;
        double z_min = 0; double z_max = 0;

        double width()  { return x_max - x_min; }
        double height() { return y_max - y_min; }
        double length() { return z_max - z_min; }

        void setBox(geo::BoxBoundedGeo const& Box) {
            x_min = Box.MinX(); x_max = Box.MaxX();
            y_min = Box.MinY(); y_max = Box.MaxY();
            z_min = Box.MinZ(); z_max = Box.MaxZ();
        }
        void setBox(double xmin, double xmax,
                    double ymin, double ymax,
                    double zmin, double zmax)
        {
            x_min = xmin; x_max = xmax;
            y_min = ymin; y_max = ymax;
            z_min = zmin; z_max = zmax;
        }
        BoundingBox() {}
        BoundingBox(double xs[2], double ys[2], double zs[2])
        {
            x_min = xs[0]; x_max = xs[1];
            y_min = ys[0]; y_max = ys[1];
            z_min = zs[0]; z_max = zs[1];
        }
        BoundingBox(double vals[6])
        {
            x_min = vals[0]; x_max = vals[1];
            y_min = vals[2]; y_max = vals[3];
            z_min = vals[4]; z_max = vals[4];
        }
        BoundingBox(double xmin, double xmax,
                    double ymin, double ymax,
                    double zmin, double zmax)
        {
            x_min = xmin; x_max = xmax;
            y_min = ymin; y_max = ymax;
            z_min = zmin; z_max = zmax;
        }
        BoundingBox(geo::BoxBoundedGeo const& Box) {
            x_min = Box.MinX(); x_max = Box.MaxX();
            y_min = Box.MinY(); y_max = Box.MaxY();
            z_min = Box.MinZ(); z_max = Box.MaxZ();
        }
    };
    ///////////////////////////////////////////////////////////////////////////////////////
    // class for storting geometry information
    ///////////////////////////////////////////////////////////////////////////////////////
    class DetectorGeometry
    {
    private:
        static DetectorGeometry * sInstance;
        static std::mutex sMutex;

    protected:
        DetectorGeometry(const std::string name);
        ~DetectorGeometry() {}
        std::string sName;
        
    public:
        // this singleton cannot be cloned
        DetectorGeometry(DetectorGeometry &other) = delete;
        // singleton should also not be assignable
        void operator=(const DetectorGeometry &) = delete;

        // static method that controls access to the singleton
        // instance
        static DetectorGeometry *getInstance(const std::string& name);

        // getters
        std::string Name() const        { return sName; }
        std::string GetWorldName()      { return fWorldName; }
        BoundingBox GetWorldBox()       { return fWorldBox; }
        std::string GetDetectorName()   { return fDetectorName; }
        BoundingBox GetDetectorBox()    { return fDetectorBox; }
        std::string GetCryostatName()   { return fCryostatName; }
        BoundingBox GetCryostatBox()    { return fCryostatBox; }
        int GetNumberOfTPCs()           { return fNumberOfTPCs; }
        std::vector<std::string> GetTPCNames()  { return fTPCNames; }
        std::vector<double> GetTPCMasses()      { return fTPCMasses; }
        std::vector<double> GetTPCDriftDistances()  { return fTPCDriftDistances; }
        BoundingBox GetTotalTPCBox()            { return fTotalTPCBox; }
        BoundingBox GetTotalActiveTPCBox()      { return fTotalActiveTPCBox; }
        double GetTotalTPCMass()                { return fTotalTPCMass; }

        // get quantities by index
        std::string GetTPCName(const size_t i);
        BoundingBox GetTPCBox(const size_t i);
        BoundingBox GetActiveTPCBox(const size_t i);
        double GetTPCMass(const size_t i);
        double GetTPCDriftDistance(const size_t i);

        // get volume information for a point
        DetectorVolume getVolume(std::vector<double> position);
        DetectorVolume getVolume(double x, double y, double z);

        // function for finding total tpc volumes
        void findTotalTPCBoxes();

        // fill the geometry ttree
        void FillTTree();
        
    private:
        ////////////////////////////////////////////////
        // Information which is automatically stored
        // meta variables
        art::ServiceHandle<geo::Geometry> fGeometryService;
        geo::GeometryCore const* fGeometryCore;

        // ROOT 
        art::ServiceHandle<art::TFileService> fTFileService;
        TTree *fGeometryTree;
        size_t fTriggerOffset;

        // map from volume names to volume type
        std::map<std::string,VolumeType> fVolumeTypeMap;
        // world volume
        std::string fWorldName;
        BoundingBox fWorldBox;
        // detector volume
        std::string fDetectorName;
        BoundingBox fDetectorBox;
        // cryostat volume
        std::string fCryostatName;
        BoundingBox fCryostatBox;
        // tpc volumes
        int fNumberOfTPCs;
        std::vector<std::string> fTPCNames;
        std::vector<BoundingBox> fTPCBoxes;
        std::vector<BoundingBox> fActiveTPCBoxes;
        std::vector<double> fTPCMasses;
        std::vector<double> fTPCDriftDistances;

        // full tpc volume
        BoundingBox fTotalTPCBox;
        BoundingBox fTotalActiveTPCBox;
        double fTotalTPCMass;

        ////////////////////////////////////////////////
        // detector material variables
        ////////////////////////////////////////////////
        // we will need to ask Geant4 about material 
        // properties for the detector volume
        // at each point of interest.  This requires holding 
        // this information in a
        // TGeoMaterial object, which is part of ROOT.
        const TGeoMaterial *fMaterial;
        geo::Point_t fMaterialPOI;
    };
}