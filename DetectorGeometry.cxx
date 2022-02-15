#include "DetectorGeometry.h"

namespace extractor 
{
    DetectorGeometry* DetectorGeometry::sInstance{nullptr};
    std::mutex DetectorGeometry::sMutex;

    DetectorGeometry *DetectorGeometry::getInstance(const std::string& name)
    {
        std::lock_guard<std::mutex> lock(sMutex);
        if (sInstance == nullptr)
        {
            sInstance = new DetectorGeometry(name);
        }
        return sInstance;
    }

    std::string DetectorGeometry::GetTPCName(const size_t i) 
    {
        if (i < fTPCNames.size()) { return fTPCNames[i]; }
        else { return fTPCNames[0]; }
    }
    BoundingBox DetectorGeometry::GetTPCBox(const size_t i) 
    {
        if (i < fTPCBoxes.size()) { return fTPCBoxes[i]; }
        else { return fTPCBoxes[0]; }
    }
    BoundingBox DetectorGeometry::GetActiveTPCBox(const size_t i) 
    {
        if (i < fActiveTPCBoxes.size()) { return fActiveTPCBoxes[i]; }
        else { return fActiveTPCBoxes[0]; }
    }
    std::vector<double> DetectorGeometry::GetTPCMasses() 
    { 
        return fTPCMasses; 
    }
    double DetectorGeometry::GetTPCMass(const size_t i) 
    {
        if (i < fTPCMasses.size()) { return fTPCMasses[i]; }
        else { return fTPCMasses[0]; }
    }
    std::vector<double> DetectorGeometry::GetTPCDriftDistances() 
    { 
        return fTPCDriftDistances; 
    }
    double DetectorGeometry::GetTPCDriftDistance(const size_t i) 
    {
        if (i < fTPCDriftDistances.size()) { return fTPCDriftDistances[i]; }
        else { return fTPCDriftDistances[0]; }
    }
    BoundingBox DetectorGeometry::GetTotalTPCBox()        
    { 
        return fTotalTPCBox; 
    }
    BoundingBox DetectorGeometry::GetTotalActiveTPCBox()  
    { 
        return fTotalActiveTPCBox; 
    }
    double DetectorGeometry::GetTotalTPCMass()            
    { 
        return fTotalTPCMass; 
    }
    ///////////////////////////////////////////////////////////////////////////////////////
    DetectorGeometry::DetectorGeometry(const std::string name)
    : sName(name)
    {
        // set up the geometry interface
        fGeometryCore = lar::providerFrom<geo::Geometry>();
        // initialize TTrees
        fGeometryTree = fTFileService->make<TTree>("Geometry", "Geometry");
        // get detector clock data
        auto const clock_data = 
            art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
        fTriggerOffset = trigger_offset(clock_data);

        // collect world info
        fWorldName = fGeometryCore->GetWorldVolumeName();
        fWorldBox.setBox(fGeometryCore->WorldBox());
        // create name-volumetype map for world
        fMaterialPOI.SetCoordinates(fWorldBox.x_min,fWorldBox.y_min,fWorldBox.z_min);
        std::string volumeName = fGeometryCore->VolumeName(fMaterialPOI);
        fVolumeTypeMap[volumeName] = VolumeType::World;
        // collect detector info
        fDetectorName = fGeometryCore->DetectorName();
        fDetectorBox.setBox(-fGeometryCore->DetHalfWidth(), fGeometryCore->DetHalfWidth(),
                            -fGeometryCore->DetHalfHeight(), fGeometryCore->DetHalfHeight(),
                            0, fGeometryCore->DetLength());
        // collect cryostat info
        // for now, assuming analysis is done over a single cryostat
        geo::CryostatGeo const& Cryo = fGeometryCore->Cryostat();
        fCryostatName = std::string(Cryo.ID());
        fCryostatBox.setBox(Cryo.Boundaries());
        // create name-volumetype map for cryostat
        fMaterialPOI.SetCoordinates(fCryostatBox.x_min,fCryostatBox.y_min,fCryostatBox.z_min);
        volumeName = fGeometryCore->VolumeName(fMaterialPOI);
        fVolumeTypeMap[volumeName] = VolumeType::Cryostat;
        // iterate over all TPCs
        fNumberOfTPCs  = fGeometryCore->TotalNTPC();
        for (geo::TPCGeo const& TPC : fGeometryCore->IterateTPCs())
        {
            fTPCNames.emplace_back(TPC.ID());
            fTPCBoxes.emplace_back(BoundingBox(TPC.BoundingBox()));
            fActiveTPCBoxes.emplace_back(BoundingBox(TPC.ActiveBoundingBox()));
            fTPCMasses.emplace_back(TPC.ActiveMass());
            fTPCDriftDistances.emplace_back(TPC.DriftDistance());
            // create name-volumetype map for this tpc
            fVolumeTypeMap[fGeometryCore->VolumeName(TPC.GetCenter())] = VolumeType::TPC;
        }
        // find the total TPC and total Active TPC volumes
        findTotalTPCBoxes();
        fTotalTPCMass = fGeometryCore->TotalMass();    
    }
    // get volume information for a point
    DetectorVolume DetectorGeometry::getVolume(std::vector<double> position)
    {

        fMaterialPOI.SetCoordinates(position[0],position[1],position[2]);
        // get the volume information
        //std::cout << "volname: " << fGeometryCore->VolumeName(fMaterialPOI) << std::endl;
        std::string volumeName = fGeometryCore->VolumeName(fMaterialPOI);
        VolumeType volumeType = fVolumeTypeMap[volumeName];
        // get the current material information
        fMaterial = fGeometryService->Material(fMaterialPOI);
        double material = fMaterial->GetZ();
        std::string materialName = fMaterial->GetName();
        //std::cout << "mat name: " << fMaterial->GetName() << std::endl;
        // return the constructed volume 
        return DetectorVolume(volumeType, volumeName, materialName, material);
    }
    // get volume information for a point
    DetectorVolume DetectorGeometry::getVolume(double x, double y, double z)
    {

        fMaterialPOI.SetCoordinates(x,y,z);
        // get the volume information
        //std::cout << "volname: " << fGeometryCore->VolumeName(fMaterialPOI) << std::endl;
        std::string volumeName = fGeometryCore->VolumeName(fMaterialPOI);
        VolumeType volumeType = fVolumeTypeMap[volumeName];
        // get the current material information
        fMaterial = fGeometryService->Material(fMaterialPOI);
        double material = fMaterial->GetZ();
        std::string materialName = fMaterial->GetName();
        //std::cout << "mat name: " << fMaterial->GetName() << std::endl;
        // return the constructed volume 
        return DetectorVolume(volumeType, volumeName, materialName, material);
    }
    // get total tpc volume information
    void DetectorGeometry::findTotalTPCBoxes()
    {
        double x_min = 0; double x_max = 0;
        double y_min = 0; double y_max = 0;
        double z_min = 0; double z_max = 0;
        for (size_t i = 0; i < fTPCBoxes.size(); i++) {
            if (fTPCBoxes[i].x_min < x_min) x_min = fTPCBoxes[i].x_min;
            if (fTPCBoxes[i].x_max > x_max) x_max = fTPCBoxes[i].x_max;
            if (fTPCBoxes[i].y_min < y_min) y_min = fTPCBoxes[i].y_min;
            if (fTPCBoxes[i].y_max > y_max) y_max = fTPCBoxes[i].y_max;
            if (fTPCBoxes[i].z_min < z_min) z_min = fTPCBoxes[i].z_min;
            if (fTPCBoxes[i].z_max > z_max) z_max = fTPCBoxes[i].z_max;
        }
        fTotalTPCBox.setBox(x_min, x_max, y_min, y_max, z_min, z_max);
        x_min = 0; x_max = 0;
        y_min = 0; y_max = 0;
        z_min = 0; z_max = 0;
        for (size_t i = 0; i < fActiveTPCBoxes.size(); i++) {
            if (fActiveTPCBoxes[i].x_min < x_min) x_min = fActiveTPCBoxes[i].x_min;
            if (fActiveTPCBoxes[i].x_max > x_max) x_max = fActiveTPCBoxes[i].x_max;
            if (fActiveTPCBoxes[i].y_min < y_min) y_min = fActiveTPCBoxes[i].y_min;
            if (fActiveTPCBoxes[i].y_max > y_max) y_max = fActiveTPCBoxes[i].y_max;
            if (fActiveTPCBoxes[i].z_min < z_min) z_min = fActiveTPCBoxes[i].z_min;
            if (fActiveTPCBoxes[i].z_max > z_max) z_max = fActiveTPCBoxes[i].z_max;
        }
        fTotalActiveTPCBox.setBox(x_min, x_max, y_min, y_max, z_min, z_max);
    }
    void DetectorGeometry::FillTTree()
    {
        // add geometry info
        fGeometryTree->Branch("world_name", &fWorldName);
        fGeometryTree->Branch("world_box_ranges", &(fWorldBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
        fGeometryTree->Branch("detector_name", &fDetectorName);
        fGeometryTree->Branch("detector_box_ranges", &(fDetectorBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
        fGeometryTree->Branch("cryostat_name", &fCryostatName);
        fGeometryTree->Branch("cryostat_box_ranges", &(fCryostatBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
        fGeometryTree->Branch("number_of_tpcs", &fNumberOfTPCs);
        fGeometryTree->Branch("tpc_names", &fTPCNames);
        for (int i = 0; i < fNumberOfTPCs; i++) {
            fGeometryTree->Branch(std::string("tpc_"+std::to_string(i)+"_name").c_str(), &(fTPCNames[i]));
            fGeometryTree->Branch(std::string("tpc_"+std::to_string(i)+"_box_ranges").c_str(), &(fTPCBoxes[i]), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
            fGeometryTree->Branch(std::string("tpc_"+std::to_string(i)+"_mass").c_str(), &(fTPCMasses[i]));
            fGeometryTree->Branch(std::string("tpc_"+std::to_string(i)+"_drift_distance").c_str(), &(fTPCDriftDistances[i]));
        }
        fGeometryTree->Branch("tpc_masses", &fTPCMasses);
        fGeometryTree->Branch("tpc_drift_distances", &fTPCDriftDistances);
        fGeometryTree->Branch("total_tpc_box_ranges", &(fTotalTPCBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
        fGeometryTree->Branch("total_active_tpc_box_ranges", &(fTotalActiveTPCBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
        fGeometryTree->Branch("total_tpc_mass", &fTotalTPCMass);
        fGeometryTree->Fill();
    }
}