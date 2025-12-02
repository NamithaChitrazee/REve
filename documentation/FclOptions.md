# Mu2eEventDisplay Configuration Documentation

There are many parameters defined within the Mu2e/EventDisplay Framework. Here we list and define all of them.

### 💡 Note on Parameter Placement (Module vs. Filler)

 In this type of framework:

1. **`filler` Block:** Is responsible for **data ingestion**. It specifies *which* collections (e.g., `ComboHitCollection`, `KalSeedPtrCollection`) should be read from the event record and *where* to find them (the module tags, like `makeSH` or `MergeKK`).

2. **Module-Level Flags:** Control the **rendering/visualization logic** of the main module. Flags like `addTrkStrawHits` or `addErrBar` tell the `Mu2eEventDisplay` module *how* to draw the data that the `filler` has already prepared. They affect display features (error bars, drawing style, etc.) rather than the data source.


## General Module Parameters

These parameters control the core behavior of the module, including logging, geometry loading, and global display switches for major detector components.

| Parameter | Type | Default (from config) | Description | 
| :--- | :--- | :--- | :--- | 
| `module_type` | String | `Mu2eEventDisplay` | **Mandatory.** Specifies the class name of the module to load. | 
| `diagLevel` | Integer | `1` | Controls the verbosity of the internal module logging. `0` is typically minimal, higher values (`1`, `2`, etc.) increase diagnostic output. | 
| `gdmlname` | String | `"Offline/gen/gdml/mu2e.gdml"` | Specifies the path to the GDML (Geometry Description Markup Language) file used to load the detector geometry for visualization. | 
| `extracted` | Boolean | `false` | flag to control if the display should use an *extracted* geometry | 
| `particles` | Array (Int) | `[11,13,2212,2112,211,22,212]` | A list of particle IDs (PDG codes) to focus the display on or track. Common IDs include electron (`11`), muon (`13`), proton (`2212`), neutron (`2112`), pion (`211`), photon (`22`), and others. | 
| `specifyTag` | Boolean | `false` | **IMPORTANT.** When set to `true`, it overrides default behavior and forces the module to use the collection tags specified explicitly within the `filler` block. | 
| `seqMode` | Boolean | `true` | **IMPORTANT.** This will allow you to navigate from the first to the last using the NextEvent() button | 


## Detector Visibility Flags

These boolean flags enable or disable the display of the major sub-detector components (geometry and related elements).

| Parameter | Type | Default (from config) | Description | 
| :--- | :--- | :--- | :--- | 
| `showCrv` | Boolean | `true` | Show the Cosmic Ray Veto (CRV) detector components. | 
| `showPS` | Boolean | `false` | Show the Production Solenoid (PS) components. | 
| `showTS` | Boolean | `false` | Show the Transport Solenoid (TS) components. | 
| `showDS` | Boolean | `false` | Show the Detector Solenoid (DS) components. | 
| `showST` | Boolean | `true` | Show the Stopping Target (ST). | 
| `showSTM` | Boolean | `false` | Show the Stopping Target Monitor (STM) components. | 
| `showCalo` | Boolean | `true` | Show the Calorimeter volume. | 
| `showCaloCrystals` | Boolean | `true` | Show the individual calorimeter crystals. | 
| `showTracker` | Boolean | `true` | Show the Tracker detector components (straw tubes). | 

## Display Functionality and Rendering Flags

These flags modify how data is rendered or what auxiliary information is included.

| Parameter | Type | Default (from config) | Description | 
| :--- | :--- | :--- | :--- | 
| `show2D` | Boolean | `true` | Enables or disables a two-dimensional projection display mode. | 
| `caloVST` | Boolean | `false` | A visualization mode intended to view the Calorimeter alone. | 
| `addErrBar` | Boolean | `true` | Add error bars to visualized hit/cluster points. | 
| `addCrystalHits` | Boolean | `false` | Display individual hits within the calorimeter crystals. | 
| `addCrvBars` | Boolean | `false` | Display the physical CRV bars/scintillators. | 
| `addKalInter` | Boolean | `false` | Add Kalmar (KalSeed) trajectory intersection points. | 
| `addTrkStrawHits` | Boolean | `false` | Add individual tracker straw hits to the display. | 
| `addTrkCaloHits` | Boolean | `false` | Add tracker hits extrapolated to the calorimeter face. | 
| `strawdisplay` | Boolean | `false` | Specific flag to control the display of individual straw tubes (often only enabled for debug/close-up views). | 


## `filler` Block Parameters

The `filler` block is responsible for configuring the **data collections** that the event display will ingest and draw. The keys are the data types, and the values are the **producer tags** specifying which module created the collection.

### Input Data Collections

| Parameter | Type | Default (from config) | Description | 
| :--- | :--- | :--- | :--- | 
| `ComboHitCollection` | Array (String) | `["makeSH"]` | Tag(s) for the combined straw hits used for tracking. | 
| `CrvRecoPulseCollection` | Array (String) | `["SelectRecoMC"]` | Tag(s) for the reconstructed pulses from the CRV. | 
| `CrvCoincidenceClusterCollection`| Array (String) | `["SelectRecoMC:CrvCoincidenceClusterFinder"]` | Tag(s) for the CRV coincidence clusters (higher-level CRV data). | 
| `TimeClusterCollection` | Array (String) | `["MHDeM"]` | Tag(s) for the primary time clusters found in the event. | 
| `CaloDigiCollection` | Array (String) | `["CaloDigiMaker"]` | Tag(s) for the digitized/raw calorimeter data (digis). | 
| `CaloClusterCollection` | Array (String) | `["CaloClusterMaker"]` | Tag(s) for the reconstructed calorimeter clusters. | 
| `KalSeedPtrCollection` | Array (String) | `["MergeKK"]` | Tag(s) for the final, merged Kalman Filter (KalSeed) tracks. | 
| `HelixSeedCollection` | Array (String) | `["MHFinderDe"]` | Tag(s) for the initial helix seeds used for tracking. | 
| `CosmicTrackSeedCollection`| String | `"CosmicTrackFinderTimeFit"` | Tag for cosmic ray tracks (if available). | 
| `MCTrajectoryCollection` | Array (String) | `["compressRecoMCs"]` | Tag(s) for the Monte Carlo (MC) generated trajectories (the truth data). | 
| `SurfaceStepCollection` | Array (String) | `[]` | Tag(s) for MC steps on detector surfaces (currently empty). | 
| `SimParticleCollection` | Array (String) | `["compressSTMDet"]` | Tag(s) for the MC simulated particles. | 

### Filler Display Control Flags

These flags control whether the data collections identified above are actually *added* to the event display's drawing list.

| Parameter | Type | Default (from config) | Description | 
| :--- | :--- | :--- | :--- | 
| `diagLevel` | Integer | `0` | Logging level specific to the `filler` block (data reading). | 
| `FillAll` | Boolean | `false` | If true, likely overrides individual `add...` flags and attempts to add all known data types. | 
| `addHits` | Boolean | `false` | Add `ComboHitCollection` (combined tracker hits). | 
| `addCrvRecoPulse` | Boolean | `false` | Add `CrvRecoPulseCollection`. | 
| `addCrvClusters` | Boolean | `false` | Add `CrvCoincidenceClusterCollection`. | 
| `addTimeClusters` | Boolean | `false` | Add `TimeClusterCollection`. | 
| `addTrkHits` | Boolean | `false` | Add tracker hits (often synonymous with `addHits` or a specialized hit collection). | 
| `addCaloDigis` | Boolean | `false` | Add `CaloDigiCollection`. | 
| `addClusters` | Boolean | `false` | Add `CaloClusterCollection`. | 
| `addKalSeeds` | Boolean | `false` | Add `KalSeedPtrCollection` (reconstructed tracks). | 
| `addCosmicTrackSeeds` | Boolean | `false` | Add `CosmicTrackSeedCollection` (cosmic tracks). | 
| `addMCTraj` | Boolean | `true` | Add `MCTrajectoryCollection` (MC truth particle trajectories). | 
| `addSurfSteps` | Boolean | `true` | Add `SurfaceStepCollection` (MC steps on surfaces). | 
| `addSimParts` | Boolean | `true` | Add `SimParticleCollection` (MC truth particles). |
