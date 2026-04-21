# CRV Z-Shift Analysis for Extracted Geometry

## Summary
A Z-shift of **1635.5 mm** (divided by 10 for cm display units) is needed to align CRV position data with CRV geometry in extracted geometry. This shift is applied in both `MainWindow::GeomDrawerExtracted()` and `DataInterface` for CRV data visualization.

## Key Findings

### 1. DetectorSystem::toDetector() Implementation
**File**: [Offline/GeometryService/inc/DetectorSystem.hh](Offline/GeometryService/inc/DetectorSystem.hh#L28-L32)

```cpp
// Takes a 3-vector in the mu2e system; return it in the detector system.
CLHEP::Hep3Vector toDetector( CLHEP::Hep3Vector const& v ) const {
    return v - _origin;  // Subtracts the detector system origin
}
```

**What it does**: Converts from Mu2e frame coordinates to detector frame by subtracting a single offset point (the detector origin).

The origin is set in [Offline/GeometryService/src/DetectorSystemMaker.cc](Offline/GeometryService/src/DetectorSystemMaker.cc#L15-L22):
```cpp
CLHEP::Hep3Vector origin( - config.getDouble("mu2e.solenoidOffset"),
                          0.,
                          config.getDouble("mu2e.detectorSystemZ0")
                        );
```

### 2. Coordinate System Origins

#### Nominal Geometry Configuration
- **mu2e.detectorSystemZ0** = 10171 mm (source: multiple geom files like `geom_2021_PhaseI.txt`)
- **tracker.z0** = 10171 mm
- Both tracker center and detector origin are at **the same Z position**
- **tracker.mother.halfLength** = 1635.1 mm (or 1635.5 mm in newer configs like `tracker_v7.txt`)

#### Extracted Geometry Configuration
- [Offline/Mu2eG4/geom/geom_common_extracted_v02.txt](Offline/Mu2eG4/geom/geom_common_extracted_v02.txt#L8-L60)
  - **inGaragePosition** = true
  - **garage.zOffset** = 14000.0 mm (Z-displacement of detector to garage/test stand)
  - **mu2e.detectorSystemZ0** = 24171 mm
  - **tracker.z0** = 24171 mm
  - **tracker.mother.z0** = 24175 mm
  - Difference from nominal: 24171 - 10171 = **14000 mm** (matches garage.zOffset exactly)

**Key insight**: In both geometries, tracker.z0 equals mu2e.detectorSystemZ0, meaning the detector coordinate system origin is centered on the tracker Z position.

### 3. TrackerMother Half-Length
**File**: [Offline/GeometryService/src/TrackerMaker.cc](Offline/GeometryService/src/TrackerMaker.cc#L106-L120)

Configuration parameter: `tracker.mother.halfLength`

Values found across geometry configs:
- tracker_v3.txt through v7.txt, cd3*.txt: **1635.1 mm** or **1635.5 mm** (most common)
- geom_run1_b_v01.txt: 1640 mm (increased to contain end rings)

**Tracker extents** (for extracted geometry with halfLength=1635.5 mm):
- Z min (upstream): 24171 - 1635.5 = **22535.5 mm**
- Z center: **24171 mm**
- Z max (downstream): 24171 + 1635.5 = **25806.5 mm**

### 4. CRV Positioning in Extracted Geometry

**File**: [Offline/Mu2eG4/geom/crv_counters_extracted_v02.txt](Offline/Mu2eG4/geom/crv_counters_extracted_v02.txt#L37-L39)

CRV first counter Z positions (in Mu2e frame, in mm):
- **CRV-EX**: 21127 mm
- **CRV-T1**: 22756 mm
- **CRV-T2**: 21955 mm

These are positioned in the Mu2e coordinate system and must be transformed through:
1. `toDetector()` - converts Mu2e frame to detector frame
2. Display offset - applies the 1635.5/10.0 cm shift

### 5. The 1635.5 mm Shift Application

#### In MainWindow (Geometry Visualization)
**File**: [MainWindow.cc](src/MainWindow.cc#L664-L689)

Applied to all CRV components in extracted geometry:
```cpp
// For CRV_EX
shift.at(2) = z_crvex - z_trk + 1635.5/10.0;  // Z-shift in cm

// For CRV_T1
shift.at(2) = z_crvt1 - z_trk + 1635.5/10.0;

// For CRV_T2  
shift.at(2) = z_crvt2 - z_trk + 1635.5/10.0;
```

Where z_crvex, z_crvt1, z_crvt2 are extracted from GDML geometry hierarchy and z_trk is tracker position.

#### In DataInterface (Data Visualization)
**File**: [DataInterface.cc](src/DataInterface.cc#L710-L835)

Applied when displaying CRV reconstruction data for extracted geometry:
```cpp
// For CRV hits
CLHEP::Hep3Vector pointInMu2e = det->toDetector(crvCounterPos);
double hit_z = pointmmTocm(pointInMu2e.z());
if(extracted) hit_z += 1635.5/10.0;  // Same correction as MainWindow for CRV display alignment

// For CRV bars and clusters
Z_cm = pointmmTocm(pointInMu2e.z());
// ... later when setting vertices ...
b->SetVertex(..., Z_cm ± length + 1635.5/10.0);
```

### 6. Why Is This Shift Needed?

The 1635.5 mm value equals the **TrackerMother half-length**. The shift appears necessary due to:

1. **Coordinate system initialization difference**: 
   - In nominal geometry, CRV and tracker are positioned consistently within the detector frame
   - In extracted geometry, CRV is "extracted" to a garage/test stand position but maintains coordinates in the Mu2e reference frame
   - The hierarchical Z-positions from GDML differ between extracted and nominal geometries

2. **Display frame alignment**:
   - Display centers on tracker (shifting it to origin)
   - Tracker extends ±1635.5 mm from its center
   - CRV data points, when transformed through `toDetector()`, need this additional offset to align visually with CRV GDML geometry elements

3. **GDML Hierarchy vs. Coordinate Transformation Mismatch**:
   - MainWindow extracts component positions from GDML hierarchy traversal and cumulative offsets
   - DataInterface uses `toDetector()` transformation on raw CRV data coordinates
   - These two approaches produce Z-displacements that differ by ~1635.5 mm in extracted geometry

### 7. Geometry Visualization vs. Data Point Alignment

The shift ensures:
- **CRV bars** (drawn from GDML geometry via MainWindow offsets) align with
- **CRV hit points** (transformed via `toDetector()` in DataInterface)

Both must display at the same physical location for correct visualization correlation.

## Code Locations

| Purpose | File | Lines |
|---------|------|-------|
| toDetector() definition | Offline/GeometryService/inc/DetectorSystem.hh | 28-32 |
| toDetector() initialization | Offline/GeometryService/src/DetectorSystemMaker.cc | 15-22 |
| Geometry display shift (MainWindow) | src/MainWindow.cc | 664-689 |
| CRV hit data shift (DataInterface) | src/DataInterface.cc | 710-835 |
| Extracted geometry config | Offline/Mu2eG4/geom/geom_common_extracted_v02.txt | 1-70 |
| CRV geometry extracted config | Offline/Mu2eG4/geom/crv_counters_extracted_v02.txt | 1-50 |
| Tracker geometry config | Offline/Mu2eG4/geom/tracker_v7.txt | 1-100 |
| TrackerMaker implementation | Offline/GeometryService/src/TrackerMaker.cc | 100-150 |

## Conclusion

The ~1635.5 mm Z-shift is **geometry-specific to extracted detector configuration** and stems from:
1. **TrackerMother.halfLength = 1635.5 mm** - a constant geometric property
2. **Different coordinate reference frames** between GDML geometry hierarchy positions and `toDetector()` transformations in extracted configuration
3. **Need to align display geometry visualization with analysis data point coordinates**

The shift is **not arbitrary** but rather reflects the fundamental difference between how tracker dimensions (±halfLength from center) relate to CRV positioning when the detector is extracted to a different Z=24171 mm coordinate origin versus the nominal Z=10171 mm location.
