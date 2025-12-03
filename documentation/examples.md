# ✨ Examples and Configuration

This section provides details on using example files and configuring the Event Display for different geometry settings.

---

## 1. The `prolog.fcl` File

The **`prolog.fcl`** file defines the **default settings** for the Event Display, covering both the geometry and event scene information (refer to the FclOption documentation for specifics).

* **Override Nominal Settings:** To use a non-default geometry or configuration, you must **update your driver `.fcl` file**. This is done in the same way you would configure any Mu2e/Offline analyzer module.

---

## 2. The `EventDisplay/examples1` Directory

This directory contains several starter `.fcl` files to help you quickly set up the Event Display for various sample types and geometry configurations.

| Example File | Description | Usage |
| :--- | :--- | :--- |
| `nominal_MDC2020.fcl` | Uses the original CRV geometry. | Recommended for **MDC2020** nominal samples. |
| `nominal_MDC2025.fcl` | Uses the updated MDC2025 CRV geometry. | Recommended for **MDC2025** nominal samples, assuming the common GDML. |
| `extracted_MDC2020.fcl` | Used for running against MDC2020 extracted samples. | **Note:** This example will be retired soon. |
| `extracted_MDC2025.fcl` | Used for extracted physics studies with MDC2025 samples. | Recommended for **MDC2025** extracted samples. |

> **Tip:** These `.fcl` files are intended as starting points and can be **edited** to add or remove specific features as needed for your study.

---

## 3. Other Geometry Configurations

If you need to run the Event Display against a geometry configuration that is **not** one of the standard examples (e.g., a custom geometry tag), you must explicitly update the paths in your driver `.fcl` file.

To configure the display for a new geometry, update the following parameters:

```
physics.analyzers.Mu2eEventDisplay.gdmlname :"Offline/gen/gdml/<new>.gdml"
services.GeometryService.inputFile: "Offline/Mu2eG4/geom/<new-geom_common.txt>"
```

* Replace `<new>.gdml` with the name of the GDML file for your configuration.

* Replace `<new-geom_common.txt>` with the name of the corresponding geometry input file.

