# Getting Started

## Initial Setup and Troubleshooting

Before running the Event Display, ensure the following setup steps are completed to prevent common errors and hangs.

### Building GDML Geometry

The Event Display requires the detector geometry to be correctly compiled into a GDML (Geometry Description Markup Language) file. Failure to do so results in a fatal ROOT error because the application cannot properly initialize its environment.

Command to Build GDML:

```muse build GDML```


Common Error if Skipped:

```
 what(): ---- FatalRootError BEGIN
 Fatal Root Error: TFile::Write
 file /dev/null not opened in write mode
 ROOT severity: 2000
---- FatalRootError END
```

This error indicates a critical failure to load the necessary geometry components.

### Setting up the .rootrc File

The $\text{.rootrc}$ configuration file is essential for ensuring the ROOT environment and the Event Display function correctly, especially regarding graphical configuration and startup behavior. Running this script creates a properly configured $\text{.rootrc}$ file in your current working directory.

Command to Create .rootrc:

```EventDisplay/config/makerootrc.sh```


### Correcting Stale Browser Sessions (Correcting hangs)

If the event display hangs or fails to start correctly upon execution, it may be due to leftover processes or stale browser calls from previous sessions.

Troubleshooting Command:

```EventDisplay/config/kill.sh```


Executing this script attempts to forcefully terminate any persistent browser connections that could interfere with the new display session, correcting unexpected hangs.


