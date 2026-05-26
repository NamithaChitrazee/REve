# 1. Set up the environment
## 1.1 Use an Analysis Musing 

The easiest way to use the display is through an analysis musing environment. For example:

```
mu2einit
muse setup AnalysisMDC2025
```

## 1.2 Build GDML (generally not needed if working from a Musing)

If you are using your own local **muse** setup and want to incorporate the display, make sure to build the GDML after cloning and building the Event Display repository. 

    ```
    muse build GDML
    ```
Skipping this step can cause a fatal ROOT error as shown below.

    ```
     what(): ---- FatalRootError BEGIN
     Fatal Root Error: TFile::Write
     file /dev/null not opened in write mode
     ROOT severity: 2000
    ---- FatalRootError END
    ```

---

# 2. Set up the $\text{.rootrc}$ file

The **$\text{.rootrc}$** file configures the ROOT environment and is required for proper start-up and graphics behavior in the display.

    ```
    EventDisplay/config/makerootrc.sh
    ```
---
# 3. Launch the event display 
Use one of the example FCL files in **EventDisplay/examples/** to visualize events in your art file:
```
    mu2e -c EventDisplay/examples/<example.fcl> -s <art file> 
```
## Launch the remote browser
If running directly on a remote mu2e machine, execute:
```
source EventDisplay/config/start_RemoteDisplay.sh --port <WXYZ> --user <username> --machine <mu2e gpvm number e.g. 04>
```

If you are on on your local machine, establish an SSH tunnel:
```
ssh -KXY -L 0<port>:localhost:0<port> <username>@mu2egpvm0<machine>.fnal.gov
```
Then open the generated URL (`http://localhost:<port>/win1/`) in your browser.
> **Note:** **Google Chrome** is recommended.

## Verify the display started up successfully
The remote mu2e machine should print logs containing the EVE URL and scene streaming information, indicating that the geometry and event scenes loaded successfully.

```
Info in <THttpEngine::Create>: Starting HTTP server on port 1234
EVE URL http://localhost:1234/win1/
...
EVEMNG ............. streaming the world scene.
...
EVEMNG ............. streaming scene [Geometry scene]
...
EVEMNG ............. streaming scene [Event scene]
...
Geometry file: /exp/mu2e/app/users/sophie/...
...
```

---

## 4. Troubleshooting stale browser sessions

If the Event Display hangs or fails to start correctly, the issue is often caused by stale browser sessions or leftover processes from a previous run.
Run the following command to terminate lingering browser and display processes:
```
EventDisplay/config/kill.sh
```
