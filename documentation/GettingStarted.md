#### 1. Set up the environment
##### 1.1 Use an Analysis Musing 

The easiest way to use the display is through an analysis musing environment. For example:

```
mu2einit
muse setup AnalysisMDC2025
```


##### 1.2 Build the GDML (generally not needed if working from a Musing)

If you are using your own local **muse** setup and want to incorporate the display, make sure to build the GDML geometry after cloning and building the Event Display repository. 

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

#### 2. Set up the $\text{.rootrc}$ file

The **$\text{.rootrc}$** file configures the ROOT environment and is required for proper startup and graphics behavior in the display.

    ```
    EventDisplay/config/makerootrc.sh
    ```
---
#### 3. Launch the event display 

Use one of the example FCL files in **EventDisplay/examples/** to visualize events from your art file:
```
    mu2e -c EventDisplay/examples/<example.fcl> -s <art file> 
```
At this stage, you may encounter issues such as missing reco data products or incorrect instance names. Refer to CommonErrors.md for common fixes.

* **Launch the remote browser**
  Set up an **SSH tunnel** from your local machine:
    ```
    ssh -KXY -L 0<port>:localhost:0<port> <username>@mu2egpvm0<machine>.fnal.gov
    ```
    Then open the generated URL (`http://localhost:<port>/win1/`) in your browser.
  > **Note:** **Google Chrome** is recommended.
  
* **Verify the display started up successfully**
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

#### 4. Troubleshooting stale browser sessions

If the Event Display hangs or fails to start correctly, the issue is often caused by stale browser sessions or leftover processes from a previous run.
Run the following command to terminate lingering browser and display processes:
    ```
    EventDisplay/config/kill.sh
    ```
