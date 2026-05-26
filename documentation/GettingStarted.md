Before running the Event Display, complete steps 1 and 2 to prevent common errors and system hangs.

#### 1a. Use a Musing 

The easiest way to use the EventDisplay is by using our Analysis Musing. For example:

```
mu2einit
muse setup AnalysisMDC2025
```


#### 1b. Building GDML Geometry (generally not needed if working from a Musing)

The Event Display requires the detector geometry to be compiled into a **GDML** (Geometry Description Markup Language) file. Skipping this step will result in a fatal ROOT error, as the application cannot initialize its environment without the geometry data. 

* **Command to Build GDML:**
    ```
    muse build GDML
    ```

* **Common Error if Skipped:**
    ```
     what(): ---- FatalRootError BEGIN
     Fatal Root Error: TFile::Write
     file /dev/null not opened in write mode
     ROOT severity: 2000
    ---- FatalRootError END
    ```
    This error indicates a critical failure to load the necessary geometry components.

---

#### 2. Setting up the $\text{.rootrc}$ File

The **$\text{.rootrc}$ configuration file** is crucial for the correct functioning of the **ROOT environment** and the Event Display, especially for graphical and startup behavior.

* **Command to Create $\text{.rootrc}$:**
    This script generates a properly configured $\text{.rootrc}$ file in your current working directory.
    ```
    EventDisplay/config/makerootrc.sh
    ```

---
#### 3. Viewing events 

Use one of the example files provided in EventDisplay/examples/ to view the events in your art file like you would do with any other mu2e process

```
    mu2e -c EventDisplay/examples/<example>.fcl -c <art file> 
```

* **Launch the remote browser**
   set up an **SSH tunnel** in your local machine or on the remote machine:
    ```
    ssh -KXY -L 0<port>:localhost:0<port> <user>@mu2egpvm0<machine>.fnal.gov
    ```
    Then, copy the URL (e.g., `http://localhost:<port>/win1/`) into your chosen browser.
  > **Note:** **Google Chrome** is the recommended browser.
  
* **Sign that the Display is Accessible:**
    The local (mu2e) machine will print a log that includes the **EVE URL** and streaming details, indicating that the geometry and event scenes have loaded successfully.

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
    The presence of the `EVE URL` and logs for streaming the `[Geometry scene]` and `[Event scene]` confirms it's safe to attempt remote access.

---

#### 5. Correcting Stale Browser Sessions (Correcting hangs)

If the event display **hangs** or fails to start correctly, it is often caused by leftover processes or stale browser calls from a previous session.

* **Troubleshooting Command:**
    Executing this script attempts to forcefully terminate any persistent browser connections that could interfere with the new display session.
    ```
    EventDisplay/config/kill.sh
    ```


