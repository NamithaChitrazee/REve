#### 1. Set up the environment
##### 1.1 Use a Musing 

The easiest way to use the EventDisplay is by using our Analysis Musing. For example:

```
mu2einit
muse setup AnalysisMDC2025
```


##### 1.2 Build GDML Geometry (generally not needed if working from a Musing)

But if you have your own local muse setup and want to incorporate the Event Display as well, do not forget the following step after cloning and building the Event Display repository. 

* **Command to Build GDML:**
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

The **$\text{.rootrc}$** file configures the ROOT environment and is required for proper Event Display startup and graphics behavior.

* **Command to Create $\text{.rootrc}$:**
    This script generates a properly configured $\text{.rootrc}$ file in your current working directory.
    ```
    EventDisplay/config/makerootrc.sh
    ```
---
#### 3. Launch the event display 

Use one of the example files provided in EventDisplay/examples/ to view the events in your art file like you would do with any other mu2e process

```
    mu2e -c EventDisplay/examples/<example.fcl> -s <art file> 
```

* **Launch the remote browser**
  Set up an **SSH tunnel** from your local machine:
    ```
    ssh -KXY -L 0<port>:localhost:0<port> <username>@mu2egpvm0<machine>.fnal.gov
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

#### 4. Troubleshooting stale browser sessions

If the event display **hangs** or fails to start correctly, it is often caused by leftover processes or stale browser calls from a previous session.

* **Troubleshooting Command:**
    Executing this script attempts to forcefully terminate any persistent browser connections that could interfere with the new display session.
    ```
    EventDisplay/config/kill.sh
    ```
