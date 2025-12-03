## 🚀 Getting Started with the Event Display

This guide covers the essential setup and troubleshooting steps required before launching the Event Display to ensure a smooth and error-free operation.

---

### 🛠️ Initial Setup and Troubleshooting

Before running the Event Display, complete the following steps to prevent common errors and system hangs.

#### 1. Building GDML Geometry

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

#### 3. Launching the Remote Browser

Once you've run the display on the local **Mu2e machines** (refer to the documentation on EventNavigation for local instructions), you can access the remote display using a browser.

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

* **Accessing the Remote Display:**
    Use the provided script on your remote machine to set up the connection.

    ```
    start_RemoteDisplay.sh --port <WXYZ> --user <your_name> --machine <mu2e gpvm number e.g. 04>
    ```

    > **Note:** **Google Chrome** is the recommended browser.

* **Manual Setup for Other Browsers:**
    If you prefer another browser, you can manually set up an **SSH tunnel**:
    ```
    ssh -KXY -L 0<port>:localhost:0<port> <user>@mu2egpvm0<machine>.fnal.gov
    ```
    Then, copy the displayed URL (e.g., `http://localhost:<port>/win1/`) into your chosen browser.

---

#### 4. Correcting Stale Browser Sessions (Correcting hangs)

If the event display **hangs** or fails to start correctly, it is often caused by leftover processes or stale browser calls from a previous session.

* **Troubleshooting Command:**
    Executing this script attempts to forcefully terminate any persistent browser connections that could interfere with the new display session.
    ```
    EventDisplay/config/kill.sh
    ```


