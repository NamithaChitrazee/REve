## 🗺️ Event Navigation

During analysis, there are several methods available to navigate to a specific event within the Event Display.

---

### 1. Default Sequential Mode

The **default navigation mode** is sequential.

* By default, the setting `seqMode` is set to `true`.
* The display begins at the **first event**.
* Users can navigate to the next event in the sequence using the **GUI's navigation controls**. 

---

### 2. Autoplay

The **Autoplay** option overrides the standard GUI navigation.

* When the "**autoplay**" checkbox is selected, the display will **stream all events sequentially** until it reaches the end of the file.

---

### 3. Go to an Event: GUI

The Event Display GUI includes **text entry boxes** that allow a user to jump directly to an event by specifying its **run** and **event number**.

* **Steps to Use (Currently a Placeholder):**
    1.  Ensure the display is loaded in `seqMode`.
    2.  Navigate to the **second event** using the `NextEvent()` function.
    3.  Enter the desired run and event number into the text boxes.
    4.  Select the **`Go()`** button, followed immediately by **`NextEvent()`**.

> **Note:** This is **not the recommended method** for navigation and is considered a placeholder feature.

---

### 4. Go to an Event: Scripts (Recommended Method)

The most robust way to navigate to a specific event is by using the custom script located in this repository: `config/EventDisplay.sh`.

#### A. Navigating with the Dataset Name

If you know the name of the **dataset** containing the event, use the following syntax:

```bash
./EventDisplay.sh --run 1201 --subrun 34 --event 15028 --dataset mcs.mu2e.ensembleMDS2cMix1BBTriggered.MDC2020ba_best_v1_3.art
