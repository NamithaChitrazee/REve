# Event Navigation

During an analysis the user will want to navigate to a given event. There are a number of ways to do this.

## Default

The default mode has `seqMode == true` meaning the display starts at the first event and the user can sequentially navigate to the next using the GUI.

## Autoplay

The "autoplay" checkbox will override the GUI navigation and simple stream all events until the file end.

## Go to an event: GUI

Text entry boxes allow the user to select to go to any event using just the run and event number. The user must do the following to allow these to function:

* Load in seqMode
* Go to second event using NextEvent()
* Enter chosen run/event and select `Go()` followed by `NextEvent()`

This is not the recommended means to navigate and is currently a placeholder.


## Go to an event: Scripts

The is a custom script in this repo ```config/EventDisplay.sh``` the usage is as follows:

```
./EventDisplay.sh --run 1201 --subrun 34 --event 15028 --dataset mcs.mu2e.ensembleMDS2cMix1BBTriggered.MDC2020ba_best_v1_3.art
```

if you know the name of the dataset.

If you are working with an ntuple, you may not know all the commands to figure out its parent mcs. In this case run:

```
./EventDisplay.sh  --run 1201 --subrun 476 --event 1  --dataset nts.mu2e.ensembleMDS2cMix1BBTriggered.MDC2020ba_best_v1_3_v06_06_00.001201_00000476.root
```

where the run, subrun and event numbers are identified from your analysis to be an event of interest in that root file.
