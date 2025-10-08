# Introduction

## Looping through events

To view events sequentially run on a given file from the dataset:

```

```

## Going to specific event

To view a specific event, you will need to establish the parent (generally mcs) dataset.

### ```pyutils```

Write a simple function in pyutils. Extract the run, subrun and event number for an event of interest and call the display using ```pydisplay```.

An example:

```
display = Display(verbosity=2)
display.pick_event(dataset, run,subrun,event)
display.launch_display(dataset, run,subrun,event)
```

### ```rooutils```

Write a simple C++ macro inspect a dataset, print out an run, subrun and event number of an interesting event. Use ```roodisplay.hh``` to help launch the display:

Here's an example:

```
launch_display("mcs.mu2e.ensembleMDS2cMix1BBTriggered.MDC2020ba_best_v1_3.art", 15028,34,1201)
```

### ```EventDisplay.sh```

```
./EventDisplay.sh --run 1201 --subrun 34 --event 15028 --datset mcs.mu2e.ensembleMDS2cMix1BBTriggered.MDC2020ba_best_v1_3.art
```
