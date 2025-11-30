# Event Navigation

## Sequential Navigation from First Event

## Autoplay

## Go to an event: GUI

## Got to an event : FCL

## Go to an event: Scripts

There are a number of ways to go from your analysis to a visualized event. Utilities are provided in both rooutils [cite] and pyutils [https://github.com/Mu2e/pyutils/blob/main/pyutils/pydisplay.py]. There is also a custom script in this repo ```config/EventDisplay.sh``` the usage is as follows:

```
./EventDisplay.sh --run 1201 --subrun 34 --event 15028 --dataset mcs.mu2e.ensembleMDS2cMix1BBTriggered.MDC2020ba_best_v1_3.art
```

if you know the name of the dataset.

If you are working with an ntuple, you may not know all the commands to figure out its parent mcs. In this case run:

```
./EventDisplay.sh  --run 1201 --subrun 476 --event 1  --dataset nts.mu2e.ensembleMDS2cMix1BBTriggered.MDC2020ba_best_v1_3_v06_06_00.001201_00000476.root
```

where the run, subrun and event numbers are identified from your analysis to be an event of interest in that root file.
