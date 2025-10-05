# EventDisplay

This repo contains the current EventDisplay source code, example fcls, tutorials and documentation.

The EventDisplay is based on REve which is an updated version of the well-known and widely used TEve ROOT based Event Visualization Software. REve uses modernized infrastcutrue and allows Web Access for remote use.

For more information about REve and the implementation for Mu2e see: https://mu2ewiki.fnal.gov/wiki/Eve7EventDisplay#Examples_of_the_Eve-7_Mu2e_Display.

# Useage

## Sequential on known file
To use the display to scroll through events sequentially utilize either the ```nominal_example``` for nominal Mu2e geometry or ```extracted_example``` for the extracted geometry and run as you would any other .fcl:

```
mu2e -c nominal_example.fcl <filename.art>
```

## Go to an given event

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
# Tutorial

We are in the process of updating our Analysis Tools Tutorial material. In the meantime, please see: https://mu2ewiki.fnal.gov/wiki/EventDisplayTutorial

# Development

This code is still under development by the Mu2e Analysis Tools Group. The development team current consists of Sophie Middleton, Andy Edmonds and Namitha Chithiraseemadam. Please contact us to contribute.

