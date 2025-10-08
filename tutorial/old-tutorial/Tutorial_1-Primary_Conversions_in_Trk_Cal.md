# Tutorial 1: Visualizing a primary conversion (with pile-up) in the detectors

In this exercise you will visualize a signal conversion electron event. 

## nominal_example.fcl

To visualize any events within the nominal Mu2e geometry you should use the examples/nominal_example.fcl as you would any .fcl file:

```
mu2e -c nominal_example.fcl <filename>.art
```

## Looking at a signal event

To visualize a conversion electron event (clean, no pile-up) run this command

```
mu2e -c nominal_example.fcl /pnfs/mu2e/tape/phy-sim/mcs/mu2e/CeEndpointMix1BBSignal/MDC2020z_best_v1_1/art/01/94/mcs.mu2e.CeEndpointMix1BBSignal.MDC2020z_best_v1_1.001210_00000787.art
```

You will see the following:

* Tracker, IPA, OPA, Stopping Target, Calorimeter Crystals, CRV
* MC Tracjectories: truth level output from G4.
* KinKal fit result: all fitting tracks within the event are shown
* CaloClusters: calorimeter clusters, size is proportional to energy deposited. Color is proportional to time. Red should coincide with the track.
* TrkStrawHits: all collections of TrkStrawHits which were used in the track fit

## Zoom in and rotate

To look more closely at parts of the event drag your mouse on the 3D view, scroll to zoom in and out.

## Extracting information from the event

Specific information about the features of a data product can be displayed by rolling your mouse over the visualized products. Try to see what's displayd for 1) the clusters, 2) the MC trajectory.

## Printing information to terminal

Select the PrintRecoInfo() button. You will see that in the terminal reco product information is printed. Do the same for PrintMCInfo().

## Going to the next event

To go to the next event simply press the "NextEvent()" button.

## Quiting

To quit simply select the "Quit()" button.
