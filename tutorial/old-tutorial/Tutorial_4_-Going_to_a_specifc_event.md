# Tutorial_4- Going to a specific event

Often during an analysis you may need to visualize one specific troublesome or interesting event. You may not want to skip sequentially through all the previous events.

## Command line

A simple command line argument can be used if you know which file the event you want is in:

```
mu2e -c nominal_example.fcl /pnfs/mu2e/tape/phy-sim/mcs/mu2e/CeEndpointMix1BBSignal/MDC2020z_best_v1_1/art/01/94/mcs.mu2e.CeEndpointMix1BBSignal.MDC2020z_best_v1_1.001210_00000787.art -e 1210:787:8
```

This will take you to run 1210, sub-run 787 and event number 8.

## Using and event filter

The user can select a specific event using the examples labeled ``select_events_<config>_exampl.fcl" where the <config> can be either nominal or extracted.

Edit the select_events_nominal_example.fcl to look for the event detailed in the above section.

## on the GUI

The GUI input is currently under development....
