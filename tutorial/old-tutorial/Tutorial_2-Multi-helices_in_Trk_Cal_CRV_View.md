# Tutorial 2: Visualize multi-track events in the detector regions

## Looking at multiple track events

A nice way to see the displays speed and usefulness is to look at events with multiple tracks. Try the following:

```
mu2e -c nominal_example.fcl /exp/mu2e/data/projects/pbar2m/datasets/pbar1b0s51r0000/mcs.mu2e.pbar1b0s51r0000.pbar2m.001210_00000000.art
```

Again zoom in and extract information.

## Removing a track

If you'd like to see tracks independently (for example if they overlap). Look on the drop down menu on the left. Roll down the EventScene, uncheck the label associated with the track you want to remove.

You should see that track removed from the display.

Add it back by rechecking the label.
