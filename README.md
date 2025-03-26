# EventDisplay

This repo contains the current EventDisplay source code, example fcls, tutorials and documentation.

The EventDisplay is based on REve which is an updated version of the well-known and widely used TEve ROOT based Event Visualization Software. REve uses modernized infrastcutrue and allows Web Access for remote use.

For more information about REve and the implementation for Mu2e see: https://mu2ewiki.fnal.gov/wiki/Eve7EventDisplay#Examples_of_the_Eve-7_Mu2e_Display.

# Run an Example

You can run the REve Mu2e code in the same way as any analyzer in Mu2e Offline. 

A number of "examples" are contained in the "examples" directory. The nominal_example.fcl is used for displaying reconstructed (and true MC) for nominal geometry:

```
mu2e -c EventDisplay/examples/nominal_example.fcl FILENAME.art  
```

Things can be added and removed on the GUI once set to "true" in the FCL. Just uncheck the boxes associated with the geometry you want to remove.

NOTE: you should run this from your working directory, not inside of EventDisplay. This is for the firefox browser path #FIXME.

For extracted position data/MC please use the extracted_example.fcl

# Go to an Event

Currently we have disabled GUI "go to any event" functionality. To filter events use the fcl files labelled "select_event....fcl".

# The GDML

Recently we upgraded the code to use the newest GDML. To get this built you need to run 

```
muse build GDML

```

Note to make this work you need a local Offline release.

# The Browser

The user is advised to add the following to their .rootrc:

```
WebEve.DisableShow:   0
WebGui.HttpPort:   01234
```
The first line enables the forefox window pop up on the Mu2e machines IGNORE THAT DISPLAY AND USE YOUR LOCAL BROWSER. The second line sets port number to 1234, this can then be used in your local machine. Assuming you have access to the Fermilab network (either on site or using VPN):

```
mu2egpvm01.fnal.gov:01234/win1/

```

There is a script in the config directory which can help create the .rootrc.

Or create a ssh tunnel:
```
ssh  -Y  -A -L 01234:localhost:01234 <username>@mu2egpvm01.fnal.gov


```

You can then view the display on your local machine without need to use a VNC. Be sure to kill the process on the Mu2e machines when finished to avoid future issues.

# Development

This code is still under development by the Mu2e Analysis Tools Group. The development team current consists of Sophie Middleton, Andy Edmonds and Namitha Chithiraseemadam. Please contact us to contribute.
