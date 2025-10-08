# Setting up REveMu2e:

## Setting up the code-base

Before beginning the tutorial please do the following:

```
cd /exp/mu2e/app/users/<username>
mkdir Tutorial
cd Tutorial
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups
setup mu2e
setup muse
muse backing TutorialBacking v03_00_01
muse setup
```

## Making the GDML

REveMu2e uses two GDML geometries to import all geometry objects: mu2e.gdml (nominal) and extracted.gdml (for extracted). These can be built by running:

```
muse build GDML
```

## Making your .rootrc

To make a .rootrc, first enter your working directory. Then run:

```
bash config/makerootrc.sh
```
You should now be ready to run the tutorial.

## Choosing a port number

Open the .rootrc in your favorite text editor. Change the numbers 01234 to 0ABCD (where ABCD are the last 4 digits of your Fermilab ID badge).

## Open an ssh tunnel

In your local computer open a ssh tunnel:

```
ssh -KXY -L 0ABCD:localhost:0ABCD <USER>@mu2egpvm0Z.fnal.gov
```

replace ABCD with your numbers inserted into the .rootrc and Z should be the Mu2e machine number you intend to launch the display from.
