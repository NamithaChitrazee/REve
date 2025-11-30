# Getting Started

## Building GDML

Before attempting to run:

```muse build GDML```

without this command you will see error messages like:

```
  what():  ---- FatalRootError BEGIN
  Fatal Root Error: TFile::Write
  file /dev/null not opened in write mode
  ROOT severity: 2000
---- FatalRootError END
```

## Setting up rootrc

To make a .rootrc in your working directory:

```EventDisplay/config/makerootrc.sh```

## Correcting hangs

If the event display hangs upon startup you can trying killing any stale browser calls using this:

```EventDisplay/config/kill.sh```


