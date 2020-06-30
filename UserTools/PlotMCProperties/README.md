# PlotMCProperties

PlotMCProperties is a tool to display basic mc truth variables in a ROOT-based output format. Properties of the simulation file like the direction, energy and positions of primary particles can be accessed via histograms. Furthermore, the detector and particle properties are also stored in a tree to enable a more flexible selection and plotting of data.

## Data

PlotMCProperties takes the properties from the `MCParticles` object. To enable a comparison with the detector response, it also accesses the properties `MCHits` (tank PMTs), `MCLAPPDHits` (tank LAPPDs) and `TDCData` (MRD PMTs).

## Configuration

PlotMCProperties needs the following configuration variables

```
verbose 1                 #the verbosity setting of the tool
OutFile mcproperties.root #define the output root filename
```
