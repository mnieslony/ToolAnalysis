# PlotMCProperties configuration toolchain

***********************
# ToolChain for PlotMCProperties
**********************

The `PlotMCProperties` tool has the purpose to plot basic properties of simulated Monte Carlo files to ensure that their properties are as expected.

The recommended sequence of tools to obtain an overview ROOT file showing the properties of a simulation file is the following: 

* LoadWCSim
* LoadWCSimLAPPD
* MCRecoEventLoader
* DigitBuilder
* EventSelector
* TimeClustering
* PlotMCProperties

All tools executed before `PlotMCProperties` are operated to extract additional information about the simulated events (e.g. number of MRD clusters, event origin in FV/PMTVol, stop in the MRD) that are not available directly from the `MCParticle` class.

