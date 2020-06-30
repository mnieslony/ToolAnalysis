# BuildTrainingData configuration toolchain

***********************
# ToolChain for BuildTrainingData
**********************

The recommended sequence of tools to generate a training sample for machine learning classifiers is the following:

* LoadGeometry
* LoadWCSim
* LoadWCSimLAPPD
* MCRecoEventLoader
* DigitBuilder
* HitCleaner
* EventSelector
* TimeClustering
* ParticleIDPDF

The `MCRecoEventLoader` selects only primary muons/electrons and loads Monte Carlo truth variables into the RecoEvent store, whereas the `EventSelector` tool selects specific event types afterwards. The `HitCleaner` tool groups digits into clusters and then passes a flag to the digits indicating whether they have been grouped into a cluster or not. The `TimeClustering` tool selects time clusters within the MRD by requiring a time coincidence cut in a minimum number of MRD channels. At the end of the toolchain, the information from all the tools is fed into the `ParticleIDPDF` tool. `ParticleIDPDF` then calculates specific properties of events that will later be used to classify events (PID / Multiple Ring Rejection).
