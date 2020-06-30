# PID

The PID tool performs a likelihood fit to distinguish between electron- and muon-like events. The tool needs access to the pdfs created by `SelectMuonSample`, which can be saved e.g. in a root-file called `likelihood_pdfs.root` and should be placed in the `supplementary` subdirectory of ToolAnalysis.

## Data

PID uses the pdfs saved in `likelihood_pdfs.root` and uses the help class `HighEReco` to perform its likelihood fits.

The result of the PID tool is the single variable `bool(muon_like)`, which will be `true` in case the fit results in the muon being the most likely hypothesis and `false` in case the event was determined to be electron-like. The variable is then stored in `ANNIEEvent` as a part of the reconstruction class.

## Configuration

The PID tool has the following configuration values:

```
likelihood_file value1
input_format value2
```

* `value1` should display the file name of the ROOT-file with the likelihood functions.
* `value2` should display the input format of the data (simulation/reco)
