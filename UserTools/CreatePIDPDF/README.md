# CreatePIDPDF

CreatePIDPDF utilizes Monte Carlo files to produce 3D PDFs of the expected behavior of a particle species in the variables
* `distance vertex - PMT`
* `angle between primary direction & PMT`
* `energy of primary particle`

## Data

CreatePIDPDF stores its result in dedicated ROOT-files in the form of TH3F histograms. There are separate PDFs for PMTs and LAPPDs, and separate histograms for storing the time and charge predictions. The output filename can be specified in the config file. The histograms that are present in the ROOT-file are the following:
* `electrons_time`
* `electrons_pe`
* `muons_time`
* `muons_pe`
* `electrons_time_lappd`
* `electrons_pe_lappd`
* `muons_time_lappd`
* `muons_pe_lappd`


## Configuration

Describe any configuration variables for CreatePIDPDF. The given example creates a muon PDF given a muon MC file.

```
OutputFile_PDF /ANNIECode/ToolAnalysis/pdf_muon.root
verbosity 3

```
