# SelectMuonSample

SelectMuonSample is a supplementary tool to create PDFs for the ANNIE Particle Identification process. There are 8 pdfs in total that are created, i.e. the photoelectrons and time pdfs for both electrons and muons. Separate pdfs are created for LAPPDs. The 3D pdfs have the variables 
* `distance PMT - interaction vertex`
* `angle PMT - particle track`
* `energy of primary particle`

## Data

SelectMuonSample takes the data obtained from the `LoadWCSim` tool and creates the pdfs based on those true values.

**MCParticles** `vector<MCParticle>`
* Takes the MC particle events from the `ANNIEEvent`store and extracts respective `MCHits`

**MCHits** `map<ChannelKey,vector<Hit>>`
* Takes the MC PMT and LAPPD hits from all `MCParticles` in the simulated electron and muon files and extracts the necessary variables to fill into PDFs


## Configuration

SelectMuonSample has the following configuration variables:

```
bins_distance value1
bins_angle value2
bins_energy value3
```

The configuration variables define how many bins the created PDFs will have in x-, y- and z-direction.
