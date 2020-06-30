# ParticleIDPDF

ParticleIDPDF stores classification-driven event properties in a csv-file to enable the training of ML classifier algorithms on this dataset.

## Data

ParticleIDPDF creates an output csv-file where all the properties are saved, delimited by commas (every event is a separate line). It is also possible to save overview plots of the written variables to a ROOT-file for manual inspection of the chosen variables.

## Configuration

ParticleIDPDF uses the following configuration variables

```
verbosity 1                   #verbosity setting of the tool
OutputFile filename           #filename where everything will be saved (.csv/.root will be appended to the name)
UseMCTruth 0                  #use additional true information from the Monte Carlo (not useful for real data)
DrawOverviewPlots 1           #save overview histogram plots of the calculated variables to enable a manual comparison of the variable properties for different event types
```
