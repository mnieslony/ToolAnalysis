#ifndef RunValidationStability_H
#define RunValidationStability_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TSystem.h"
#include "TFitResult.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TCanvas.h"

/**
 * \class RunValidationStability
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
* Contact: b.richards@qmul.ac.uk
*/
class RunValidationStability: public Tool {


 public:

  RunValidationStability(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.


 private:


   int verbosity;
   std::string input_directory;
   std::string output_file;
   int current_run = -1;
   int number_runs = 50;


};


#endif
