#ifndef FindTrackLengthInWater_H
#define FindTrackLengthInWater_H

#include <string>
#include <iostream>
#include <sstream>
#include "ANNIEalgorithms.h"

#include "Tool.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "ExampleRoot.h"

class FindTrackLengthInWater: public Tool {


 public:

  FindTrackLengthInWater();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  double find_lambda(double xmu_rec,double ymu_rec,double zmu_rec,double xrecDir,double yrecDir,double zrecDir,double x_pmtpos,double y_pmtpos,double z_pmtpos,double theta_cher);
  bool Finalise();


 private:
  TTree* regTree;
  TTree * nu_eneNEW;
  ExampleRoot* Data;

  int maxhits0=1100;
  long currententry;
  long NumEvents;
  bool first=1; bool deny_access=0;
  double diffDirAbs2=0; double diffDirAbs=0;
  double recoDWallR2=0; double recoDWallZ2=0;
  int count1=0;
  
  std::ofstream csvfile;
  std::string myfile;
  std::string outputdir="";
  bool writefile=false;
  TFile* outputFile;
  
  // #####
  Geometry* anniegeom=nullptr;
  int totalPMTs;
  int totalLAPPDs;
  double tank_radius;
  double tank_halfheight;
  
  // Variables to be passed to DNNFindTrackLenghInWater
  std::string TrackLengthTrainingDataFile;
  std::string TrackLengthTestingDataFile;
  std::string TrackLengthCheckpointDir;
  std::string DNNTrackLengthPredictionsFile;
  
	// verbosity levels: if 'verbosity' < this level, the message type will be logged.
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	std::string logmessage;
	int get_ok;
};


#endif
