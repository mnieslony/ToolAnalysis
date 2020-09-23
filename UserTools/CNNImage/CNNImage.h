#ifndef CNNImage_H
#define CNNImage_H

#include <string>
#include <fstream>

#include "TH2F.h"
#include "TMath.h"
#include "TFile.h"

#include "Tool.h"


/**
 * \class CNNImage
 *
 * The tool CNNImage is supposed to create custom csv input files for CNN classification processes. The framework relies on the EventDisplay tool and basically feeds that geometric information of the hit pattern into a matrix format
*
* $Author: M.Nieslony $
* $Date: 2019/08/09 10:44:00 $
* Contact: mnieslon@uni-mainz.de
*/

class CNNImage: public Tool {


 public:

  CNNImage(); ///< constructor for CNNImage class
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.

  void draw_cnn_image(); ///< Fill the hit pattern information in a TH2
  void ConvertPositionTo2D(Position xyz_pos, double &x, double &y); ///<Convert 3D position into 2D rolled up coordinates
  void ConvertPositionTo2D_Top(Position xyz_pos, double &x, double &y); ///<Convert 3D position into 2D rolled up coordinates (top PMTs only)
  void ConvertPositionTo2D_Bottom(Position xyz_pos, double &x, double &y); ///<Convert 3D position into 2D rolled up coordinates (bottom PMTs only)


 private:

  //configuration variables
  std::string cnn_outpath; //path where to save the CNN image information
  std::string detector_config;
  std::string data_mode;  //Charge, Time
  std::string save_mode;  //How is the PMT information supposed to be written out? Geometric/PMT-wise
  int dimensionX;        //dimension of the CNN image in x-direction (Geometric mode)
  int dimensionY;        //dimension of the CNN image in y-direction (Geometric mode)
  bool includeTopBottom;
  int verbosity;
  //ANNIEEvent variables
  int runnumber;
  int subrunnumber;
  int evnum;
  TimeClass* EventTime=nullptr;
  std::vector<MCParticle>* mcparticles = nullptr;
  std::map<unsigned long, std::vector<MCHit>>* MCHits = nullptr;
  std::vector<TriggerClass>* TriggerData = nullptr;
  uint16_t MCTriggernum=0;
  TimeClass TriggerTime;
  Geometry *geom = nullptr;

  // RecoEvent
  int nrings;

  //geometry variables
  double tank_radius;
  double tank_height;
  double tank_center_x;
  double tank_center_y;
  double tank_center_z;
  double min_z, max_z;
  int n_tank_pmts;
  double size_top_drawing = 0.1;
  std::vector<double> vec_pmt2D_x, vec_pmt2D_x_Top, vec_pmt2D_y_Top, vec_pmt2D_x_Bottom, vec_pmt2D_y_Bottom, vec_pmt2D_y;
  int npmtsX, npmtsY;
  std::vector<double> phi_positions;

  //I/O variables
  ofstream outfile, outfile_time, outfile_firsttime, outfile_Rings;
  ofstream outfile_abs, outfile_abs_time, outfile_abs_firsttime;
  TFile *file = nullptr;
  TH1F *h_trigger_neutron = nullptr;
  TH1F *h_charge_prompt = nullptr;
  TH1F *h_charge_delayed = nullptr;
  TH1F *h_time_total = nullptr;
  TH1F *h_ly_prompt = nullptr;

  //PMT information
  std::map<int, double> x_pmt, y_pmt, z_pmt;
  std::map<unsigned long,double> charge, time, first_time;
  double maximum_pmts;
  double total_charge_pmts;
  int total_hits_pmts;
  double min_time_pmts, max_time_pmts;
  double max_firsttime_pmts, min_firsttime_pmts;

  //detectorkey layout organization
  std::map<unsigned long, int> channelkey_to_pmtid;
  std::map<int,unsigned long> pmt_tubeid_to_channelkey;
  std::vector<unsigned long> pmt_detkeys, pmt_chankeys;
  std::vector<unsigned long> hitpmt_detkeys;

  //Verbosity levels
  int v_error = 0;
  int v_warning = 1;
  int v_message = 2;
  int v_debug = 3;
  int vv_debug = 4;


};


#endif
