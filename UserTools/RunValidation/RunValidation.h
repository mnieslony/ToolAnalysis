#ifndef RunValidation_H
#define RunValidation_H

#include <string>
#include <iostream>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TROOT.h"
#include "TCanvas.h"

#include "Tool.h"


/**
 * \class RunValidation
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $Author: M.Nieslony $
* $Date: 2020/03/02 10:44:00 $
* Contact: mnieslon@uni-mainz.de
*/
class RunValidation: public Tool {


 public:

  RunValidation(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.

  void DefineHistograms();  ///< Define histograms for RunValidation tool

 private:

  //configuration variables
  int verbosity;
  std::string outfile_path;
  std::string filename_suffix;
  bool invert_mrd_times;
  std::string singlePEgains;
  int user_runnumber;
  int user_subrunnumber;
  int user_runtype;
  bool createImages;
  std::string savePath; 

  //Data storing variables
  std::map<double,std::vector<Hit>>* m_all_clusters;  //from ClusterFinder tool
  std::map<double,std::vector<unsigned long>>* m_all_clusters_detkey;  //from ClusterFinder tool
  std::map<double,double> ClusterChargeBalances;  //from ClusterClassifiers tool
  std::vector<std::vector<int>> MrdTimeClusters;  //from TimeClustering tool
  std::vector<double> MrdDigitTimes;  //from TimeClustering tool
  std::vector<unsigned long> mrddigitchankeysthisevent;  //from TimeClustering tool
  std::map<unsigned long,double> pmt_gains; //from txt file
  std::map<unsigned long,std::vector<Hit>>* TDCData=nullptr;
  Geometry *geom = nullptr;

  //Counter variables to get rough rates
  int n_entries;
  int n_pmt_clusters;
  int n_pmt_clusters_threshold;
  int n_mrd_clusters;
  int n_pmt_mrd_clusters;
  int n_pmt_mrd_nofacc;
  int n_facc;
  int n_facc_pmt;
  int n_facc_mrd;
  int n_facc_pmt_mrd;
  bool first_entry;
  int n_pmt_mrd_time;
  int n_pmt_mrd_time_facc;
  int n_pmt_mrd_time_nofacc;

  //RunStartTime variables
  int RunNumber, SubRunNumber, RunType, EventNumber;
  ULong64_t RunStartTime, EventTimeTank;
  ULong64_t start_time;
  ULong64_t current_time;
  int GlobalRunNumber;
  int GlobalSubRunNumber;

  //ROOT-related variables
  TFile *outfile = nullptr;
  TH1D *MRD_t_clusters = nullptr;
  TH1D *MRD_t_clusters_cosmic = nullptr;
  TH1D *PMT_t_clusters = nullptr;
  TH1D *PMT_t_clusters_2pe = nullptr;
  TH1D *PMT_t_clusters_5pe = nullptr;
  TH1D *PMT_t_clusters_10pe = nullptr;
  TH1D *PMT_t_clusters_30pe = nullptr;
  TH1D *PMT_t_clusters_cosmics = nullptr;
  TH1D *PMT_t_clusters_2pe_cosmics = nullptr;
  TH1D *PMT_t_clusters_5pe_cosmics = nullptr;
  TH1D *PMT_t_clusters_10pe_cosmics = nullptr;
  TH1D *PMT_t_clusters_30pe_cosmics = nullptr;
  TH1D *PMT_t_clusters_led = nullptr;
  TH1D *PMT_t_clusters_2pe_led = nullptr;
  TH1D *PMT_t_clusters_5pe_led = nullptr;
  TH1D *PMT_t_clusters_10pe_led = nullptr;
  TH1D *PMT_t_clusters_30pe_led = nullptr;
  TH1D *PMT_t_clusters_full = nullptr;
  TH1D *PMT_t_clusters_2pe_full = nullptr;
  TH1D *PMT_t_clusters_5pe_full = nullptr;
  TH1D *PMT_t_clusters_10pe_full = nullptr;
  TH1D *PMT_t_clusters_30pe_full = nullptr;
  TH2D *MRD_PMT_t = nullptr;
  TH2D *MRD_PMT_t_100pe = nullptr;
  TH1D *MRD_PMT_Deltat = nullptr;
  TH1D *MRD_PMT_Deltat_100pe = nullptr;
  TH2D *MRD_PMT_t_Cosmic = nullptr;
  TH2D *MRD_PMT_t_100pe_Cosmic = nullptr;
  TH1D *MRD_PMT_Deltat_Cosmic = nullptr;
  TH1D *MRD_PMT_Deltat_100pe_Cosmic = nullptr;
  TH1D *FMV_PMT_Deltat = nullptr;
  TH1D *FMV_PMT_Deltat_100pe = nullptr;
  TH1D *MRD_FMV_Deltat = nullptr;
  TH1D *PMT_prompt_charge = nullptr;
  TH1D *PMT_prompt_charge_10hits = nullptr;
  TH1D *PMT_prompt_charge_zoom = nullptr;
  TH1D *PMT_prompt_charge_MRDCoinc = nullptr;
  TH1D *PMT_prompt_charge_FMV = nullptr;
  TH1D *PMT_prompt_charge_Beam = nullptr;
  TH1D *PMT_prompt_charge_Cosmic = nullptr;
  TH1D *PMT_prompt_charge_LED = nullptr;
  TH1D *PMT_prompt_charge_MRDCoinc_NoFMV = nullptr;
  TH2D *PMT_prompt_charge_CB = nullptr;
  TH1D *PMT_chargeperpmt = nullptr;
  TH1D *PMT_chargeperpmt_100pe = nullptr;
  TH1D *PMT_delayed_charge = nullptr;
  TH1D *PMT_delayed_charge_10hits = nullptr;
  TH1D *PMT_delayed_charge_zoom = nullptr;
  TH2D *PMT_delayed_charge_CB = nullptr;
  TH1D *PMT_Channelkeys = nullptr;
  TH1D *MRD_Channelkeys = nullptr;
  TH1D *FMV_Channelkeys = nullptr;
  TH1D *ANNIE_counts = nullptr;
  TH1D *ANNIE_rates = nullptr;
  TH1D *ANNIE_fractions = nullptr;
  TH1D *PMT_DelayedMult = nullptr;
  TH1D *PMT_DelayedMult_Coinc = nullptr;
  TH1D *PMT_DelayedMult_Coinc_NoFMV = nullptr;
  TH1D *PMT_DelayedMult_Coinc_NoFMV_CB = nullptr;
  TH1D *PMT_DelayedTime = nullptr;
  TH1D *PMT_DelayedTime_CB = nullptr;
  TH1D *PMT_DelayedTime_Coinc = nullptr;
  TH1D *PMT_DelayedTime_Coinc_NoFMV = nullptr;
  TH1D *PMT_DelayedTime_Coinc_NoFMV_CB = nullptr;
  TH1D *ADCWaveform_Samples = nullptr;
  TH1D *Triggerwords = nullptr;
 
  //Canvas
  TCanvas *canvas_beamspill = nullptr;
  TCanvas *canvas_mrd_pmt = nullptr;
  TCanvas *canvas_triggerwords = nullptr;
  TCanvas *canvas_prompt_charge = nullptr;
  TCanvas *canvas_neutron_mult = nullptr;
  TCanvas *canvas_channelkey = nullptr;
  TCanvas *canvas_rates = nullptr;

  //verbosity-related variables
  int v_error = 0;
  int v_warning = 1;
  int v_message = 2;
  int v_debug = 3;
  int vv_debug = 4;


};


#endif
