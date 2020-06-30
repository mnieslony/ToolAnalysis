#ifndef PID_H
#define PID_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "HighEReco.h"
#include "TFile.h"
#include "TH3.h"
#include "TH1.h"
#include "TimeClass.h"


class PID: public Tool {

 public:

  PID();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

  int verbosity = 0;
  std::string inputfile;
  unsigned long NumEvents;
  int loop;
  int time_cut_veto;
  int n_electrons, n_muons, n_electrons_cuts, n_muons_cuts, n_reco_electrons, n_reco_muons;
  bool draw_debug;
  bool parametric;

  //geometry variables
  std::string MCFile;
  uint32_t RunNumber;
  uint32_t SubRunNumber;
  uint32_t EventNumber;
  uint64_t MCEventNum;
  uint16_t MCTriggernum;
  int n_tank_pmts, n_veto_pmts, n_mrd_pmts, n_lappds;
  double tank_center_x, tank_center_y, tank_center_z;
  Geometry *geom = nullptr;
  std::map<int, double> pmts_x, pmts_y, pmts_z;
  std::map<int, double> pmts_orientation_x, pmts_orientation_y, pmts_orientation_z;
  bool is_electron;

  //mc variables
  std::vector<MCParticle>* MCParticles=nullptr;
  std::map<unsigned long,std::vector<Hit>>* MCHits=nullptr;
  std::map<unsigned long,std::vector<LAPPDHit>>* MCLAPPDHits=nullptr;
  std::map<unsigned long,std::vector<Hit>>* TDCData=nullptr;
  TimeClass *EventTime = nullptr;
  std::map<unsigned long,int> detectorkey_to_lappdid;
  std::map<unsigned long,int> channelkey_to_pmtid;
  std::map<unsigned long,int> channelkey_to_mrdpmtid;
  std::map<unsigned long,int> channelkey_to_faccpmtid;
  std::map<int, unsigned long> pmt_tubeid_to_channelkey;
  static const unsigned long n_channels_per_lappd = 60;

  //pid variables
  int is_simulation;		//different read out because different store structure reco/mctruth
  bool is_muon;
  std::string pdf_filename_e, pdf_filename_mu;
  std::string out_filename, out_filename_text;
  double likelihood_value_electron, likelihood_value_muon, likelihood_value_electron_pe, likelihood_value_electron_time, likelihood_value_muon_pe, likelihood_value_muon_time;
  double vtx_x, vtx_y, vtx_z;
  double vtx_stop_x, vtx_stop_y, vtx_stop_z;
  double vtx_dir_x, vtx_dir_y, vtx_dir_z;
  double dir_x, dir_y, dir_z, dir_theta, dir_phi;
  double dir_stop_x, dir_stop_y, dir_stop_z;
  double primary_energy;
  double time_true;
  double t_reco;
  int primary_pid;
  double pmt_charge[200], pmt_time[200];
  double total_pe;
  std::vector<double> UnUsedRho, UnUsedY, UnUsedT, UnUsedEvNum, UnUsedMCEvNum, UnUsedNHits, UnUsedEnergy;

  //pid histos
  TH3F *electronPhotons, *muonPhotons, *electronTimes, *muonTimes;
  TH3F *electronPhotons_lappd, *muonPhotons_lappd, *electronTimes_lappd, *muonTimes_lappd;
  TH1F *likelihood_electron, *likelihood_muon, *likelihood_electron_pdfE, *likelihood_muon_pdfE, *likelihood_electron_pdfMu, *likelihood_muon_pdfMu;
  TH1F *likelihood_electron_pe, *likelihood_electron_time, *likelihood_muon_pe, *likelihood_muon_time;
  TH2F *likelihood_nhits;
  TFile *f_pdf_e, *f_pdf_mu, *file_out;

  //checking variables/histograms
  TH1F *obs_pe, *exp_pe, *diff_pe, *obs_t, *exp_t, *diff_t;
  TApplication *app_pid;

  //status variables (verbosity, bool for file checking)
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  std::string logmessage;
  int get_ok;

};


#endif
