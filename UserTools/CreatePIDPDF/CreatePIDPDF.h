#ifndef CreatePIDPDF_H
#define CreatePIDPDF_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TRandom3.h"

#include "Hit.h"
#include "LAPPDHit.h"
#include "TimeClass.h"
#include "Position.h"
#include "Direction.h"

class CreatePIDPDF: public Tool {

 public:

  CreatePIDPDF();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

  static const int num_energy = 10;
  static const int num_distance = 20;
  static const int num_angle = 100;

 private:

  int verbosity=0;
  std::string inputfile;
  unsigned long NumEvents;
  bool parametric;

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
  bool is_electron;
  double tank_radius, detector_h;
  double inner_radius = 1.275;		//values for inner structure, outer PMT mountings
  double y_min = 0.;
  double y_max = 0.;

  std::map<unsigned long,int> detectorkey_to_lappdid;
  std::map<unsigned long,int> channelkey_to_pmtid;
  std::map<int,unsigned long> pmt_tubeid_to_channelkey;
  std::map<unsigned long,int> channelkey_to_mrdpmtid;
  std::map<unsigned long,int> channelkey_to_faccpmtid;
  static const unsigned long n_channels_per_lappd = 60;

  std::vector<MCParticle>* MCParticles=nullptr;
  std::map<unsigned long,std::vector<Hit>>* MCHits=nullptr;
  std::map<unsigned long,std::vector<LAPPDHit>>* MCLAPPDHits=nullptr;
  std::map<unsigned long,std::vector<Hit>>* TDCData=nullptr;
  TimeClass *EventTime = nullptr;

  double vtx_x, vtx_y, vtx_z, vtx_rho;
  double dir_x, dir_y, dir_z, dir_theta, dir_phi;
  double distWallVert, distWallHor, distInnerStrVert, distInnerStrHor;
  double total_charge, nhits, barycenter_angle, average_time, average_dist;
  double nhits_lappd, average_time_lappd, barycenter_angle_lappd;
  double muonenergy, electronenergy;
  double time_true, t_reco;
  TRandom3 frand; 

  TH3F *h_electron_vertex;
  TH3F *h_electron_dir;
  TH3F *muons_time;
  TH3F *muons_pe;
  TH3F *muons_time_lappd;
  TH3F *muons_pe_lappd;
  TH3F *n_muons;
  TH3F *n_muons_total;
  TH3F *n_muons_lappd;
  TH3F *electrons_time;
  TH3F *electrons_pe;
  TH3F *electrons_time_lappd;
  TH3F *electrons_pe_lappd;
  TH3F *n_electrons;
  TH3F *n_electrons_total;
  TH3F *n_electrons_lappd;
  TH3F *muons_prob_unhit;
  TH3F *electrons_prob_unhit;
  TH1F *log_tubenr;
 
  TH1F *all_theta;
  TH1F *all_npmts;
  TH1F *all_pe;
  TH1F *all_pe_total;
  TH1F *all_time;
  TH1F *all_theta_bary;
  TH1F *all_theta_rms;
  TH1F *all_theta_var;
  TH1F *all_theta_skew;
  TH1F *all_theta_kurt;
  TH1F *all_ratio_ring;
  TH1F *all_ratio_downstream;
  TH1F *all_ratio_ring_noweight;
  TH1F *all_theta_lappd;
  TH1F *all_pe_total_lappd;
  TH1F *all_time_lappd;
  TH1F *all_theta_bary_lappd;
  TH1F *all_theta_rms_lappd;
  TH1F *all_theta_var_lappd;
  TH1F *all_theta_skew_lappd;
  TH1F *all_theta_kurt_lappd;
  TH1F *all_ratio_ring_lappd;
  TH1F *all_ratio_downstream_lappd;
  TH1F *all_pmts_mrd;
  

  //setup vector structure to fill TH3D Likelihoof pdfs in the end
  
  double pe_electrons[num_energy][num_distance][num_angle]={{{0}}};
  double times_electrons[num_energy][num_distance][num_angle]={{{0}}};
  long n_pe_electrons[num_energy][num_distance][num_angle]={{{0}}};
  long n_pe_electrons_total[num_energy][num_distance][num_angle]={{{0}}};
  long n_times_electrons[num_energy][num_distance][num_angle]={{{0}}};
  double pe_muons[num_energy][num_distance][num_angle]={{{0}}};
  double times_muons[num_energy][num_distance][num_angle]={{{0}}};
  long n_pe_muons[num_energy][num_distance][num_angle]={{{0}}};
  long n_pe_muons_total[num_energy][num_distance][num_angle]={{{0}}};
  long n_times_muons[num_energy][num_distance][num_angle]={{{0}}};
  double pe_electrons_lappd[num_energy][num_distance][num_angle]={{{0}}};
  double times_electrons_lappd[num_energy][num_distance][num_angle]={{{0}}};
  long n_pe_electrons_lappd[num_energy][num_distance][num_angle]={{{0}}};
  long n_times_electrons_lappd[num_energy][num_distance][num_angle]={{{0}}};
  double pe_muons_lappd[num_energy][num_distance][num_angle]={{{0}}};
  double times_muons_lappd[num_energy][num_distance][num_angle]={{{0}}};
  long n_pe_muons_lappd[num_energy][num_distance][num_angle]={{{0}}};
  long n_times_muons_lappd[num_energy][num_distance][num_angle]={{{0}}};
  double distance_step, angle_step, energy_step;

  //vectors to save single event information
  std::vector<std::vector<double>> singleEvent_distance, singleEvent_angle, singleEvent_charge, singleEvent_time; 
  std::vector<std::vector<int>> singleEvent_id;
  std::vector<double> singleEvent_distHor, singleEvent_distVert, singleEvent_dist, singleEvent_energy, singleEvent_nhits, singleEvent_total_charge, singleEvent_average_time, singleEvent_average_distance, singleEvent_barycenter, singleEvent_rms, singleEvent_variance, singleEvent_skewness, singleEvent_kurtosis, singleEvent_distInnerStrVert, singleEvent_distInnerStrHor, singleEvent_ratioChargeRing, singleEvent_rms_bary, singleEvent_ratioChargeBary;
  std::vector<double> singleEvent_nhits_lappd, singleEvent_average_time_lappd, singleEvent_barycenter_lappd, singleEvent_rms_lappd, singleEvent_variance_lappd, singleEvent_skewness_lappd, singleEvent_kurtosis_lappd, singleEvent_ratioChargeRing_lappd, singleEvent_rms_bary_lappd, singleEvent_ratioChargeBary_lappd;
  std::vector<int> singleEvent_num_mrd_hits;
  std::vector<int> singleEvent_pid, singleEvent_evnum;
  std::vector<double> singleEvent_highestQ;
  std::vector<double> singleEvent_vtxX, singleEvent_vtxY, singleEvent_vtxZ, singleEvent_dirX, singleEvent_dirY, singleEvent_dirZ;

  //arrays of histograms to save the underlying information for the pdfs --> quite verbose, for debugging purposes
  TH1F* charge_electrons[num_energy][num_distance][num_angle];
  TH1F* time_electrons[num_energy][num_distance][num_angle];
  TH1F* charge_muons[num_energy][num_distance][num_angle];
  TH1F* time_muons[num_energy][num_distance][num_angle];

  TFile *out_file_pdf;			//separate output file for the generated likelihood pdfs
  TFile *out_file_root;
  std::string filename_pdf;
  std::string filename_csv;
  std::string filename_csv_average;
  std::string filename_csv_data;
  std::string filename_root;
  std::string filename_root_overview;

  double cherenkov_angle = 0.719889;	//arccos(1/1.33)

  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  std::string logmessage;
  int get_ok;

};


#endif
