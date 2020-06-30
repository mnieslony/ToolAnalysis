#ifndef SelectMuonSample_H
#define SelectMuonSample_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TFile.h"

#include "Hit.h"
#include "LAPPDHit.h"
#include "TimeClass.h"
#include "Position.h"
#include "Direction.h"

class SelectMuonSample: public Tool {

 public:

  SelectMuonSample();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();


 private:

  int verbosity=0;
  std::string inputfile;
  unsigned long NumEvents;

  std::string MCFile;
  uint32_t RunNumber;
  uint32_t SubRunNumber;
  uint32_t EventNumber;
  uint64_t MCEventNum;
  uint16_t MCTriggernum;
  int n_tank_pmts, n_veto_pmts, n_mrd_pmts;
  double tank_center_x, tank_center_y, tank_center_z;
  Geometry geom;
  std::map<int, double> pmts_x, pmts_y, pmts_z;



  std::vector<MCParticle>* MCParticles=nullptr;
  std::map<ChannelKey,std::vector<Hit>>* MCHits=nullptr;
  std::map<ChannelKey,std::vector<LAPPDHit>>* MCLAPPDHits=nullptr;
  std::map<ChannelKey,std::vector<Hit>>* TDCData=nullptr;
  TimeClass *EventTime = nullptr;

  TH1F *h_muon_energy;
  TH1F *h_muon_pe;
  TH1F *h_muon_lappd_pe;
  TH3F *h_muon_vertex;
  TH3F *h_muon_dir;
  TH1F *h_electron_energy;
  TH1F *h_electron_pe;
  TH1F *h_electron_lappd_pe;
  TH3F *h_electron_vertex;
  TH3F *h_electron_dir;
  TH3F *muons_time;
  TH3F *muons_pe;
  TH3F *muons_time_lappd;
  TH3F *muons_pe_lappd;
  TH3F *electrons_time;
  TH3F *electrons_pe;
  TH3F *electrons_time_lappd;
  TH3F *electrons_pe_lappd;


  TTree *tree_beam;
  double mu_energy, mu_pe, mu_pe_lappd;
  double mu_pos[3],mu_stop_pos[3],mu_dir[3],mu_start_stop[3],mu_stop_dir[3];
  int is_electron;
  std::vector<double> time_pmt, time_lappd, mu_angle_lappd, mu_distance_lappd, mu_angle_single, mu_distance_single;
  std::vector<double> mu_pe_single;
  std::vector<int> mu_pmt_single, mu_lappd;
  double vtx_x, vtx_y, vtx_z;
  double vtx_stop_x, vtx_stop_y, vtx_stop_z;
  double vtx_dir_x, vtx_dir_y, vtx_dir_z;
  double dir_x, dir_y, dir_z;
  double dir_stop_x, dir_stop_y, dir_stop_z;
  double muonenergy, electronenergy;

  // setup vector structure to fill TH3D Likelihoof pdfs in the end
  double pe_electrons[21][100][100]={{{0}}};
  double times_electrons[21][100][100]={{{0}}};
  long n_pe_electrons[21][100][100]={{{0}}};
  long n_times_electrons[21][100][100]={{{0}}};
  double pe_muons[21][100][100]={{{0}}};
  double times_muons[21][100][100]={{{0}}};
  long n_pe_muons[21][100][100]={{{0}}};
  long n_times_muons[21][100][100]={{{0}}};
  double pe_electrons_lappd[21][100][100]={{{0}}};
  double times_electrons_lappd[21][100][100]={{{0}}};
  long n_pe_electrons_lappd[21][100][100]={{{0}}};
  long n_times_electrons_lappd[21][100][100]={{{0}}};
  double pe_muons_lappd[21][100][100]={{{0}}};
  double times_muons_lappd[21][100][100]={{{0}}};
  long n_pe_muons_lappd[21][100][100]={{{0}}};
  long n_times_muons_lappd[21][100][100]={{{0}}};
  
  double distance_step, angle_step, energy_step;

  TFile *out_file;
  std::string filename;
  TFile *out_file_pdf;			//separate output file for the generated likelihood pdfs
  std::string filename_pdf;

  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  std::string logmessage;
  int get_ok;

};


#endif
