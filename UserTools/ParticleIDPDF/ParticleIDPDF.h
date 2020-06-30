#ifndef ParticleIDPDF_H
#define ParticleIDPDF_H

#include <string>
#include <iostream>
#include <vector>
#include <map>

#include "Tool.h"
#include "TH1F.h"
#include "TFile.h"
#include "TMath.h"

#include "Hit.h"
#include "LAPPDHit.h"
#include "Position.h"
#include "Direction.h"
#include "RecoVertex.h"
#include "RecoDigit.h"


class ParticleIDPDF: public Tool {

 public:

  ParticleIDPDF();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

  // Configuration variables
  int verbosity=0;
  std::string filename;
  bool use_mctruth;
  bool draw_overview;

  // ANNIEEvent / RecoStore variables
  int evnum, mcevnum;
  std::map<unsigned long,std::vector<Hit>>* MCHits=nullptr;
  std::map<unsigned long,std::vector<LAPPDHit>>* MCLAPPDHits=nullptr;
  std::map<unsigned long,std::vector<Hit>>* TDCData=nullptr;
  RecoVertex *TrueVertex = nullptr;
  RecoVertex *TrueStopVertex = nullptr;
  bool EventCutStatus;
  std::vector<RecoDigit> RecoDigits; 
  double TrueMuonEnergy;
  int NumMrdTimeClusters;
  int nrings;
  bool no_pik; 

  // Geometry variables
  Geometry *geom = nullptr;
  int n_tank_pmts, n_veto_pmts, n_mrd_pmts, n_lappds;
  double tank_center_x, tank_center_y, tank_center_z;
  std::map<unsigned long, double> pmts_x, pmts_y, pmts_z;
  double tank_R, tank_H;
  double tank_innerR = 1.275;		//values for inner structure, outer PMT mountings
  double tank_ymin = 0.;
  double tank_ymax = 0.;

  // CSV variables control histograms (1D)
  
  TH1F *hist_pmtPE = nullptr;
  TH1F *hist_pmtTime = nullptr;
  TH1F *hist_pmtAngle = nullptr;
  TH1F *hist_pmtAngle2 = nullptr;
  TH1F *hist_pmtPhi = nullptr;
  TH1F *hist_pmtY = nullptr;
  TH1F *hist_pmtDist = nullptr;
  TH1F *hist_pmtHits = nullptr;
  TH1F *hist_pmtPEtotal = nullptr;
  TH1F *hist_pmtAvgTime = nullptr;
  TH1F *hist_pmtAvgDist = nullptr;
  TH1F *hist_pmtAngleBary = nullptr;
  TH1F *hist_pmtAngleRMS = nullptr;
  TH1F *hist_pmtAngleVar = nullptr;
  TH1F *hist_pmtAngleSkew = nullptr;
  TH1F *hist_pmtAngleKurt = nullptr;
  TH1F *hist_pmtBaryRMS = nullptr;
  TH1F *hist_pmtBaryVar = nullptr;
  TH1F *hist_pmtBarySkew = nullptr;
  TH1F *hist_pmtBaryKurt = nullptr;
  TH1F *hist_pmtBaryFracLargeAngle = nullptr;
  TH1F *hist_pmtFracRing = nullptr;
  TH1F *hist_pmtFracDownstream = nullptr;
  TH1F *hist_pmtFracRingNoWeight = nullptr;
  TH1F *hist_pmtFracHighestQ = nullptr;
  TH1F *hist_pmtFracClustered = nullptr;
  TH1F *hist_pmtFracLowCharge = nullptr;
  TH1F *hist_pmtFracLateTime = nullptr;
  TH1F *hist_pmtFracEarlyTime = nullptr;
  TH1F *hist_pmtAngleBary_all = nullptr;
  TH1F *hist_pmtAngleBary_all_ChargeWeighted = nullptr;
  TH1F *hist_pmtAngle2Bary_all = nullptr;
  TH1F *hist_pmtAngle2Bary_all_ChargeWeighted = nullptr;
  TH1F *hist_pmtYBary_all = nullptr;
  TH1F *hist_pmtYBary_all_ChargeWeighted = nullptr;
  TH1F *hist_pmtPhiBary_all = nullptr;
  TH1F *hist_pmtPhiBary_all_ChargeWeighted = nullptr;
  TH1F *hist_pmtPhiBaryRMS = nullptr;
  TH1F *hist_pmtPhiBaryVar = nullptr;
  TH1F *hist_pmtPhiBaryFracLargeAngle = nullptr;

  TH1F *hist_lappdPE = nullptr;
  TH1F *hist_lappdTime = nullptr;
  TH1F *hist_lappdAngle = nullptr;
  TH1F *hist_lappdDist = nullptr;
  TH1F *hist_lappdHits = nullptr;
  TH1F *hist_lappdPEtotal = nullptr;
  TH1F *hist_lappdAvgTime = nullptr;
  TH1F *hist_lappdAvgDist = nullptr;
  TH1F *hist_lappdAngleBary = nullptr;
  TH1F *hist_lappdAngleRMS = nullptr;
  TH1F *hist_lappdAngleVar = nullptr;
  TH1F *hist_lappdAngleSkew = nullptr;
  TH1F *hist_lappdAngleKurt = nullptr;
  TH1F *hist_lappdBaryRMS = nullptr;
  TH1F *hist_lappdBaryVar = nullptr;
  TH1F *hist_lappdBarySkew = nullptr;
  TH1F *hist_lappdBaryKurt = nullptr;
  TH1F *hist_lappdFracRing = nullptr;
  TH1F *hist_lappdFracRingNoWeight = nullptr;
  
  TH1F *hist_mrdPaddles = nullptr;
  TH1F *hist_mrdLayers = nullptr;
  TH1F *hist_mrdconsLayers = nullptr;  
  TH1F *hist_mrdClusters = nullptr;
  TH1F *hist_mrdXSpread = nullptr;
  TH1F *hist_mrdYSpread = nullptr;
  TH1F *hist_mrdAdjacentHits = nullptr;
  TH1F *hist_mrdPaddlesPerLayer = nullptr; 

  TH1F *hist_energy = nullptr;
  TH1F *hist_nrings = nullptr;
  TH1F *hist_multiplerings = nullptr;

  //vectors to save csv event information

  //mctruth variables
  std::vector<double> vec_energy, vec_distHor, vec_distVert, vec_distInnerHor, vec_distInnerVert, vec_pmt_fracRing, vec_pmt_fracRingNoWeight;
  std::vector<double> vec_lappd_fracRing, vec_lappd_fracRingNoWeight; 

  //general variables
  std::vector<double> vec_pmt_avgDist, vec_pmt_hits, vec_pmt_totalQ, vec_pmt_avgT, vec_pmt_baryAngle, vec_pmt_rmsAngle, vec_pmt_varAngle, vec_pmt_skewAngle, vec_pmt_kurtAngle, vec_pmt_rmsBary, vec_pmt_varBary, vec_pmt_skewBary, vec_pmt_kurtBary, vec_pmt_fracDownstream, vec_pmt_fracHighestQ, vec_pmt_highestQ, vec_pmt_fracClustered;
  std::vector<double> vec_lappd_avgDist, vec_lappd_hits, vec_lappd_totalQ, vec_lappd_avgT, vec_lappd_baryAngle, vec_lappd_rmsAngle, vec_lappd_varAngle, vec_lappd_skewAngle, vec_lappd_kurtAngle, vec_lappd_rmsBary, vec_lappd_varBary, vec_lappd_skewBary, vec_lappd_kurtBary; 
  std::vector<int> vec_mrd_paddles, vec_mrd_layers, vec_mrd_conslayers, vec_mrd_clusters;
  std::vector<int> vec_mrd_adjacenthits;
  std::vector<double> vec_mrd_xspread, vec_mrd_yspread, vec_mrd_paddlesperlayer;
  std::vector<int> vec_nrings, vec_multiplerings;
  std::vector<int> vec_evnum, vec_mcevnum, vec_pid;
  std::vector<double> vec_pmt_fracLowQ, vec_pmt_fracLateT, vec_pmt_fracEarlyT;
  std::vector<double> vec_pmt_rmsPhi, vec_pmt_varPhi, vec_pmt_largeanglePhi, vec_pmt_largeangleAngle;


  double cherenkov_angle = 0.719889;	//arccos(1/1.33), if assuming relevant primaries move with velocity c

  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  std::string logmessage;
  int get_ok;

};


#endif
