#include "RunValidationStability.h"

RunValidationStability::RunValidationStability():Tool(){}


bool RunValidationStability::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  m_variables.Get("verbosity",verbosity);
  m_variables.Get("InputDirectory",input_directory);
  m_variables.Get("OutputFile",output_file);
  m_variables.Get("CurrentRun",current_run);
  m_variables.Get("NumberOfRuns",number_runs);

  return true;
}


bool RunValidationStability::Execute(){

  std::vector<int> run_numbers;
  std::vector<double> rate, rate_pmt, rate_pmtmrd, rate_pmtmrdfmv, rate_pmtmrdnofmv, rate_fmv, rate_mrd;
  std::vector<double> frac, frac_pmt, frac_pmtmrd, frac_pmtmrdfmv, frac_pmtmrdnofmv, frac_fmv, frac_mrd;
  std::vector<double> align_pmt_mrd, align_pmt_fmv, align_mrd_fmv;
  std::vector<double> frac_prompt, frac_delayed, frac_beam, frac_led, frac_cosmic;
  std::vector<double> avg_mult, avg_mult_coinc, avg_mult_coinc_nofmv_cb;

  TH1D *PMT_t_clusters_combined = new TH1D("PMT_t_clusters_combined","PMT cluster times (combined)",250,0,2000);
  TH1D *PMT_t_clusters_2pe_combined = new TH1D("PMT_t_clusters_2pe_combined","PMT cluster times (#geq 2.p.e, combined)",250,0,2000);
  TH1D *PMT_t_clusters_cosmics_combined = new TH1D("PMT_t_clusters_cosmics_combined","PMT cluster times Cosmics (combined)",250,0,2000);
  TH1D *PMT_t_clusters_2pe_cosmics_combined = new TH1D("PMT_t_clusters_2pe_cosmics_combined","PMT cluster times Cosmics (#geq 2.p.e, combined)",250,0,2000);
  TH1D *PMT_t_clusters_led_combined = new TH1D("PMT_t_clusters_led_combined","PMT cluster times LED (combined)",250,0,2000);
  TH1D *PMT_t_clusters_2pe_led_combined = new TH1D("PMT_t_clusters_2pe_led_combined","PMT cluster times LED (#geq 2.p.e, combined)",250,0,2000);
  TH1D *MRD_t_clusters_combined = new TH1D("MRD_t_clusters_combined","MRD cluster times (combined)",250,0,4000);
  TH1D *PMT_DelayedMult_combined = new TH1D("PMT_DelayedMult_combined","Neutron multiplicity (combined)",20,0,20);
  TH1D *PMT_DelayedMult_Coinc_combined = new TH1D("PMT_DelayedMult_Coinc_combined","Neutron multiplicity (PMT/MRD Coincidence, combined)",20,0,20);
  TH1D *PMT_DelayedMult_Coinc_NoFMV_CB_combined = new TH1D("PMT_DelayedMult_Coinc_NoFMV_CB_combined","Neutron multiplicity (PMT/MRD Coincidence, No FMV, CB < 0.4, combined)",20,0,20);
  TH1D *Triggerwords_combined = new TH1D("Triggerwords_combined","Triggerwords (combined)",60,0,60);  
  TH1D *ADCWaveform_Samples_combined = new TH1D("ADCWaveformSamples_combined","ADC Waveform Samples",5000,0,50000);
  TH1D *PMT_DelayedTime_combined = new TH1D("PMT_DelayedTime_combined","Neutron candidate time distribution (combined)",100,12000,67000);
  TH1D *PMT_DelayedTime_CB_combined = new TH1D("PMT_DelayedTime_CB_combined","Neutron candidate time distribution (CB < 0.4,combined)",100,12000,67000);
  TH1D *PMT_DelayedTime_Coinc_combined = new TH1D("PMT_DelayedTime_Coinc_combined","Neutron candidate time distribution (PMT/MRD Coinc, combined)",100,12000,67000);
  TH1D *PMT_DelayedTime_Coinc_NoFMV_combined = new TH1D("PMT_DelayedTime_Coinc_NoFMV_combined","Neutron candidate time distribution (PMT/MRD Coinc, No FMV, combined)",100,12000,67000);
  TH1D *PMT_DelayedTime_Coinc_NoFMV_CB_combined = new TH1D("PMT_DelayedTime_Coinc_NoFMV_CB_combined","Neutron candidate time distribution (PMT/MRD Coinc, No FMV, CB < 0.4, combined)",100,12000,67000);


  for (int i_run = current_run - number_runs; i_run <= current_run; i_run++){
    std::stringstream runvalidation_name;
    runvalidation_name << input_directory << "/RunValidation_R" << i_run << "S0T3.root";
    if (!gSystem->AccessPathName(runvalidation_name.str().c_str())){
      TFile *fval = new TFile(runvalidation_name.str().c_str(),"READ");
      TH1D *PMT_t_clusters = (TH1D*) fval->Get("PMT_t_clusters");
      PMT_t_clusters_combined->Add(PMT_t_clusters);
      TH1D *PMT_t_clusters_2pe = (TH1D*) fval->Get("PMT_t_clusters_2pe");
      PMT_t_clusters_2pe_combined->Add(PMT_t_clusters_2pe);
      TH1D *PMT_t_clusters_cosmics = (TH1D*) fval->Get("PMT_t_clusters_cosmics");
      PMT_t_clusters_cosmics_combined->Add(PMT_t_clusters_cosmics);
      TH1D *PMT_t_clusters_2pe_cosmics = (TH1D*) fval->Get("PMT_t_clusters_2pe_cosmics");
      PMT_t_clusters_2pe_cosmics_combined->Add(PMT_t_clusters_2pe_cosmics);
      TH1D *PMT_t_clusters_led = (TH1D*) fval->Get("PMT_t_clusters_led");
      PMT_t_clusters_led_combined->Add(PMT_t_clusters_led);
      TH1D *PMT_t_clusters_2pe_led = (TH1D*) fval->Get("PMT_t_clusters_2pe_led");
      PMT_t_clusters_2pe_led_combined->Add(PMT_t_clusters_2pe_led);
      TH1D *MRD_t_clusters = (TH1D*) fval->Get("MRD_t_clusters");
      MRD_t_clusters_combined->Add(MRD_t_clusters);
      TH1D *MRD_PMT_Deltat_100pe = (TH1D*) fval->Get("MRD_PMT_Deltat_100pe");
      TH1D *FMV_PMT_Deltat_100pe = (TH1D*) fval->Get("FMV_PMT_Deltat_100pe");
      TH1D* MRD_FMV_Deltat = (TH1D*) fval->Get("MRD_FMV_Deltat");
      TH1D *PMT_DelayedMult = (TH1D*) fval->Get("PMT_DelayedMult");
      PMT_DelayedMult_combined->Add(PMT_DelayedMult);
      TH1D *PMT_DelayedMult_Coinc = (TH1D*) fval->Get("PMT_DelayedMult_Coinc");
      PMT_DelayedMult_Coinc_combined->Add(PMT_DelayedMult_Coinc);
      TH1D *PMT_DelayedMult_Coinc_NoFMV_CB = (TH1D*) fval->Get("PMT_DelayedMult_Coinc");
      PMT_DelayedMult_Coinc_NoFMV_CB_combined->Add(PMT_DelayedMult_Coinc_NoFMV_CB);
      TH1D *ADCWaveform_Samples = (TH1D*) fval->Get("ADCWaveform_Samples");
      ADCWaveform_Samples_combined->Add(ADCWaveform_Samples);
      TH1D *Triggerwords = (TH1D*) fval->Get("Triggerwords");
      Triggerwords_combined->Add(Triggerwords);
      TH1D *ANNIE_rates = (TH1D*) fval->Get("ANNIE_rates");
      TH1D *ANNIE_fractions = (TH1D*) fval->Get("ANNIE_fractions");
      TH1D *PMT_DelayedTime = (TH1D*) fval->Get("PMT_DelayedTime");
      PMT_DelayedTime_combined->Add(PMT_DelayedTime);
      TH1D *PMT_DelayedTime_CB = (TH1D*) fval->Get("PMT_DelayedTime_CB");
      PMT_DelayedTime_CB_combined->Add(PMT_DelayedTime_CB);
      TH1D *PMT_DelayedTime_Coinc = (TH1D*) fval->Get("PMT_DelayedTime_Coinc");
      PMT_DelayedTime_Coinc_combined->Add(PMT_DelayedTime_Coinc);
      TH1D *PMT_DelayedTime_Coinc_NoFMV = (TH1D*) fval->Get("PMT_DelayedTime_Coinc_NoFMV");
      PMT_DelayedTime_Coinc_NoFMV_combined->Add(PMT_DelayedTime_Coinc_NoFMV);
      TH1D *PMT_DelayedTime_Coinc_NoFMV_CB = (TH1D*) fval->Get("PMT_DelayedTime_Coinc_NoFMV_CB");
      PMT_DelayedTime_Coinc_NoFMV_CB_combined->Add(PMT_DelayedTime_Coinc_NoFMV_CB);

      TFitResultPtr ptr = MRD_PMT_Deltat_100pe->Fit("gaus","S","",700,800);
      align_pmt_mrd.push_back(ptr->Parameter(1));
      TFitResultPtr ptr2 = FMV_PMT_Deltat_100pe->Fit("gaus","S","",750,850);
      align_pmt_fmv.push_back(ptr2->Parameter(1));
      TFitResultPtr ptr3 = MRD_FMV_Deltat->Fit("gaus","S","",0,100);
      align_mrd_fmv.push_back(ptr3->Parameter(1));

      double temp_avg_mult=0;
      double temp_avg_mult_coinc=0;
      double temp_avg_mult_coinc_nofmv_cb=0;
      for (int i_bin=0; i_bin < PMT_DelayedMult->GetXaxis()->GetNbins(); i_bin++){
        temp_avg_mult+=(i_bin*PMT_DelayedMult->GetBinContent(i_bin+1));
        temp_avg_mult_coinc+=(i_bin*PMT_DelayedMult_Coinc->GetBinContent(i_bin+1));
        temp_avg_mult_coinc_nofmv_cb+=(i_bin*PMT_DelayedMult_Coinc_NoFMV_CB->GetBinContent(i_bin+1));
      }
      int nentries = PMT_DelayedMult->GetEntries();
      int nentries_coinc = PMT_DelayedMult_Coinc_NoFMV_CB->GetEntries();
      int nentries_coinc_nofmv_cb = PMT_DelayedMult_Coinc_NoFMV_CB->GetEntries();
      if (nentries > 0) temp_avg_mult /= nentries;
      if (nentries_coinc > 0) temp_avg_mult_coinc /= nentries_coinc;
      if (nentries_coinc_nofmv_cb > 0) temp_avg_mult_coinc_nofmv_cb /= nentries_coinc_nofmv_cb;
      avg_mult.push_back(temp_avg_mult);
      avg_mult_coinc.push_back(temp_avg_mult_coinc);
      avg_mult_coinc_nofmv_cb.push_back(temp_avg_mult_coinc_nofmv_cb); 

      double n_beam = Triggerwords->GetBinContent(6);
      double n_led = Triggerwords->GetBinContent(32);
      double n_cosmic = Triggerwords->GetBinContent(37);
      double n_total = n_beam+n_led+n_cosmic;
      if (n_total>0){
        n_beam/=n_total;
        n_led/=n_total;
        n_cosmic/=n_total;
      }
      frac_beam.push_back(n_beam);
      frac_led.push_back(n_led);
      frac_cosmic.push_back(n_cosmic);
 
      double n_prompt = ADCWaveform_Samples->GetBinContent(100);
      double n_delayed = ADCWaveform_Samples->GetBinContent(3500);
      double n_prompt_delayed = n_prompt+n_delayed;
      if (n_prompt_delayed){
        n_prompt /= n_prompt_delayed;
        n_delayed /= n_prompt_delayed;
      }
      frac_prompt.push_back(n_prompt);
      frac_delayed.push_back(n_delayed);

      rate.push_back(ANNIE_rates->GetBinContent(1));
      rate_pmt.push_back(ANNIE_rates->GetBinContent(2));
      rate_pmtmrd.push_back(ANNIE_rates->GetBinContent(5));
      rate_pmtmrdfmv.push_back(ANNIE_rates->GetBinContent(10));
      rate_pmtmrdnofmv.push_back(ANNIE_rates->GetBinContent(6));
      rate_fmv.push_back(ANNIE_rates->GetBinContent(7));
      rate_mrd.push_back(ANNIE_rates->GetBinContent(8));
      
      frac.push_back(ANNIE_fractions->GetBinContent(1));
      frac_pmt.push_back(ANNIE_fractions->GetBinContent(2));
      frac_pmtmrd.push_back(ANNIE_fractions->GetBinContent(5));
      frac_pmtmrdfmv.push_back(ANNIE_fractions->GetBinContent(10));
      frac_pmtmrdnofmv.push_back(ANNIE_fractions->GetBinContent(6));
      frac_fmv.push_back(ANNIE_fractions->GetBinContent(7));
      frac_mrd.push_back(ANNIE_fractions->GetBinContent(8));
      run_numbers.push_back(i_run);
      fval->Close();
      delete fval;
    }
  }

  //Open output file and fill TGraphs
  std::stringstream filename_out;
  filename_out << "RunValidationStability_R" << current_run - number_runs << "-" << current_run << ".root";
  TFile *file_output = new TFile(filename_out.str().c_str(),"RECREATE");
  
  TGraph *gr_rate = new TGraph();
  TGraph *gr_rate_pmt = new TGraph();
  TGraph *gr_rate_pmtmrd = new TGraph(); 
  TGraph *gr_rate_pmtmrdfmv = new TGraph(); 
  TGraph *gr_rate_pmtmrdnofmv = new TGraph(); 
  TGraph *gr_rate_fmv = new TGraph();
  TGraph *gr_rate_mrd = new TGraph();
  TGraph *gr_frac = new TGraph();
  TGraph *gr_frac_pmt = new TGraph();
  TGraph *gr_frac_pmtmrd = new TGraph();
  TGraph *gr_frac_pmtmrdfmv = new TGraph();
  TGraph *gr_frac_pmtmrdnofmv = new TGraph();
  TGraph *gr_frac_fmv = new TGraph();
  TGraph *gr_frac_mrd = new TGraph();
  TGraph *gr_align_pmt_mrd = new TGraph();
  TGraph *gr_align_pmt_fmv = new TGraph();
  TGraph *gr_align_mrd_fmv = new TGraph();
  TGraph *gr_frac_prompt = new TGraph();
  TGraph *gr_frac_delayed = new TGraph();
  TGraph *gr_frac_beam = new TGraph();
  TGraph *gr_frac_led = new TGraph();
  TGraph *gr_frac_cosmic = new TGraph();
  TGraph *gr_avg_mult = new TGraph();
  TGraph *gr_avg_mult_coinc = new TGraph();
  TGraph *gr_avg_mult_coinc_nofmv_cb = new TGraph();

  gr_rate->GetXaxis()->SetTitle("Run Number");
  gr_rate_pmt->GetXaxis()->SetTitle("Run Number");
  gr_rate_pmtmrd->GetXaxis()->SetTitle("Run Number");
  gr_rate_pmtmrdfmv->GetXaxis()->SetTitle("Run Number");
  gr_rate_pmtmrdnofmv->GetXaxis()->SetTitle("Run Number");
  gr_rate_fmv->GetXaxis()->SetTitle("Run Number");
  gr_rate_mrd->GetXaxis()->SetTitle("Run Number");
  gr_frac->GetXaxis()->SetTitle("Run Number");
  gr_frac_pmt->GetXaxis()->SetTitle("Run Number");
  gr_frac_pmtmrd->GetXaxis()->SetTitle("Run Number");
  gr_frac_pmtmrdfmv->GetXaxis()->SetTitle("Run Number");
  gr_frac_pmtmrdnofmv->GetXaxis()->SetTitle("Run Number");
  gr_frac_fmv->GetXaxis()->SetTitle("Run Number");
  gr_frac_mrd->GetXaxis()->SetTitle("Run Number");
  gr_align_pmt_mrd->GetXaxis()->SetTitle("Run Number");
  gr_align_pmt_fmv->GetXaxis()->SetTitle("Run Number");
  gr_align_mrd_fmv->GetXaxis()->SetTitle("Run Number");
  gr_frac_prompt->GetXaxis()->SetTitle("Run Number");
  gr_frac_delayed->GetXaxis()->SetTitle("Run Number");
  gr_frac_beam->GetXaxis()->SetTitle("Run Number");
  gr_frac_led->GetXaxis()->SetTitle("Run Number");
  gr_frac_cosmic->GetXaxis()->SetTitle("Run Number");
  gr_avg_mult->GetXaxis()->SetTitle("Run Number");
  gr_avg_mult_coinc->GetXaxis()->SetTitle("Run Number");
  gr_avg_mult_coinc_nofmv_cb->GetXaxis()->SetTitle("Run Number");

  gr_rate->GetYaxis()->SetTitle("R [Hz]");
  gr_rate_pmt->GetYaxis()->SetTitle("R [Hz]");
  gr_rate_pmtmrd->GetYaxis()->SetTitle("R [Hz]");
  gr_rate_pmtmrdfmv->GetYaxis()->SetTitle("R [Hz]");
  gr_rate_pmtmrdnofmv->GetYaxis()->SetTitle("R [Hz]");
  gr_rate_fmv->GetYaxis()->SetTitle("R [Hz]");
  gr_rate_mrd->GetYaxis()->SetTitle("R [Hz]");
  gr_frac->GetYaxis()->SetTitle("Fraction");
  gr_frac_pmt->GetYaxis()->SetTitle("Fraction");
  gr_frac_pmtmrd->GetYaxis()->SetTitle("Fraction");
  gr_frac_pmtmrdfmv->GetYaxis()->SetTitle("Fraction");
  gr_frac_pmtmrdnofmv->GetYaxis()->SetTitle("Fraction");
  gr_frac_fmv->GetYaxis()->SetTitle("Fraction");
  gr_frac_mrd->GetYaxis()->SetTitle("Fraction");
  gr_align_pmt_mrd->GetYaxis()->SetTitle("#Delta t_{PMT-MRD}");
  gr_align_pmt_fmv->GetYaxis()->SetTitle("#Delta t_{PMT-FMV}");
  gr_align_mrd_fmv->GetYaxis()->SetTitle("#Delta t_{MRD-FMV}");
  gr_frac_prompt->GetYaxis()->SetTitle("Fraction");
  gr_frac_delayed->GetYaxis()->SetTitle("Fraction");
  gr_frac_beam->GetYaxis()->SetTitle("Fraction");
  gr_frac_led->GetYaxis()->SetTitle("Fraction");
  gr_frac_cosmic->GetYaxis()->SetTitle("Fraction");
  gr_avg_mult->GetYaxis()->SetTitle("Neutron multiplicity");
  gr_avg_mult_coinc->GetYaxis()->SetTitle("Neutron multiplicity");
  gr_avg_mult_coinc_nofmv_cb->GetYaxis()->SetTitle("Neutron multiplicity");
  
  gr_rate->SetLineStyle(2);
  gr_rate_pmt->SetLineStyle(2);
  gr_rate_pmtmrd->SetLineStyle(2);
  gr_rate_pmtmrdfmv->SetLineStyle(2);
  gr_rate_pmtmrdnofmv->SetLineStyle(2);
  gr_rate_fmv->SetLineStyle(2);
  gr_rate_mrd->SetLineStyle(2);
  gr_frac->SetLineStyle(2);
  gr_frac_pmt->SetLineStyle(2);
  gr_frac_pmtmrd->SetLineStyle(2);
  gr_frac_pmtmrdfmv->SetLineStyle(2);
  gr_frac_pmtmrdnofmv->SetLineStyle(2);
  gr_frac_fmv->SetLineStyle(2);
  gr_frac_mrd->SetLineStyle(2);
  gr_align_pmt_mrd->SetLineStyle(2);
  gr_align_pmt_fmv->SetLineStyle(2);
  gr_align_mrd_fmv->SetLineStyle(2);
  gr_frac_prompt->SetLineStyle(2);
  gr_frac_delayed->SetLineStyle(2);
  gr_frac_beam->SetLineStyle(2);
  gr_frac_led->SetLineStyle(2);
  gr_frac_cosmic->SetLineStyle(2);
  gr_avg_mult->SetLineStyle(2);
  gr_avg_mult_coinc->SetLineStyle(2);
  gr_avg_mult_coinc_nofmv_cb->SetLineStyle(2);
  
  gr_rate->SetLineWidth(2);
  gr_rate_pmt->SetLineWidth(2);
  gr_rate_pmtmrd->SetLineWidth(2);
  gr_rate_pmtmrdfmv->SetLineWidth(2);
  gr_rate_pmtmrdnofmv->SetLineWidth(2);
  gr_rate_fmv->SetLineWidth(2);
  gr_rate_mrd->SetLineWidth(2);
  gr_frac->SetLineWidth(2);
  gr_frac_pmt->SetLineWidth(2);
  gr_frac_pmtmrd->SetLineWidth(2);
  gr_frac_pmtmrdfmv->SetLineWidth(2);
  gr_frac_pmtmrdnofmv->SetLineWidth(2);
  gr_frac_fmv->SetLineWidth(2);
  gr_frac_mrd->SetLineWidth(2);
  gr_align_pmt_mrd->SetLineWidth(2);
  gr_align_pmt_fmv->SetLineWidth(2);
  gr_align_mrd_fmv->SetLineWidth(2);
  gr_frac_prompt->SetLineWidth(2);
  gr_frac_delayed->SetLineWidth(2);
  gr_frac_beam->SetLineWidth(2);
  gr_frac_led->SetLineWidth(2);
  gr_frac_cosmic->SetLineWidth(2);
  gr_avg_mult->SetLineWidth(2);
  gr_avg_mult_coinc->SetLineWidth(2);
  gr_avg_mult_coinc_nofmv_cb->SetLineWidth(2);
  
  gr_rate->SetMarkerStyle(21);
  gr_rate_pmt->SetMarkerStyle(21);
  gr_rate_pmtmrd->SetMarkerStyle(21);
  gr_rate_pmtmrdfmv->SetMarkerStyle(21);
  gr_rate_pmtmrdnofmv->SetMarkerStyle(21);
  gr_rate_fmv->SetMarkerStyle(21);
  gr_rate_mrd->SetMarkerStyle(21);
  gr_frac->SetMarkerStyle(21);
  gr_frac_pmt->SetMarkerStyle(21);
  gr_frac_pmtmrd->SetMarkerStyle(21);
  gr_frac_pmtmrdfmv->SetMarkerStyle(21);
  gr_frac_pmtmrdnofmv->SetMarkerStyle(21);
  gr_frac_fmv->SetMarkerStyle(21);
  gr_frac_mrd->SetMarkerStyle(21);
  gr_align_pmt_mrd->SetMarkerStyle(21);
  gr_align_pmt_fmv->SetMarkerStyle(21);
  gr_align_mrd_fmv->SetMarkerStyle(21);
  gr_frac_prompt->SetMarkerStyle(21);
  gr_frac_delayed->SetMarkerStyle(21);
  gr_frac_beam->SetMarkerStyle(21);
  gr_frac_led->SetMarkerStyle(21);
  gr_frac_cosmic->SetMarkerStyle(21);
  gr_avg_mult->SetMarkerStyle(21);
  gr_avg_mult_coinc->SetMarkerStyle(21);
  gr_avg_mult_coinc_nofmv_cb->SetMarkerStyle(21);
  
  gr_rate->SetMarkerSize(0.7);
  gr_rate_pmt->SetMarkerSize(0.7);
  gr_rate_pmtmrd->SetMarkerSize(0.7);
  gr_rate_pmtmrdfmv->SetMarkerSize(0.7);
  gr_rate_pmtmrdnofmv->SetMarkerSize(0.7);
  gr_rate_fmv->SetMarkerSize(0.7);
  gr_rate_mrd->SetMarkerSize(0.7);
  gr_frac->SetMarkerSize(0.7);
  gr_frac_pmt->SetMarkerSize(0.7);
  gr_frac_pmtmrd->SetMarkerSize(0.7);
  gr_frac_pmtmrdfmv->SetMarkerSize(0.7);
  gr_frac_pmtmrdnofmv->SetMarkerSize(0.7);
  gr_frac_fmv->SetMarkerSize(0.7);
  gr_frac_mrd->SetMarkerSize(0.7);
  gr_align_pmt_mrd->SetMarkerSize(0.7);
  gr_align_pmt_fmv->SetMarkerSize(0.7);
  gr_align_mrd_fmv->SetMarkerSize(0.7);
  gr_frac_prompt->SetMarkerSize(0.7);
  gr_frac_delayed->SetMarkerSize(0.7);
  gr_frac_beam->SetMarkerSize(0.7);
  gr_frac_led->SetMarkerSize(0.7);
  gr_frac_cosmic->SetMarkerSize(0.7);
  gr_avg_mult->SetMarkerSize(0.7);
  gr_avg_mult_coinc->SetMarkerSize(0.7);
  gr_avg_mult_coinc_nofmv_cb->SetMarkerSize(0.7);

  std::stringstream ss_title_pmt, ss_title_pmt_2pe, ss_title_pmt_cosmics, ss_title_pmt_2pe_cosmics, ss_title_pmt_led, ss_title_pmt_2pe_led, ss_title_mrd;
  std::stringstream ss_title_mult, ss_title_mult_coinc, ss_title_mult_coinc_nofmv_cb, ss_title_trigword, ss_title_adc, ss_title_delayedt, ss_title_delayedt_coinc, ss_title_delayedt_coinc_nofmv, ss_title_delayedt_coinc_nofmv_cb, ss_title_delayedt_cb, ss_title_mult_coinc_nofmv;
  ss_title_pmt << "PMT Cluster times (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_pmt_2pe << "PMT Cluster times (#geq 2p.e.,Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_pmt_cosmics << "PMT Cluster times Cosmics (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_pmt_2pe_cosmics << "PMT Cluster times Cosmics (#geq 2p.e., Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_pmt_led << "PMT Cluster times LED (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_pmt_2pe_led << "PMT Cluster times LED (#geq 2p.e., Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_mrd << "MRD Cluster times (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_mult << "Neutron multiplicity (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_mult_coinc << "Neutron multiplicity PMT/MRD Coincidence (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_mult_coinc_nofmv << "Neutron multiplicity PMT/MRD Coincidence, No FMV, CB < 0.4 (Run "<<current_run - number_runs << "- Run "<<current_run << ")";
  ss_title_trigword << "Triggerwords (Run "<<current_run - number_runs << "- Run "<<current_run <<")";
  ss_title_adc << "ADC Waveform Samples (Run "<<current_run - number_runs <<"- Run "<<current_run <<")";
  ss_title_delayedt << "Neutron candidate time distribution (Run "<<current_run - number_runs << "- Run "<<current_run <<")";
  ss_title_delayedt_cb << "Neutron candidate time distribution CB < 0.4 (Run "<<current_run - number_runs << "- Run "<<current_run <<")";
  ss_title_delayedt_coinc << "Neutron candidate time distribution PMT/MRD Coincidence (Run "<<current_run - number_runs << "- Run "<<current_run <<")";
  ss_title_delayedt_coinc_nofmv << "Neutron candidate time distribution PMT/MRD Coinc, No FMV (Run "<<current_run - number_runs << "- Run "<<current_run <<")";
  ss_title_delayedt_coinc_nofmv_cb << "Neutron candidate time distribution PMT/MRD Coinc, No FMV, CB < 0.4 (Run "<<current_run - number_runs << "- Run "<<current_run <<")";
  PMT_t_clusters_combined->SetTitle(ss_title_pmt.str().c_str());
  PMT_t_clusters_2pe_combined->SetTitle(ss_title_pmt_2pe.str().c_str());
  PMT_t_clusters_cosmics_combined->SetTitle(ss_title_pmt_cosmics.str().c_str());
  PMT_t_clusters_2pe_cosmics_combined->SetTitle(ss_title_pmt_2pe_cosmics.str().c_str());
  PMT_t_clusters_led_combined->SetTitle(ss_title_pmt_led.str().c_str());
  PMT_t_clusters_2pe_led_combined->SetTitle(ss_title_pmt_2pe_led.str().c_str());
  MRD_t_clusters_combined->SetTitle(ss_title_mrd.str().c_str());
  PMT_DelayedMult_combined->SetTitle(ss_title_mult.str().c_str());
  PMT_DelayedMult_Coinc_combined->SetTitle(ss_title_mult_coinc.str().c_str());
  PMT_DelayedMult_Coinc_NoFMV_CB_combined->SetTitle(ss_title_mult_coinc_nofmv_cb.str().c_str());
  Triggerwords_combined->SetTitle(ss_title_trigword.str().c_str());
  ADCWaveform_Samples_combined->Write();
  PMT_DelayedTime_combined->SetTitle(ss_title_delayedt.str().c_str());
  PMT_DelayedTime_CB_combined->SetTitle(ss_title_delayedt_cb.str().c_str());
  PMT_DelayedTime_Coinc_combined->SetTitle(ss_title_delayedt_coinc.str().c_str());
  PMT_DelayedTime_Coinc_NoFMV_combined->SetTitle(ss_title_delayedt_coinc_nofmv.str().c_str());
  PMT_DelayedTime_Coinc_NoFMV_CB_combined->SetTitle(ss_title_delayedt_coinc_nofmv_cb.str().c_str());
  
  for (int i_run=0; i_run < (int) run_numbers.size(); i_run++){

    gr_rate->SetPoint(i_run,run_numbers.at(i_run),rate.at(i_run));
    gr_rate_pmt->SetPoint(i_run,run_numbers.at(i_run),rate_pmt.at(i_run));
    gr_rate_pmtmrd->SetPoint(i_run,run_numbers.at(i_run),rate_pmtmrd.at(i_run));
    gr_rate_pmtmrdfmv->SetPoint(i_run,run_numbers.at(i_run),rate_pmtmrdfmv.at(i_run));
    gr_rate_pmtmrdnofmv->SetPoint(i_run,run_numbers.at(i_run),rate_pmtmrdnofmv.at(i_run));
    gr_rate_fmv->SetPoint(i_run,run_numbers.at(i_run),rate_fmv.at(i_run));
    gr_rate_mrd->SetPoint(i_run,run_numbers.at(i_run),rate_mrd.at(i_run));
    gr_frac->SetPoint(i_run,run_numbers.at(i_run),frac.at(i_run));
    gr_frac_pmt->SetPoint(i_run,run_numbers.at(i_run),frac_pmt.at(i_run));
    gr_frac_pmtmrd->SetPoint(i_run,run_numbers.at(i_run),frac_pmtmrd.at(i_run));
    gr_frac_pmtmrdfmv->SetPoint(i_run,run_numbers.at(i_run),frac_pmtmrdfmv.at(i_run));
    gr_frac_pmtmrdnofmv->SetPoint(i_run,run_numbers.at(i_run),frac_pmtmrdnofmv.at(i_run));
    gr_frac_fmv->SetPoint(i_run,run_numbers.at(i_run),frac_fmv.at(i_run));
    gr_frac_mrd->SetPoint(i_run,run_numbers.at(i_run),frac_mrd.at(i_run));
    gr_align_pmt_mrd->SetPoint(i_run,run_numbers.at(i_run),align_pmt_mrd.at(i_run));
    gr_align_pmt_fmv->SetPoint(i_run,run_numbers.at(i_run),align_pmt_fmv.at(i_run));
    gr_align_mrd_fmv->SetPoint(i_run,run_numbers.at(i_run),align_mrd_fmv.at(i_run));
    gr_frac_prompt->SetPoint(i_run,run_numbers.at(i_run),frac_prompt.at(i_run));
    gr_frac_delayed->SetPoint(i_run,run_numbers.at(i_run),frac_delayed.at(i_run));
    gr_frac_beam->SetPoint(i_run,run_numbers.at(i_run),frac_beam.at(i_run));
    gr_frac_led->SetPoint(i_run,run_numbers.at(i_run),frac_led.at(i_run));
    gr_frac_cosmic->SetPoint(i_run,run_numbers.at(i_run),frac_cosmic.at(i_run));
    gr_avg_mult->SetPoint(i_run,run_numbers.at(i_run),avg_mult.at(i_run));
    gr_avg_mult_coinc->SetPoint(i_run,run_numbers.at(i_run),avg_mult_coinc.at(i_run));
    gr_avg_mult_coinc_nofmv_cb->SetPoint(i_run,run_numbers.at(i_run),avg_mult_coinc_nofmv_cb.at(i_run));
  }

  TMultiGraph *multi_rate = new TMultiGraph();
  TMultiGraph *multi_frac = new TMultiGraph();
  TMultiGraph *multi_types = new TMultiGraph();
  TMultiGraph *multi_mult = new TMultiGraph();
  TLegend *leg_rate = new TLegend(0.7,0.7,0.88,0.88);
  TLegend *leg_frac = new TLegend(0.7,0.7,0.88,0.88);
  TLegend *leg_types = new TLegend(0.7,0.7,0.88,0.88);
  TLegend *leg_mult = new TLegend(0.7,0.7,0.88,0.88);
  TCanvas *canvas_rate = new TCanvas("canvas_rate","Rate Multigraph",900,600);
  TCanvas *canvas_frac = new TCanvas("canvas_frac","Fractions Multigraph",900,600);
  TCanvas *canvas_types = new TCanvas("canvas_types","Trigger Types Multigraph",900,600);
  TCanvas *canvas_mult = new TCanvas("canvas_mult","Multiplicity Multigraph",900,600);

  gr_rate->SetLineColor(0);
  gr_rate->SetMarkerColor(0);
  multi_rate->Add(gr_rate);
  leg_rate->AddEntry(gr_rate,"Event rate","l");
  gr_rate_pmt->SetLineColor(kOrange);
  gr_rate_pmt->SetMarkerColor(kOrange);
  multi_rate->Add(gr_rate_pmt);
  leg_rate->AddEntry(gr_rate_pmt,"PMT Cluster rate","l");
  gr_rate_pmtmrd->SetLineColor(8);
  gr_rate_pmtmrd->SetMarkerColor(8);
  multi_rate->Add(gr_rate_pmtmrd);
  leg_rate->AddEntry(gr_rate_pmtmrd,"PMT+MRD Cluster rate","l");
  gr_rate_pmtmrdfmv->SetLineColor(kRed);
  gr_rate_pmtmrdfmv->SetMarkerColor(kRed);
  multi_rate->Add(gr_rate_pmtmrdfmv);
  leg_rate->AddEntry(gr_rate_pmtmrdfmv,"PMT+MRD+FMV Cluster rate","l");
  gr_rate_pmtmrdnofmv->SetLineColor(kBlue);
  gr_rate_pmtmrdnofmv->SetMarkerColor(kBlue);
  multi_rate->Add(gr_rate_pmtmrdnofmv);
  leg_rate->AddEntry(gr_rate_pmtmrdnofmv,"PMT_MRD (No FMV) Cluster rate","l");
  canvas_rate->cd();
  multi_rate->Draw("apl");
  multi_rate->GetYaxis()->SetTitle("Rate [Hz]");
  multi_rate->GetXaxis()->SetTitle("Run number");
  leg_rate->Draw();

  gr_frac_pmt->SetLineColor(kOrange);
  gr_frac_pmt->SetMarkerColor(kOrange);
  multi_frac->Add(gr_frac);
  leg_frac->AddEntry(gr_frac,"Event fraction","l");
  multi_frac->Add(gr_frac_pmt);
  leg_frac->AddEntry(gr_frac_pmt,"PMT Cluster fraction","l");
  gr_frac_pmtmrd->SetLineColor(8);
  gr_frac_pmtmrd->SetMarkerColor(8);
  multi_frac->Add(gr_frac_pmtmrd);
  leg_frac->AddEntry(gr_frac_pmtmrd,"PMT+MRD Cluster fraction","l");
  gr_frac_pmtmrdfmv->SetLineColor(kRed);
  gr_frac_pmtmrdfmv->SetMarkerColor(kRed);
  multi_frac->Add(gr_frac_pmtmrdfmv);
  leg_frac->AddEntry(gr_frac_pmtmrdfmv,"PMT+MRD+FMV Cluster fraction","l");
  gr_frac_pmtmrdnofmv->SetLineColor(kBlue);
  gr_frac_pmtmrdnofmv->SetMarkerColor(kBlue);
  multi_frac->Add(gr_frac_pmtmrdnofmv);
  leg_frac->AddEntry(gr_frac_pmtmrdnofmv,"PMT_MRD (No FMV) Cluster fraction","l");
  canvas_frac->cd();
  multi_frac->Draw("apl");
  multi_frac->GetYaxis()->SetTitle("Fraction");
  multi_frac->GetXaxis()->SetTitle("Run Number");
  leg_frac->Draw();

  gr_frac_prompt->SetLineColor(kBlue);
  gr_frac_prompt->SetMarkerColor(kBlue);
  gr_frac_delayed->SetLineColor(kViolet);
  gr_frac_delayed->SetMarkerColor(kViolet);
  gr_frac_beam->SetLineColor(kOrange);
  gr_frac_beam->SetMarkerColor(kOrange);
  gr_frac_led->SetLineColor(kRed);
  gr_frac_led->SetMarkerColor(kRed);
  gr_frac_cosmic->SetLineColor(8);
  gr_frac_cosmic->SetMarkerColor(8);
  multi_types->Add(gr_frac_prompt);
  leg_types->AddEntry(gr_frac_prompt,"prompt window","l");
  multi_types->Add(gr_frac_delayed);
  leg_types->AddEntry(gr_frac_delayed,"delayed window","l");
  multi_types->Add(gr_frac_beam);
  leg_types->AddEntry(gr_frac_beam,"Beam Trigger","l");
  multi_types->Add(gr_frac_led);
  leg_types->AddEntry(gr_frac_led,"LED Trigger","l");
  multi_types->Add(gr_frac_cosmic);
  leg_types->AddEntry(gr_frac_cosmic,"Cosmic Trigger","l");
  canvas_types->cd();
  multi_types->Draw("apl");
  multi_types->GetYaxis()->SetTitle("Fraction");
  multi_types->GetXaxis()->SetTitle("Run Number");
  leg_types->Draw();

  gr_avg_mult->SetLineColor(0);
  gr_avg_mult->SetMarkerColor(0);
  gr_avg_mult_coinc->SetLineColor(kBlue);
  gr_avg_mult_coinc->SetMarkerColor(kBlue);
  gr_avg_mult_coinc_nofmv_cb->SetLineColor(8);
  gr_avg_mult_coinc_nofmv_cb->SetMarkerColor(8);
  multi_mult->Add(gr_avg_mult);
  leg_mult->AddEntry(gr_avg_mult,"Average neutron multiplicity","l");
  multi_mult->Add(gr_avg_mult_coinc);
  leg_mult->AddEntry(gr_avg_mult_coinc,"Average neutron multiplicity (MRD+Tank coincidence)","l");
  multi_mult->Add(gr_avg_mult_coinc);
  leg_mult->AddEntry(gr_avg_mult_coinc_nofmv_cb,"Average neutron multiplicity (MRD+Tank, no FMV)","l");
  canvas_mult->cd();
  multi_mult->Draw("apl");
  multi_mult->GetYaxis()->SetTitle("Multiplicity");
  multi_mult->GetXaxis()->SetTitle("Run Number");
  leg_mult->Draw();

  file_output->cd();
  gr_rate->Write("gr_rate");
  gr_rate_pmt->Write("gr_rate_pmt");
  gr_rate_pmtmrd->Write("gr_rate_pmtmrd");
  gr_rate_pmtmrdfmv->Write("gr_rate_pmtmrdfmv");
  gr_rate_pmtmrdnofmv->Write("gr_rate_pmtmrdnofmv");
  gr_rate_fmv->Write("gr_rate_fmv");
  gr_rate_mrd->Write("gr_rate_mrd");
  gr_frac->Write("gr_frac");
  gr_frac_pmt->Write("gr_frac_pmt");
  gr_frac_pmtmrd->Write("gr_frac_pmtmrd");
  gr_frac_pmtmrdfmv->Write("gr_frac_pmtmrdfmv");
  gr_frac_pmtmrdnofmv->Write("gr_frac_pmtmrdnofmv");
  gr_frac_fmv->Write("gr_frac_fmv");
  gr_align_pmt_mrd->Write("gr_align_pmt_mrd");
  gr_align_pmt_fmv->Write("gr_align_pmt_fmv");
  gr_align_mrd_fmv->Write("gr_align_mrd_fmv");
  gr_frac_prompt->Write("gr_frac_prompt");
  gr_frac_delayed->Write("gr_frac_delayed");
  gr_frac_beam->Write("gr_frac_beam");
  gr_frac_led->Write("gr_frac_led");
  gr_frac_cosmic->Write("gr_frac_cosmic");
  gr_avg_mult->Write("gr_avg_mult");
  gr_avg_mult_coinc->Write("gr_avg_mult_coinc");
  gr_avg_mult_coinc_nofmv_cb->Write("gr_avg_mult_coinc_nofmv_cb");
  PMT_t_clusters_combined->Write();
  PMT_t_clusters_2pe_combined->Write();
  PMT_t_clusters_cosmics_combined->Write();
  PMT_t_clusters_2pe_cosmics_combined->Write();
  PMT_t_clusters_led_combined->Write();
  PMT_t_clusters_2pe_led_combined->Write();
  MRD_t_clusters_combined->Write();
  PMT_DelayedMult_combined->Write();
  PMT_DelayedMult_Coinc_combined->Write();
  PMT_DelayedMult_Coinc_NoFMV_CB_combined->Write();
  Triggerwords_combined->Write();
  ADCWaveform_Samples_combined->Write();
  PMT_DelayedTime_combined->Write();
  PMT_DelayedTime_CB_combined->Write();
  PMT_DelayedTime_Coinc_combined->Write();
  PMT_DelayedTime_Coinc_NoFMV_combined->Write();
  PMT_DelayedTime_Coinc_NoFMV_CB_combined->Write();
  canvas_rate->Write();
  canvas_frac->Write();
  canvas_types->Write();
  canvas_mult->Write();
  file_output->Close();
  delete file_output;

  return true;
}


bool RunValidationStability::Finalise(){

  return true;
}
