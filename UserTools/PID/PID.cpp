#include "PID.h"
#include "HighEReco.h"

PID::PID():Tool(){}

    bool PID::Initialise(std::string configfile, DataModel &data){

    std::cout <<"PID: Initialising..."<<std::endl;

    //-------------------------------------------------------------------------------
    //--------------------------- Useful header -------------------------------------
    //-------------------------------------------------------------------------------

    if(configfile!="")  m_variables.Initialise(configfile);
    //m_variables.Print();
    m_data= &data;

    m_variables.Get("verbosity",verbosity);
    m_variables.Get("OutputFile",out_filename);
    m_variables.Get("OutputFileText",out_filename_text);
    m_variables.Get("PDFFile_electron",pdf_filename_e);
    m_variables.Get("PDFFile_muon",pdf_filename_mu);
    m_variables.Get("UseParametric",parametric);
    m_variables.Get("IsSimulation",is_simulation);
    m_variables.Get("DrawDebug",draw_debug);

    m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry",geom);

    int get_ok = m_data->Stores.at("ANNIEEvent")->Header->Get("AnnieGeometry",geom);
    if(get_ok==0) cerr<<"couldn't find the AnnieGeometry! :("<<endl;
    n_tank_pmts = geom->GetNumDetectorsInSet("Tank");
    n_lappds = geom->GetNumDetectorsInSet("LAPPD");
    n_mrd_pmts = geom->GetNumDetectorsInSet("MRD");
    n_veto_pmts = geom->GetNumDetectorsInSet("Veto");

    m_data->CStore.Get("detectorkey_to_lappdid",detectorkey_to_lappdid);
    m_data->CStore.Get("channelkey_to_pmtid",channelkey_to_pmtid);
    m_data->CStore.Get("channelkey_to_mrdpmtid",channelkey_to_mrdpmtid);
    m_data->CStore.Get("channelkey_to_faccpmtid",channelkey_to_faccpmtid);
    m_data->CStore.Get("pmt_tubeid_to_channelkey",pmt_tubeid_to_channelkey);

    //------------------------------------------------------------------------------
    //-------SETUP all variables----------------------------------------------------
    //------------------------------------------------------------------------------

    Position detector_center = geom->GetTankCentre();
    tank_center_x = detector_center.X();
    tank_center_y = detector_center.Y();
    tank_center_z = detector_center.Z();
    std::cout <<"PID: Num Tank PMTs: "<<n_tank_pmts<<", num MRD PMTs: "<<n_mrd_pmts<<", num Veto PMTs: "<<n_veto_pmts<<", num LAPPDs: "<<n_lappds<<std::endl;
    std::vector<std::map<unsigned long,Detector>* >* Detectors = geom->GetDetectors();
    std::cout <<"Detectors size: "<<Detectors->size()<<std::endl;

    for (int i_detector = 0; i_detector < Detectors->size(); i_detector++){

  	std::map<unsigned long, Detector>* temp_map = Detectors->at(i_detector);

  	for(auto&& apmt : (*temp_map)){
  		unsigned long detkey = apmt.first;
  		unsigned long chankey;
  		if (i_detector == 0) chankey = detkey*n_channels_per_lappd;
  		else if (i_detector == 1) chankey = pmt_tubeid_to_channelkey[detkey-n_lappds+1];
  		else chankey = n_channels_per_lappd*n_lappds+(detkey-n_lappds);
		Detector* thedetector = geom->ChannelToDetector(chankey);

		if (i_detector == 1) {
			//std::cout <<"PMT id: "<<channelkey_to_pmtid[chankey];
			Position position_PMT = thedetector->GetDetectorPosition();
			pmts_x.insert(std::pair<int,double>(channelkey_to_pmtid[chankey]-1,position_PMT.X()-tank_center_x));
			pmts_y.insert(std::pair<int,double>(channelkey_to_pmtid[chankey]-1,position_PMT.Y()-tank_center_y));
			pmts_z.insert(std::pair<int,double>(channelkey_to_pmtid[chankey]-1,position_PMT.Z()-tank_center_z));
		} 
  	}
    }

    //------------------------------------------------------------------------------
    ////-------Get e- and mu- PDFs----------------------------------------------------
    //------------------------------------------------------------------------------

    f_pdf_e = new TFile(pdf_filename_e.c_str(),"READ");
    electronPhotons = (TH3F*) f_pdf_e->Get("electrons_pe");
    electronTimes = (TH3F*) f_pdf_e->Get("electrons_time");
    f_pdf_mu = new TFile(pdf_filename_mu.c_str(),"READ");
    muonPhotons = (TH3F*) f_pdf_mu->Get("muons_pe");
    muonTimes = (TH3F*) f_pdf_mu->Get("muons_time");

    //------------------------------------------------------------------------------
    ////-------Define TH1s to save likelihood distributions---------------------------
    //------------------------------------------------------------------------------

    likelihood_electron = new TH1F("likelihood_electron","Likelihood distribution for electron events",1000,1,0);
    likelihood_electron->SetLineColor(4);
    likelihood_electron->GetXaxis()->SetTitle("Log Likelihood (e-) - Log Likelihood (mu-)");
    likelihood_muon = new TH1F("likelihood_muon","Likelihood distribution for muon events",1000,1,0);
    likelihood_muon->SetLineColor(2);
    likelihood_muon->GetXaxis()->SetTitle("Log Likelihood (e-) - Log Likelihood (mu-)");
    likelihood_electron_pdfE = new TH1F("likelihood_electron_pdfE","Likelihood distribution for electron events",1000,1,0);
    likelihood_electron_pdfE->SetLineColor(4);
    likelihood_electron_pdfE->GetXaxis()->SetTitle("Log Likelihood (e-)");
    likelihood_electron_pdfMu = new TH1F("likelihood_electron_pdfMu","Likelihood distribution for electron events",1000,1,0);
    likelihood_electron_pdfMu->SetLineColor(9);
    likelihood_electron_pdfMu->GetXaxis()->SetTitle("Log Likelihood (mu-)");
    likelihood_muon_pdfE = new TH1F("likelihood_muon_pdfE","Likelihood distribution for muon events",1000,1,0);
    likelihood_muon_pdfE->SetLineColor(2);
    likelihood_muon_pdfE->GetXaxis()->SetTitle("Log Likelihood (e-)");
    likelihood_muon_pdfMu = new TH1F("likelihood_muon_pdfMu","Likelihood distribution for muon events",1000,1,0);
    likelihood_muon_pdfMu->SetLineColor(9);
    likelihood_muon_pdfMu->GetXaxis()->SetTitle("Log Likelihood (mu-)");
    likelihood_electron_pe = new TH1F("likelihood_electron_pe","Likelihood distribution for electron events [charge]",1000,1,0);
    likelihood_electron_pe->SetLineColor(4);
    likelihood_electron_pe->GetXaxis()->SetTitle("Log Likelihood (e-) - Log Likelihood (mu-) [charge]");
    likelihood_muon_pe = new TH1F("likelihood_muon_pe","Likelihood distribution for muon events [charge]",1000,1,0);
    likelihood_muon_pe->SetLineColor(2);
    likelihood_muon_pe->GetXaxis()->SetTitle("Log Likelihood (e-) - Log Likelihood (mu-) [charge]");
    likelihood_electron_time = new TH1F("likelihood_electron_time","Likelihood distribution for electron events [time]",1000,1,0);
    likelihood_electron_time->SetLineColor(4);
    likelihood_electron_time->GetXaxis()->SetTitle("Log Likelihood (e-) - Log Likelihood (mu-) [time]");
    likelihood_muon_time = new TH1F("likelihood_muon_time","Likelihood distribution for muon events [time]",1000,1,0);
    likelihood_muon_time->SetLineColor(2);
    likelihood_muon_time->GetXaxis()->SetTitle("Log Likelihood (e-) - Log Likelihood (mu-) [time]");
    likelihood_nhits = new TH2F("likelihood_nhits","Likelihood (e - mu) : nhit PMTs",150,0,150,100000,-5000,5000);
    likelihood_nhits->GetXaxis()->SetTitle("nhit PMTs");
    likelihood_nhits->GetYaxis()->SetTitle("Likelihood (e-) - Likelihood (mu-)");
    

    file_out = new TFile(out_filename.c_str(),"RECREATE");

    //-------------------------------------------------------------------------------
    //---------define TApplication to debug while running----------------------------
    //-------------------------------------------------------------------------------

    //only if one needs to debug
    if (draw_debug){
      int myargc=0;
      char *myargv[] = {(const char*)"somestring"};
      app_pid = new TApplication("AppPID",&myargc,myargv);
    }

    loop = 0;
    time_cut_veto=0;
    n_electrons = 0;
    n_muons = 0;
    n_electrons_cuts = 0;
    n_muons_cuts = 0;
    n_reco_electrons = 0;
    n_reco_muons = 0;

    Log("PID tool: Initializing",v_message, verbosity);

    return true;

}


bool PID::Execute(){


  std::cout <<"PID tool: Executing..."<<std::endl;;

  get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCTriggernum",MCTriggernum);
  if(not get_ok){ Log("PID Tool: Error retrieving MCTriggernum from ANNIEEvent!",v_error,verbosity); return false; }
  if(MCTriggernum>0){ Log("PID Tool: Skipping delayed trigger",v_debug,verbosity); return true;}
  get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCFile",MCFile);
  if(not get_ok){ Log("PID Tool: Error retrieving MCFile from ANNIEEvent!",v_error,verbosity); return false; }
  get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCEventNum",MCEventNum);
  if(not get_ok){ Log("PID Tool: Error retrieving MCEventNum from ANNIEEvent!",v_error,verbosity); return false; }
  get_ok = m_data->Stores.at("ANNIEEvent")->Get("EventNumber",EventNumber);
  if(not get_ok){ Log("PID Tool: Error retrieving EventNumber from ANNIEEvent!",v_error,verbosity); return false; }
  get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCParticles",MCParticles);
  if(not get_ok){ Log("PID Tool: Error retrieving MCParticles,true from ANNIEEvent!",v_error,verbosity); return false; }
  get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCHits",MCHits);
  if(not get_ok){ Log("PID Tool: Error retrieving MCHits,true from ANNIEEvent!",v_error,verbosity); return false; }
  get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCLAPPDHits",MCLAPPDHits);
  if(not get_ok){ Log("PID Tool: Error retrieving MCLAPPDHits,true from ANNIEEvent!",v_error,verbosity); return false; }
  get_ok = m_data->Stores.at("ANNIEEvent")->Get("TDCData",TDCData);
  if(not get_ok){ Log("PID Tool: Error retrieving TDCData,true from ANNIEEvent!",v_error,verbosity); return false; }
  get_ok = m_data->Stores.at("ANNIEEvent")->Get("EventTime",EventTime);
  if(not get_ok){ Log("PID Tool: Error retrieving EventTime,true from ANNIEEvent!",v_error,verbosity); return false; }

  //----------------------------------------------------------------------------
  //-------------DEFINE histograms to show performance of application-------------- 
  //----------------------------------------------------------------------------

  std::stringstream ss_loop;
  ss_loop<<loop;
  /*std::string str_obs_pe = "obs_pe_";
  std::string str_exp_pe = "exp_pe_";
  std::string str_diff_pe = "diff_pe_";
  std::string str_obs_t = "obs_t_";
  std::string str_exp_t = "exp_t_";
  std::string str_diff_t = "diff_t_";
  std::string obs_pe_current = str_obs_pe+ss_loop.str();
  std::string exp_pe_current = str_exp_pe+ss_loop.str();
  std::string diff_pe_current = str_diff_pe+ss_loop.str();
  std::string obs_t_current = str_obs_t+ss_loop.str();
  std::string exp_t_current = str_exp_t+ss_loop.str();
  std::string diff_t_current = str_diff_t+ss_loop.str();

  obs_pe = new TH1F(obs_pe_current.c_str(),"Observed p.e.",100,1,0);
  exp_pe = new TH1F(exp_pe_current.c_str(),"Expected p.e.",100,1,0);
  diff_pe = new TH1F(diff_pe_current.c_str(),"Difference p.e.",100,1,0);
  obs_t = new TH1F(obs_t_current.c_str(),"Observed times",100,1,0);
  exp_t = new TH1F(exp_t_current.c_str(),"Expected times",100,1,0);
  diff_t = new TH1F(diff_t_current.c_str(),"Difference times",100,1,0);*/

  //----------------------------------------------------------------------------
  //--------LOAD reconstructed (true) properties from Reco Store (MCTruth)------
  //----------------------------------------------------------------------------

  MCParticle primary; 

  if (is_simulation){	           //load true variables in case of loading simulation files directly

    logmessage = "PID Tool: Processing truth tracks and digits for "+MCFile
                                +", MCEvent "+to_string(MCEventNum)+", MCTrigger "+to_string(MCTriggernum);
    Log(logmessage,v_debug,verbosity);
   	std::cout <<"Start of Event!"<<std::endl;
		std::cout <<"Trigger Number: "<<MCTriggernum<<", Event Number: "<<EventNumber<<std::endl;
    bool mufound=false;
    bool efound=false;

    if(MCParticles){
      Log("PID Tool: Num MCParticles = "+to_string(MCParticles->size()),v_message,verbosity);
      for(int particlei=0; particlei<MCParticles->size(); particlei++){
        MCParticle aparticle = MCParticles->at(particlei);
		    if(aparticle.GetParentPdg()!=0) continue; 
		    if (aparticle.GetFlag()!=0) continue;			//	flag only needed for non-neutrino sample to produce PDFs
		    if(aparticle.GetPdgCode()==13){	
		      primary = aparticle;
                      mufound=true;
		      is_electron=0;
		      break;
		    }
		    else if (aparticle.GetPdgCode()==11){
		      primary = aparticle;
		      efound=true;
		      is_electron=1;
		      break;
		    }
        else continue;
		  }
		}
		else {
			Log("PID Tool: No MCParticles in this event!",v_error,verbosity);
		}
		if (not mufound && not efound){
			Log("PID Tool: Neither electrons nor muons in this event!",v_error,verbosity);
			return true;
		}

	}
	else {

		//define here the readout behavior in the case of using a Reconstruction Store instead of MCTruth
		std::cout <<"RECO fit not implemented yet. :) [coming soon]"<<std::endl;

	}

	if (is_electron) n_electrons++;
	else n_muons++;

	//----------------------------------------------------------------------------
	//--------DEFINE HighEReco instance ------------------------------------------
	//----------------------------------------------------------------------------

	HighEReco HE_Event;

	if (is_simulation){ //read true values for vertex, direction and energy of primary

		vtx_x = primary.GetStartVertex().X()-tank_center_x;
		vtx_y = primary.GetStartVertex().Y()-tank_center_y;
		vtx_z = primary.GetStartVertex().Z()-tank_center_z;
		dir_x = primary.GetStartDirection().X();
		dir_y = primary.GetStartDirection().Y();
		dir_z = primary.GetStartDirection().Z();

		dir_theta = TMath::ACos(dir_z /TMath::Sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z));
    		dir_phi = TMath::ATan2(dir_y, dir_x);
		primary_energy = primary.GetStartEnergy();
		primary_pid = (is_electron)? 11 : 13;
		time_true = primary.GetStartTime();
		t_reco = time_true-EventTime->GetNs();
		std::cout <<"Particle ID: "<<primary_pid<<", true energy: "<<primary_energy<<", start time: "<<time_true<<", EventTime: "<<EventTime->GetNs()<<", event time: "<<t_reco<<std::endl;

		if (MCHits){

			std::vector<int> PMT_Hit(n_tank_pmts,0);

			//std::cout <<"MCHits present! Looping over entries..."<<std::endl;
			for (int i_pmt=0; i_pmt<n_tank_pmts;i_pmt++){
				pmt_charge[i_pmt]=0;
				pmt_time[i_pmt]=0;
			}
      			total_pe = 0;
			for(std::pair<unsigned long,std::vector<Hit>>&& apair : *MCHits){
        		unsigned long chankey = apair.first;
        		int wcsim_pmt_id = channelkey_to_pmtid[chankey];
        		//std::cout <<"chankey: "<<chankey<<", PMT id: "<<wcsim_pmt_id<<std::endl;
			Detector* thistube = geom->ChannelToDetector(chankey);
			if (thistube->GetDetectorElement()=="Tank"){
            			std::vector<Hit>& hits = apair.second;
            			double charge_pmt=0;
            			double t_pmt;

            		if (parametric){
              			//use parametric model from DigitBuilder to calculate smeared charge and time values
              			std::vector<double> hitTimes;
          			std::vector<double> hitCharges;
          			for(Hit& ahit : hits){
          				PMT_Hit.at(wcsim_pmt_id-1)=1;
                			hitTimes.push_back(ahit.GetTime()*1.0); 
            				hitCharges.push_back(ahit.GetCharge());
          			}
          			// Calculate median and sum
          			std::sort(hitTimes.begin(), hitTimes.end());
          			size_t timesize = hitTimes.size();
          			if (timesize % 2 == 0){
            				t_pmt = (hitTimes.at(timesize/2 - 1) + hitTimes.at(timesize/2))/2;
          			} else {
            				t_pmt = hitTimes.at(timesize/2);
          			}
          			for(std::vector<double>::iterator it = hitCharges.begin(); it != hitCharges.end(); ++it){
            	 			charge_pmt += *it;
          			}	

              			pmt_charge[wcsim_pmt_id-1] = charge_pmt;
              			pmt_time[wcsim_pmt_id-1] = t_pmt;
             			total_pe+=charge_pmt;
            		} else {  //only use first hit in case of non-parametric model

              		int loop_hits=0;
            		for (Hit& ahit : hits){
            			if (loop_hits==0){
                  			PMT_Hit.at(wcsim_pmt_id-1)=1;
                  			charge_pmt = ahit.GetCharge();
            		  		t_pmt = ahit.GetTime();
                		}
                		loop_hits++;
            		}

  			pmt_charge[wcsim_pmt_id-1] = charge_pmt;
  			pmt_time[wcsim_pmt_id-1] = t_pmt;
  			total_pe+=charge_pmt;
            		}
			}
		}
	} else {
		std::cout <<"No MCHits! Simulation mode is invalid!"<<std::endl;
	}

		HE_Event.electronPhotons = (TH3D*) electronPhotons;
		HE_Event.muonPhotons = (TH3D*) muonPhotons;
		HE_Event.electronTimes = (TH3D*) electronTimes;
		HE_Event.muonTimes = (TH3D*) muonTimes;

		double maxTime = 90.;		//implement later that all hits after 90ns will be ommited (fine tune exact number for the time cut)
		HE_Event.nPEs = total_pe;    //so far, no cut on time, but probably needed in future
		//find out which PMTs were hit
		HE_Event.nHitPMT = 0;
		for (int i_PMT=0;i_PMT<n_tank_pmts;i_PMT++){
			HE_Event.nHitPMT += (pmt_charge[i_PMT] > 0);
		}
		HE_Event.nUnhitPMT = n_tank_pmts - HE_Event.nHitPMT;
		std::cout <<"Hit PMTs: "<<HE_Event.nHitPMT<<", unHit PMTs: "<<HE_Event.nUnhitPMT<<std::endl;

		//----------------------------------------------------------------------------
		//--------Fill Hit + Unhit arrays for HighEReco object------------------------
		//----------------------------------------------------------------------------

		HE_Event.hitPMTx = new double[HE_Event.nHitPMT];
		HE_Event.hitPMTy = new double[HE_Event.nHitPMT];
		HE_Event.hitPMTz = new double[HE_Event.nHitPMT];
		HE_Event.hitPMTDirX = new double[HE_Event.nHitPMT];			//was used in TITUS to calculate detection efficiencies --> include also here?
		HE_Event.hitPMTDirY = new double[HE_Event.nHitPMT];
		HE_Event.hitPMTDirZ = new double[HE_Event.nHitPMT];
		HE_Event.hitPMTTimeRes = new double[HE_Event.nHitPMT];
		HE_Event.hitPMTPEs = new int[HE_Event.nHitPMT];
		HE_Event.hitT = new double[HE_Event.nHitPMT];
		HE_Event.unhitPMTx = new double[HE_Event.nUnhitPMT];
		HE_Event.unhitPMTy = new double[HE_Event.nUnhitPMT];
		HE_Event.unhitPMTz = new double[HE_Event.nUnhitPMT];
		HE_Event.unhitPMTDirX = new double[HE_Event.nUnhitPMT];
		HE_Event.unhitPMTDirY = new double[HE_Event.nUnhitPMT];
		HE_Event.unhitPMTDirZ = new double[HE_Event.nUnhitPMT];
		HE_Event.unhitPMTTimeRes = new double[HE_Event.nUnhitPMT];

		int i_hit=0;
		int i_unhit =0;

		for (int i_pmt = 0;i_pmt < n_tank_pmts; i_pmt++){
			int det_id = i_pmt;
			double x_pmt = pmts_x.at(det_id);
			double y_pmt = pmts_y.at(det_id);
			double z_pmt = pmts_z.at(det_id);
			double x_dir_pmt = 1.;  //orientation of PMT (x)
			double y_dir_pmt = 0.;  //orientation of PMT (y)
			double z_dir_pmt = 0.;  //orientation of PMT (z)
			if (pmt_charge[i_pmt] > 0) {
				HE_Event.hitPMTx[i_hit] = x_pmt;
				HE_Event.hitPMTy[i_hit] = y_pmt;
				HE_Event.hitPMTz[i_hit] = z_pmt;
				HE_Event.hitPMTDirX[i_hit] = x_dir_pmt;
				HE_Event.hitPMTDirY[i_hit] = y_dir_pmt;
				HE_Event.hitPMTDirZ[i_hit] = z_dir_pmt;
				HE_Event.hitPMTTimeRes[i_hit] = 1.0;	//assume 1 ns for all PMTs
				HE_Event.hitT[i_hit] = pmt_time[i_pmt];
				HE_Event.hitPMTPEs[i_hit] = int(pmt_charge[i_pmt]);
				i_hit++;
			}else {
				HE_Event.unhitPMTx[i_unhit] = x_pmt;
				HE_Event.unhitPMTy[i_unhit] = y_pmt;
				HE_Event.unhitPMTz[i_unhit] = z_pmt;
				HE_Event.unhitPMTDirX[i_unhit] = x_dir_pmt;
				HE_Event.unhitPMTDirY[i_unhit] = y_dir_pmt;
				HE_Event.unhitPMTDirZ[i_unhit] = z_dir_pmt;
				HE_Event.unhitPMTTimeRes[i_unhit] = 1.0;	//assume 1 ns for all PMTs
				i_unhit++;
			}
		}

		std::cout <<"Check: Unhit PMTs: "<<i_unhit<<", hit PMTs: "<<i_hit<<std::endl;

	} else {//read reco variables for vertex, direction, energy

		std::cout <<"RECO fit not implemented yet. :) [coming soon]"<<std::endl;

	}

	//----------------------------------------------------------------------------
	//--------PERFORM Likelihood fit----------------------------------------------
	//----------------------------------------------------------------------------

	//no distinction between simulation and reco store anymore, the read in variables are now treated exactly the same
	//structure of the HighEReco Likelihood fit in the following lines
	 /*void LikelihoodFit(double &trackCorrection, double &recoVtxX, double &recoVtxY,
	double &recoVtxZ, double &recoT, double &recoDirPhi,
	double &recoDirTheta, double &recoKE, double &recoLnL, int ipnu);*/


    double lambda_min = 10000000;  double lambda_max = -99999999.9; double vertical_dist_max = -99999999.9; 

    if ((t_reco < 100. && t_reco > -100.) && HE_Event.nHitPMT != 0){					//how many events are removed through that?

    	if (is_electron) n_electrons_cuts++;
    	else n_muons_cuts++;

    	double track_correction = 0.;
    	std::cout << "------------------------------- LIKELIHOOD FIT -------------------------------------" << std::endl;
    	std::cout <<"Start parameters for Likelihood fit: track correction: "<<track_correction<<", vtx = ("<<vtx_x<<" , "<<vtx_y<<" , "<<vtx_z<<"), time: "<<t_reco<<", dir phi: "<<dir_phi<<", dir theta: "<<dir_theta<<", energy: "<<primary_energy<<std::endl;
    
    	//electron hypothesis fit
    	std::cout <<"---------------------------- Electron hypothesis fit --------------------------------"<<std::endl;
    	HE_Event.LikelihoodFit(track_correction,vtx_x,vtx_y,vtx_z,t_reco,dir_phi, dir_theta,primary_energy,likelihood_value_electron,11);
    	HE_Event.LikelihoodFitCharge(track_correction,vtx_x,vtx_y,vtx_z,t_reco,dir_phi, dir_theta,primary_energy,likelihood_value_electron_pe,11);
    	HE_Event.LikelihoodFitTime(track_correction,vtx_x,vtx_y,vtx_z,t_reco,dir_phi, dir_theta,primary_energy,likelihood_value_electron_time,11);
  	//muon hypothesis fit
  	std::cout <<"------------------------------- Muon hypothesis fit ---------------------------------"<<std::endl;
  	HE_Event.LikelihoodFit(track_correction,vtx_x,vtx_y,vtx_z,t_reco,dir_phi, dir_theta, primary_energy,likelihood_value_muon,13);
    	HE_Event.LikelihoodFitCharge(track_correction,vtx_x,vtx_y,vtx_z,t_reco,dir_phi, dir_theta, primary_energy,likelihood_value_muon_pe,13);
   	HE_Event.LikelihoodFitTime(track_correction,vtx_x,vtx_y,vtx_z,t_reco,dir_phi, dir_theta, primary_energy,likelihood_value_muon_time,13);

	likelihood_nhits->Fill(HE_Event.nHitPMT,likelihood_value_electron - likelihood_value_muon);

  	std::string print_pid;
  	print_pid = (is_electron)? "electron" : "muon";
  	std::cout << "Likelihood Fit: "<<print_pid<<", likelihood e-: "<<likelihood_value_electron<<", likelihood mu-: "<<likelihood_value_muon<<", likelihood: "<<likelihood_value_electron - likelihood_value_muon<<std::endl;
	std::cout << "------------------------------------------------------------------------------------"<<std::endl;
  	if (is_electron) {
  		likelihood_electron->Fill(likelihood_value_electron-likelihood_value_muon);
  		likelihood_electron_pdfE->Fill(likelihood_value_electron);
  		likelihood_electron_pdfMu->Fill(likelihood_value_muon);
      		likelihood_electron_pe->Fill(likelihood_value_electron_pe-likelihood_value_muon_pe);
      		likelihood_electron_time->Fill(likelihood_value_electron_time-likelihood_value_muon_time);
  		if (likelihood_value_electron-likelihood_value_muon < 0) n_reco_electrons++;
  	}
  	else {
  		likelihood_muon->Fill(likelihood_value_electron-likelihood_value_muon); 
  		likelihood_muon_pdfE->Fill(likelihood_value_electron);
  		likelihood_muon_pdfMu->Fill(likelihood_value_muon);
      		likelihood_muon_pe->Fill(likelihood_value_electron_pe-likelihood_value_muon_pe);
      		likelihood_muon_time->Fill(likelihood_value_electron_time-likelihood_value_muon_time);
  		if (likelihood_value_electron-likelihood_value_muon > 0) n_reco_muons++;
  	}

 	 } else {

		time_cut_veto++;
		std::cout <<"-------------------------NO LIKELIHOOD FIT-------------------------------------------"<<std::endl;
		std::cout <<"Time parameter out of range: t_event = "<<t_reco<<std::endl;
		std::cout << "------------------------------------------------------------------------------------"<<std::endl;

    		UnUsedEvNum.push_back(EventNumber);
    		UnUsedMCEvNum.push_back(MCEventNum);
    		UnUsedRho.push_back(sqrt(vtx_x*vtx_x+vtx_z*vtx_z));
    		UnUsedY.push_back(vtx_y);
    		UnUsedT.push_back(t_reco);
    		UnUsedNHits.push_back(HE_Event.nHitPMT);
    		UnUsedEnergy.push_back(primary_energy);
	}

	loop++;

	/* delete obs_pe;
	delete exp_pe;
	delete diff_pe;
	delete obs_t;
	delete exp_t;
	delete diff_t; */

	return true;

}

bool PID::Finalise(){

	//----------------------------------------------------------------------------
	//-------------WRITE histos to output root file ------------------------------
	//----------------------------------------------------------------------------

	file_out->cd();
	likelihood_electron->Write();
	likelihood_electron_pdfE->Write();
	likelihood_electron_pdfMu->Write();
	likelihood_muon->Write();
	likelihood_muon_pdfE->Write();
	likelihood_muon_pdfMu->Write();
  	likelihood_electron_pe->Write();
  	likelihood_electron_time->Write();
  	likelihood_muon_pe->Write();
  	likelihood_muon_time->Write();
	likelihood_nhits->Write();
	file_out->Close();
	std::cout <<"Number of times for time cut veto: "<<time_cut_veto<<std::endl;

	std::cout <<"-------------------------EFFICIENCIES---------------------------"<<std::endl;
	std::cout <<"Number of tagged electrons (after cuts): "<< double(n_reco_electrons)/n_electrons_cuts<<std::endl;
	std::cout <<"Number of tagged muons (after cuts): "<< double(n_reco_muons)/n_muons_cuts<<std::endl;
	std::cout <<"Number of tagged electrons (total): "<<  double(n_reco_electrons)/n_electrons<<std::endl;
	std::cout <<"Number of tagged muons (total): "<<  double(n_reco_muons)/n_muons<<std::endl;
	std::cout <<"Numbers: "<<std::endl;
	std::cout <<"Tagged electrons: "<<n_reco_electrons<<", tagged muons: "<<n_reco_muons<<std::endl;
	std::cout <<"# electrons (total): "<<n_electrons<<", # muons (total): "<<n_muons<<std::endl;
	std::cout <<"# electrons (cuts): "<<n_electrons_cuts<<", # muons (cuts): "<<n_muons_cuts<<std::endl;
  	std::cout <<"-----------------------------------------------------------------"<<std::endl;

  	ofstream outfile_efficiencies(out_filename_text.c_str());
  	outfile_efficiencies <<"-- n_e --"<<"  "<<"--n_e_cuts--"<<"  "<<"--n_e_tagged--"<<"  "<<"--n_mu--"<<"  "<<"--n_mu_cuts--"<<"  "<<"--n_mu_tagged--"<<std::endl;
  	outfile_efficiencies << n_electrons << "  " << n_electrons_cuts << "  "<< n_reco_electrons<<"  "<<n_muons<<"  "<<n_muons_cuts<<"  "<<n_reco_muons<<std::endl;
  	outfile_efficiencies.close();

  	ofstream outfile_skipped("skipped_events.txt");
  	outfile_skipped << "EvNum -- EvNum (MC) -- Rho -- Y -- T -- NHits -- Energy --"<<std::endl;
  	for (int i_skip=0;i_skip<UnUsedEnergy.size();i_skip++){
    		outfile_skipped << UnUsedEvNum.at(i_skip) <<"  "<<UnUsedMCEvNum.at(i_skip)<<"  "<<UnUsedRho.at(i_skip)<<"  "<<UnUsedY.at(i_skip)<<"  "<<UnUsedT.at(i_skip)<<"  "<<UnUsedNHits.at(i_skip)<<"  "<<UnUsedEnergy.at(i_skip)<<std::endl;
  	}
 	outfile_skipped.close();

	//----------------------------------------------------------------------------
	//---------DELETE all variables (close files, etc,pp...) ---------------------
	//----------------------------------------------------------------------------

	f_pdf_e->Close();
	f_pdf_mu->Close();

	return true;

}

