#include "ParticleIDPDF.h"

ParticleIDPDF::ParticleIDPDF():Tool(){}

bool ParticleIDPDF::Initialise(std::string configfile, DataModel &data){

	//---------------------------------------------------------------
	//----------------- Useful header -------------------------------
	//---------------------------------------------------------------

	if(configfile!="")  m_variables.Initialise(configfile); 
	m_data= &data;

	//---------------------------------------------------------------
	//----------------- Configuration variables ---------------------
	//---------------------------------------------------------------
	
	m_variables.Get("verbosity",verbosity);
	m_variables.Get("OutputFilePDF",filename);
	m_variables.Get("UseMCTruth",use_mctruth);
	m_variables.Get("DrawOverviewPlots",draw_overview);

	std::string logmessage = "ParticleIDPDF Tool: Output Filename is "+filename+".csv";
	Log(logmessage,v_message,verbosity);

	//---------------------------------------------------------------
	//----------------- Geometry variables --------------------------
	//---------------------------------------------------------------

	m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry",geom);
	Position detector_center = geom->GetTankCentre();
	tank_center_x = detector_center.X();
	tank_center_y = detector_center.Y();
	tank_center_z = detector_center.Z();
	tank_R = geom->GetTankRadius();
	n_tank_pmts = geom->GetNumDetectorsInSet("Tank");
	n_lappds = geom->GetNumDetectorsInSet("LAPPD");
	n_mrd_pmts = geom->GetNumDetectorsInSet("MRD");
	n_veto_pmts = geom->GetNumDetectorsInSet("Veto");
  	std::map<std::string,std::map<unsigned long,Detector*> >* Detectors = geom->GetDetectors();

	//----------------------------------------------------------------------
	//-----------read in PMT x/y/z positions into vectors-------------------
	//----------------------------------------------------------------------

	tank_ymin = 9999.;
	tank_ymax = -9999.;

	for (std::map<unsigned long,Detector*>::iterator it  = Detectors->at("Tank").begin();
                                                    it != Detectors->at("Tank").end();
                                                  ++it){
    		Detector* apmt = it->second;
    		unsigned long detkey = it->first;
    		unsigned long chankey = apmt->GetChannels()->begin()->first;
    		Position position_PMT = apmt->GetDetectorPosition();
    		pmts_x.insert(std::pair<unsigned long,double>(chankey,position_PMT.X()-tank_center_x));
    		pmts_y.insert(std::pair<unsigned long,double>(chankey,position_PMT.Y()-tank_center_y));
    		pmts_z.insert(std::pair<unsigned long,double>(chankey,position_PMT.Z()-tank_center_z));
    		if (pmts_y[chankey]>tank_ymax) tank_ymax = pmts_y.at(chankey);
    		if (pmts_y[chankey]<tank_ymin) tank_ymin = pmts_y.at(chankey);
  	}

	//-----------------------------------------------------------------------
	//------- Histograms showing the variables used in the classification----
	//-----------------------------------------------------------------------

  	hist_pmtPE = new TH1F("hist_pmtPE","Single PMT Charge overview",500,0,500);
  	hist_pmtTime = new TH1F("hist_pmtTime","Single PMT Time overview",2000,0,2000);
  	hist_pmtAngle = new TH1F("hist_pmtAngle","Single PMT Angle overview",100,0,TMath::Pi());
  	hist_pmtAngle2 = new TH1F("hist_pmtAngle2","Single PMT Angle2 overview",100,0,2*TMath::Pi());
  	hist_pmtPhi = new TH1F("hist_pmtPhi","Single PMT Phi overview",100,0,2*TMath::Pi());
	hist_pmtY = new TH1F("hist_pmtY","Single PMT Y overview",100,-2.,2.);
  	hist_pmtDist = new TH1F("hist_pmtDist","Single PMT Distance overview",100,0.,4.);
  	hist_pmtHits = new TH1F("hist_pmtHits","Hit PMTs overview",150,0,150);
  	hist_pmtPEtotal = new TH1F("hist_pmtPEtotal","PMT Total Charge overview",500,0,12000);
  	hist_pmtAvgTime = new TH1F("hist_pmtAvgTime","PMT Average Time overview",100,0,50);
  	hist_pmtAvgDist = new TH1F("hist_pmtAvgDist","PMT Average Distance overview",100,0,4.);
  	hist_pmtAngleBary = new TH1F("hist_pmtAngleBary","PMT Barycenter Angle overview",100,0,TMath::Pi());
  	hist_pmtAngleRMS = new TH1F("hist_pmtAngleRMS","PMT Angular RMS overview",100,0,TMath::Pi());
  	hist_pmtAngleVar = new TH1F("hist_pmtAngleVar","PMT Angular Variance overview",100,0,TMath::Pi());
  	hist_pmtAngleSkew = new TH1F("hist_pmtAngleSkew","PMT Angular Skewness overview",100,-10,2);
  	hist_pmtAngleKurt = new TH1F("hist_pmtAngleKurt","PMT Angular Kurtosis overview",100,-10,2);
  	hist_pmtBaryRMS = new TH1F("hist_pmtBaryRMS","PMT Angular RMS (Barycenter) overview",100,0,TMath::Pi());
  	hist_pmtBaryVar = new TH1F("hist_pmtBaryVar","PMT Angular Variance (Barycenter) overview",100,0,TMath::Pi());;
 	hist_pmtBarySkew = new TH1F("hist_pmtBarySkew","Angular Skewness (Barycenter) overview",100,-10,2);
  	hist_pmtBaryKurt = new TH1F("hist_pmtBaryKurt","Angular Kurtosis (Barycenter) overview",100,-10,2);
	hist_pmtBaryFracLargeAngle = new TH1F("hist_pmtBaryFracLargeAngle","Angle (Barycenter) fraction large angles",100,0,1.);
  	hist_pmtFracRing = new TH1F("hist_pmtFracRing","Fraction of PMT hits in Ring (weighted)",100,0,1);
  	hist_pmtFracDownstream = new TH1F("hist_pmtFracDownstream","Fraction of PMT hits downstream",100,0,1);
  	hist_pmtFracRingNoWeight = new TH1F("hist_pmtFracRingNoWeight","Fraction of PMT hits in Ring",100,0,1);
  	hist_pmtFracHighestQ = new TH1F("hist_pmtFracHighestQ","Charge fraction highest PMT",100,0,1);
	hist_pmtFracClustered = new TH1F("hist_pmtFracClustered","Fraction of charge (clustered)",100,0,1);
	hist_pmtFracEarlyTime = new TH1F("hist_pmtFracEarlyTime","Fraction of early PMT hits",100,0,1.);		//t (PMT) < 4 ns
	hist_pmtFracLateTime = new TH1F("hist_pmtFracLateTime","Fraction of late PMT hits",100,0,1);		//t (PMT) > 10 ns
	hist_pmtFracLowCharge = new TH1F("hist_pmtFracLowCharge","Fraction of low charge PMT hits",100,0,1);	//Q (PMT) < 30 p.e.
	hist_pmtAngleBary_all = new TH1F("hist_pmtAngleBary_all","Angular distribution w.r.t. barycenter",100,0,TMath::Pi());
	hist_pmtAngleBary_all_ChargeWeighted = new TH1F("hist_pmtAngleBary_all_ChargeWeighted","Charge weighted angular distribution w.r.t. barycenter",100,0,TMath::Pi());
	hist_pmtAngle2Bary_all = new TH1F("hist_pmtAngle2Bary_all","Angular distribution w.r.t. barycenter (II)",100,-TMath::Pi(),TMath::Pi());
	hist_pmtAngle2Bary_all_ChargeWeighted = new TH1F("hist_pmtAngle2Bary_all_ChargeWeighted","Charge weighted angular distribution w.r.t. barycenter (II)",100,-TMath::Pi(),TMath::Pi());
	hist_pmtYBary_all = new TH1F("hist_pmtYBary_all","Y distribution w.r.t. barycenter",100,-2.,2.);
	hist_pmtYBary_all_ChargeWeighted = new TH1F("hist_pmtYBary_all_ChargeWeighted","Charge weighted Y distribution w.r.t. barycenter",100,-2.,2.);
	hist_pmtPhiBary_all = new TH1F("hist_pmtPhiBary_all","Phi distribution w.r.t. barycenter",100,-TMath::Pi(),TMath::Pi());
	hist_pmtPhiBary_all_ChargeWeighted = new TH1F("hist_pmtPhiBary_all_ChargeWeighted","Charge weighted phi distribution w.r.t. barycenter",100,-TMath::Pi(),TMath::Pi());
	hist_pmtPhiBaryRMS = new TH1F("hist_pmtPhiBaryRMS","Phi RMS (Barycenter)",100,0,TMath::Pi());
	hist_pmtPhiBaryVar = new TH1F("hist_pmtPhiBaryVar","Phi Variance (Barycenter)",100,0,2*TMath::Pi());
	hist_pmtPhiBaryFracLargeAngle = new TH1F("hist_pmtPhiBaryFracLargeAngle","Phi Large Angles (Barycenter)",100,0,1.);


  	hist_lappdPE = new TH1F("hist_lappdPE","Single LAPPD Hit Charge overview",100,0,5);
  	hist_lappdTime = new TH1F("hist_lappdTime","Single LAPPD Hit Time overview",1000,0,1000);
  	hist_lappdAngle = new TH1F("hist_lappdAngle","Single LAPPD Hit Angular overview",100,0,TMath::Pi());
  	hist_lappdDist = new TH1F("hist_lappdDist","Single LAPPD Hit Distance overview",100,0,4.);
  	hist_lappdHits = new TH1F("hist_lappdHits","Number LAPPD Hits overview",1000,0,1000);
  	hist_lappdPEtotal = new TH1F("hist_lappdPEtotal","LAPPD Total Charge overview",1000,0,1000);
  	hist_lappdAvgTime = new TH1F("hist_lappdAvgTime","LAPPD Average Time overview",100,0,50);
	hist_lappdAvgDist = new TH1F("hist_lappdAvgDist","LAPPD Average Distance overview",100,0,4.);
  	hist_lappdAngleBary = new TH1F("hist_lappdAngleBary","LAPPD Barycenter Angle overview",100,0,TMath::Pi());
  	hist_lappdAngleRMS = new TH1F("hist_lappdAngleRMS","LAPPD Angular RMS overview",100,0,TMath::Pi());
  	hist_lappdAngleVar = new TH1F("hist_lappdAngleVar","LAPPD Angular Variance overview",100,0,TMath::Pi());
  	hist_lappdAngleSkew = new TH1F("hist_lappdAngleSkew","LAPPD Angular Skewness overview",100,-2,10);
  	hist_lappdAngleKurt = new TH1F("hist_lappdAngleKurt","LAPPD Angular Kurtosis overview",100,-2,10);
  	hist_lappdBaryRMS = new TH1F("hist_lappdBaryRMS","LAPPD Angular RMS (Barycenter) overview",100,0,TMath::Pi());
  	hist_lappdBaryVar = new TH1F("hist_lappdBaryVar","LAPPD Angular Variance (Barycenter) overview",100,0,TMath::Pi());
  	hist_lappdBarySkew = new TH1F("hist_lappdBarySkew","LAPPD Angular Skewness (Barycenter) overview",100,-2,10);
  	hist_lappdBaryKurt = new TH1F("hist_lappdBaryKurt","LAPPD Angular Kurtosis (Barycenter) overview",100,-2,10);
  	hist_lappdFracRing = new TH1F("hist_lappdFracRing","Fraction of LAPPD Hits in Ring",100,0,1);
  	hist_lappdFracRingNoWeight = new TH1F("hist_lappdFracRingNoWeight","Fraction of LAPPD Hits in Ring",100,0,1);

  	hist_mrdPaddles = new TH1F("hist_mrdPaddles","Num MRD Paddles hit",310,0,310);
  	hist_mrdLayers = new TH1F("hist_mrdLayers","Num MRD Layers hit",11,0,11);
  	hist_mrdconsLayers = new TH1F("hist_mrdconsLayers","Num consecutive MRD Layers hit",11,0,11);
	hist_mrdClusters = new TH1F("hist_mrdClusters","Num MRD Clusters",20,0,20);
	hist_energy = new TH1F("hist_energy","Energy primary particle",5000,0,5000);
	hist_nrings = new TH1F("hist_nrings","Number of rings",20,0,20);
	hist_multiplerings = new TH1F("hist_multiple_rings","Multiple Rings",2,0,2);

	hist_mrdXSpread = new TH1F("hist_mrdXSpread","MRD Layer-wise X spread",100,0,5.);
	hist_mrdYSpread = new TH1F("hist_mrdYSpread","MRD Layer-wise Y spread",100,0,5.);
	hist_mrdAdjacentHits = new TH1F("hist_mrdAdjacentHits","MRD # Adjacent hits",100,0,30);
	hist_mrdPaddlesPerLayer = new TH1F("hist_mrdPaddlesPerLayer","MRD Paddles per Layer",100,0,10.);


	Log("ParticleIDPDF Tool: Initializing",v_message,verbosity);

	return true;
}


bool ParticleIDPDF::Execute(){

	Log("ParticleIDPDF Tool: Executing",v_debug,verbosity);

	//-----------------------------------------------------------------------
	//--------------------- Getting variables from stores -------------------
	//-----------------------------------------------------------------------
	
	get_ok = m_data->Stores.count("ANNIEEvent");
	if(!get_ok){
		Log("ParticleIDPDF Tool: No ANNIEEvent store!",v_error,verbosity);
		return false;
	};

	get_ok = m_data->Stores.count("RecoEvent");
	if(!get_ok){
		Log("ParticleIDPDF Tool: No RecoEvent store!",v_error,verbosity);
		return false;
	}

	get_ok = m_data->Stores.at("ANNIEEvent")->Get("EventNumber",evnum);
	if(not get_ok){ Log("ParticleIDPDF Tool: Error retrieving EventNumber,true from ANNIEEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCEventNum",mcevnum);
	if(not get_ok){ Log("ParticleIDPDF Tool: Error retrieving MCEventNum,true from ANNIEEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCHits",MCHits);
	if(not get_ok){ Log("ParticleIDPDF Tool: Error retrieving MCHits,true from ANNIEEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCLAPPDHits",MCLAPPDHits);
	if(not get_ok){ Log("ParticleIDPDF Tool: Error retrieving MCLAPPDHits,true from ANNIEEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("TDCData",TDCData);
	if(not get_ok){ Log("ParticleIDPDF Tool: Error retrieving TDCData,true from ANNIEEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("NumMrdTimeClusters",NumMrdTimeClusters);
	if(not get_ok){ Log("ParticleIDPDF Tool: Error retrieving NumMrdTimeClusters,true from ANNIEEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("RecoEvent")->Get("TrueVertex",TrueVertex);
	if(not get_ok){ Log("ParticleIDPDF Tool: Error retrieving TrueVertex,true from RecoEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("RecoEvent")->Get("TrueStopVertex",TrueStopVertex);
	if(not get_ok){ Log("ParticleIDPDF Tool: Error retrieving TrueStopVertex,true from RecoEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("RecoEvent")->Get("EventCutStatus",EventCutStatus);
	if(not get_ok){ Log("ParticleIDPDF Tool: Error retrieving EventCutStatus,true from RecoEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("RecoEvent")->Get("RecoDigit",RecoDigits);
	if(not get_ok){ Log("ParticleIDPDF Tool: Error retrieving RecoDigit,true from RecoEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("RecoEvent")->Get("TrueMuonEnergy",TrueMuonEnergy);
	if(not get_ok){ Log("ParticleIDPDF Tool: Error retrieving TrueMuonEnergy,true from RecoEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("RecoEvent")->Get("NRings",nrings);
	if(not get_ok){ Log("ParticleIDPDF Tool: Error retrieving NRings, true from RecoEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("RecoEvent")->Get("NoPiK",no_pik);
	if(not get_ok){ Log("ParticleIDPDF Tool: Error retrieving NoPiK, true from RecoEvent!",v_error,verbosity); return false; }


	//-----------------------------------------------------------------------
	//------------- Getting true info from primary particle -----------------
	//-----------------------------------------------------------------------

	//check first if the cuts were passed
	//
	bool sensible_energy = true;
	if (fabs(TrueMuonEnergy+9999)<0.01) sensible_energy = false;
	
	if (!EventCutStatus || !sensible_energy){
		Log("ParticleIDPDF Tool: EventCutStatus is false! Selected event cuts were not passed...",v_message,verbosity);
	} else {

	//std::cout <<"TrueMuonEnergy: "<<TrueMuonEnergy<<std::endl;  	
	//if yes, get the relevant information from the data
	Position pos = TrueVertex->GetPosition();
	//std::cout <<"ParticleIDPDF tool: Primary particle position (cm): ("<<pos.X()<<","<<pos.Y()<<","<<pos.Z()<<std::endl;
	pos.UnitToMeter();
	double pos_x = pos.X();
	double pos_y = pos.Y();
	double pos_z = pos.Z();
	Direction dir = TrueVertex->GetDirection();
	double dir_x = dir.X();
	double dir_y = dir.Y();
	double dir_z = dir.Z();

	//std::cout <<"ParticleIDPDF tool: Primary particle position: ("<<pos_x<<","<<pos_y<<","<<pos_z<<")"<<std::endl;
	//std::cout <<"ParticleIDPDF tool: Primary particle direction: ("<<dir_x<<","<<dir_y<<","<<dir_z<<")"<<std::endl;
	
	double distWallVert = fabs(1.98-pos_y);
	double distWallHor = 1.524 - sqrt(pos_x*pos_x+pos_z*pos_z);
	distWallVert /= 1.98;
	distWallHor /= 1.524;
	double distInnerStrHor = tank_innerR - sqrt(pos_x*pos_x+pos_z*pos_z);
	double distInnerStrVert = (pos_y > 0)? (tank_ymax - pos_y) : (pos_y - tank_ymin);
	double true_time = TrueVertex->GetTime();

	logmessage = "ParticleIDPDF Tool: Interaction Vertex is at ("+to_string(pos_x)
		    +", "+to_string(pos_y)+", "+to_string(pos_z)+")\n"
		    +"ParticleIDPDF Tool: Primary particle has energy "+to_string(TrueMuonEnergy)+"MeV and direction ("
		    +to_string(dir_x)+", "+to_string(dir_y)+", "+to_string(dir_z)+" )";
	Log(logmessage,v_debug,verbosity);


	//-----------------------------------------------------------------------
	//--------------- Read out Digits (PMT/LAPPD hits) ----------------------
	//-----------------------------------------------------------------------

	//information available directly from data
	double pmt_QDownstream=0.;
	double pmt_avgT=0.;
	Position pmtBaryQ(0.,0.,0.);
	double pmt_totalQ=0.;
	double pmt_totalQ_Clustered=0.;
	double pmt_highestQ=0.;
	int pmt_hits=0;
	int pmt_hits_late=0;
	int pmt_hits_early=0;
	int pmt_hits_lowq=0;

	Position lappdBaryQ(0.,0.,0.);
	double lappd_totalQ=0.;
	double lappd_avgT=0.;
	int lappd_hits=0;
	
	std::vector<double> pmtQ, lappdQ, pmtT, lappdT;
	std::vector<Position> pmtPos, lappdPos;

	//information only available when using mctruth information
	double pmt_rmsAngle=0.; 
	double pmt_avgDist=0.;
	double pmt_QRing=0.;
	int pmt_hitsRing=0;

	double lappd_rmsAngle=0.;
	double lappd_avgDist=0.;
	double lappd_QRing=0.;
	int lappd_hitsRing=0;
	
	std::vector<double> pmtDist, pmtAngle, pmtAngle2, pmtPhi, pmtY, lappdDist, lappdAngle;

	//std::cout <<"ParticleID PDF: RecoDigits size: "<<RecoDigits.size()<<std::endl;

	for (unsigned int i_digit = 0; i_digit < RecoDigits.size(); i_digit++){

		//std::cout <<"----------------------------------------------------"<<std::endl;
		RecoDigit thisdigit = RecoDigits.at(i_digit);
		Position detector_pos = thisdigit.GetPosition();
		detector_pos.UnitToMeter();
		//std::cout <<"detector_pos: ("<<detector_pos.X()<<","<<detector_pos.Y()<<","<<detector_pos.Z()<<")"<<std::endl;
		Direction detector_dir(detector_pos.X()-pos_x,detector_pos.Y()-pos_y,detector_pos.Z()-pos_z);
		double detector_dirX = detector_pos.X()-pos_x;
		double detector_dirY = detector_pos.Y()-pos_y;
		double detector_dirZ = detector_pos.Z()-pos_z;
		//std::cout <<"manual dir: ("<<detector_pos.X()-pos_x<<","<<detector_pos.Y()-pos_y<<","<<detector_pos.Z()-pos_z<<")"<<std::endl;
		int digittype = thisdigit.GetDigitType();		//0 - PMTs, 1 - LAPPDs	
		double digitQ = thisdigit.GetCalCharge();
		double digitT = thisdigit.GetCalTime();
		//std::cout <<"ParticleIDPDF Tool: Q = "<<digitQ<<", T = "<<digitT<<", digittype = "<<digittype<<std::endl;
		double detDist, detAngle, detAngle2, detPhi;
		if (use_mctruth){	
			//do calculations with respect to interaction vertex
			detDist = sqrt(pow(detector_pos.X()-pos_x,2)+pow(detector_pos.Y()-pos_y,2)+pow(detector_pos.Z()-pos_z,2));
			//std::cout <<"nominator for acos: "<<detector_dirX*dir_x+detector_dirY*dir_y+detector_dirZ*dir_z<<std::endl;
			detAngle = acos((detector_dirX*dir_x+detector_dirY*dir_y+detector_dirZ*dir_z)/detDist);
		} else {
			//otherwise calculations in relation to the center of the tank
			detDist = sqrt(pow(detector_pos.X(),2)+pow(detector_pos.Y(),2)+pow(detector_pos.Z(),2));
			detAngle = acos((detector_pos.Z()*1.)/detDist);
			if (fabs(detector_dirX) < 0.001 || fabs(detector_dirY) < 0.001)  detAngle2 = 0.;
			if (detector_dirX > 0. && detector_dirY > 0. ) detAngle2 = acos(detector_dirY/detector_dirX);
			if (detector_dirX < 0. && detector_dirY > 0. ) detAngle2 = TMath::Pi()/2.+acos(-detector_dirX/detector_dirY);
			if (detector_dirX < 0. && detector_dirY < 0. ) detAngle2 = TMath::Pi()+acos(detector_dirY/detector_dirX);
			if (detector_dirX > 0. && detector_dirY < 0. ) detAngle2 = 3./2*TMath::Pi()+acos(detector_dirX/-detector_dirY);
	                if (detector_dirX > 0. && detector_dirZ > 0. ) detPhi = atan(detector_dirZ/detector_dirX)+TMath::Pi()/2.;
                        if (detector_dirX > 0. && detector_dirZ < 0. ) detPhi = atan(detector_dirX/-detector_dirZ);
                        if (detector_dirX < 0. && detector_dirZ < 0. ) detPhi = 3*TMath::Pi()/2.+atan(detector_dirZ/detector_dirX);
                        if (detector_dirX < 0. && detector_dirZ > 0. ) detPhi = TMath::Pi()+atan(-detector_dirX/detector_dirZ);			
		}
		pmt_avgDist+=detDist;
		//std::cout <<"ParticleIDPDF: detDist: "<<detDist<<", detAngle: "<<detAngle<<std::endl;

		if (digittype == 0){

			hist_pmtPE->Fill(digitQ);
			hist_pmtTime->Fill(digitT);
			hist_pmtAngle->Fill(detAngle);
			hist_pmtDist->Fill(detDist);
			hist_pmtAngle2->Fill(detAngle2);
			hist_pmtPhi->Fill(detPhi);
			hist_pmtY->Fill(detector_dirY);
			pmtQ.push_back(digitQ);
			pmtT.push_back(digitT);
			pmtPos.push_back(detector_pos);
			pmt_totalQ+=digitQ;
			if (thisdigit.GetFilterStatus() == 1) pmt_totalQ_Clustered+=digitQ;
			pmt_avgT+=digitT;
			pmt_hits++;
			if (digitQ < 30) pmt_hits_lowq++;
			if (digitT > 10) pmt_hits_late++;	
			else if (digitT < 4) pmt_hits_early++;
			if (digitQ > pmt_highestQ) pmt_highestQ = digitQ;
			if (detAngle < TMath::Pi()/2.) pmt_QDownstream+=digitQ;
			pmt_rmsAngle+=(detAngle*detAngle);
			pmtDist.push_back(detDist);
			pmtAngle.push_back(detAngle);
			pmtAngle2.push_back(detAngle2);
			pmtPhi.push_back(detPhi);
			pmtY.push_back(detector_dirY);
                        pmtBaryQ += digitQ*detector_pos;

			if (use_mctruth){
                        	if (detAngle < cherenkov_angle) {
                                	pmt_QRing+=digitQ;
                                	pmt_hitsRing++;
                        	}
			}


		} else if (digittype == 1){

			hist_lappdPE->Fill(digitQ);
			hist_lappdTime->Fill(digitT);
			hist_lappdAngle->Fill(detAngle);
			hist_lappdDist->Fill(detDist);
			lappdQ.push_back(digitQ);
			lappdT.push_back(digitT);
			lappdPos.push_back(detector_pos);
			lappd_totalQ+=digitQ;
			lappd_hits++;
			lappd_avgT+=digitT;
			lappd_rmsAngle+=(detAngle*detAngle);
			lappdDist.push_back(detDist);
			lappdAngle.push_back(detAngle);
			lappdBaryQ += digitQ*detector_pos;

			if (use_mctruth){
				if (detAngle < cherenkov_angle) {
					lappd_QRing+=digitQ;
					lappd_hitsRing++;
				}
			}

		} else {
			
			std::string logmessage = "Digit type " + std::to_string(digittype) + "was not recognized. Omit reading in entry from RecoEvent store.";
			Log(logmessage,v_message,verbosity);
			
		}

	}
	
	//-----------------------------------------------------------------------
	//--------------- Calculate properties for csv file ---------------------
	//-----------------------------------------------------------------------

	double pmt_fracRing=0.;
	double pmt_fracRingNoWeight=0.;	
	double pmt_frachighestQ=0.;
	double pmt_fracQDownstream=0.;
	double pmt_fracClustered=0.;
	double pmt_fracLowQ=0.;
	double pmt_fracLate=0.;
	double pmt_fracEarly=0.;
	double lappd_fracRing=0.;
	double lappd_fracRingNoWeight=0.;

	//std::cout <<"ParticleIDPDF tool: Calculate properties for csv file: "<<std::endl;

	if (pmtQ.size()!=0) {

		pmt_avgT /= pmtQ.size();
		pmt_avgDist /= pmtQ.size();
		pmt_rmsAngle = sqrt(pmt_rmsAngle/pmtQ.size());
		pmtBaryQ = (1./pmt_totalQ)*pmtBaryQ;
		m_data->CStore.Set("pmtBaryQ",pmtBaryQ);
		pmt_fracQDownstream = pmt_QDownstream/pmt_totalQ;
		pmt_frachighestQ = pmt_highestQ/pmt_totalQ;
		std::cout <<"pmt_totalQ_Clustered = "<<pmt_totalQ_Clustered<<", pmt_totalQ = "<<pmt_totalQ;
		pmt_fracClustered = pmt_totalQ_Clustered/pmt_totalQ;
		std::cout <<", pmt_fracClustered = "<<pmt_fracClustered<<std::endl;
		pmt_fracLowQ = double(pmt_hits_lowq)/pmt_hits;
		pmt_fracLate = double(pmt_hits_late)/pmt_hits;		
		pmt_fracEarly = double(pmt_hits_early)/pmt_hits;

		if (use_mctruth){
			pmt_fracRing = pmt_QRing/pmt_totalQ;
			pmt_fracRingNoWeight = double(pmt_hitsRing)/pmt_hits;
		}
	}
	if (lappdQ.size()!=0){
		
		lappd_avgT /= lappdQ.size();
		lappd_avgDist /= lappdQ.size();
		lappd_rmsAngle = sqrt(lappd_rmsAngle/lappdQ.size());
		lappdBaryQ = (1./lappd_totalQ)*lappdBaryQ;
		if (use_mctruth){
			lappd_fracRing = lappd_QRing/lappd_totalQ;
			lappd_fracRingNoWeight = double(lappd_hitsRing)/lappd_hits;
		}
	}

	//std::cout <<"ParticleIDPDF tool: Calculate barycenter for csv file: "<<std::endl;
	double pmtBaryDist, pmtBaryAngle, pmtBaryAngle2, pmtBaryY, pmtBaryPhi, lappdBaryDist, lappdBaryAngle;

	if (use_mctruth){
		//angle and distance of barycenter calculated with respect to interaction point and primary particle direction
		Position dirBaryQ = pmtBaryQ-pos;
		pmtBaryDist = dirBaryQ.Mag();
		pmtBaryAngle = acos((dirBaryQ.X()*dir_x+dirBaryQ.Y()*dir_y+dirBaryQ.Z()*dir_z)/pmtBaryDist);
		Position lappd_dirBaryQ = lappdBaryQ-pos;
		lappdBaryDist = lappd_dirBaryQ.Mag();
		lappdBaryAngle = acos((lappd_dirBaryQ.X()*dir_x+lappd_dirBaryQ.Y()*dir_y+lappd_dirBaryQ.Z()*dir_z)/lappdBaryDist);
	} else {
		//angle and distance of barycenter calculated with respect to (0,0,0) and (0,0,1)-direction
		pmtBaryDist = pmtBaryQ.Mag();
		if (fabs(pmtBaryDist)<0.0001) {
			pmtBaryAngle = 0;
			pmtBaryAngle2 = 0.;
			pmtBaryPhi = 0.;
			pmtBaryY = 0.;
		}
		else {
			pmtBaryAngle = acos(pmtBaryQ.Z()/pmtBaryDist);
                        if (pmtBaryQ.X() > 0. && pmtBaryQ.Y() > 0. ) pmtBaryAngle2 = atan(pmtBaryQ.Y()/pmtBaryQ.X());
                        if (pmtBaryQ.X() < 0. && pmtBaryQ.Y() > 0. ) pmtBaryAngle2 = TMath::Pi()/2.+atan(-pmtBaryQ.X()/pmtBaryQ.Y());
                        if (pmtBaryQ.X() < 0. && pmtBaryQ.Y() < 0. ) pmtBaryAngle2 = TMath::Pi()+atan(pmtBaryQ.Y()/pmtBaryQ.X());
                        if (pmtBaryQ.X() > 0. && pmtBaryQ.Y() < 0. ) pmtBaryAngle2 = 3./2*TMath::Pi()+atan(pmtBaryQ.X()/-pmtBaryQ.Y());
                        if (pmtBaryQ.X() > 0. && pmtBaryQ.Z() > 0. ) pmtBaryPhi = atan(pmtBaryQ.Z()/pmtBaryQ.X())+TMath::Pi()/2.;
                        if (pmtBaryQ.X() > 0. && pmtBaryQ.Z() < 0. ) pmtBaryPhi = atan(pmtBaryQ.X()/-pmtBaryQ.Z());
                        if (pmtBaryQ.X() < 0. && pmtBaryQ.Z() < 0. ) pmtBaryPhi = 3*TMath::Pi()/2.+atan(pmtBaryQ.Z()/pmtBaryQ.X());
                        if (pmtBaryQ.X() < 0. && pmtBaryQ.Z() > 0. ) pmtBaryPhi = TMath::Pi()+atan(-pmtBaryQ.X()/pmtBaryQ.Z());
			pmtBaryY = pmtBaryQ.Y();
		}
		lappdBaryDist = lappdBaryQ.Mag();
		//std::cout <<"lappdBaryQ.Mag: "<<lappdBaryQ.Mag()<<", z component: "<<lappdBaryQ.Z()<<", acos: "<<acos(lappdBaryQ.Z()/lappdBaryDist)<<std::endl;
		if (fabs(lappdBaryDist)<0.001) lappdBaryAngle=0.;
		else lappdBaryAngle = acos(lappdBaryQ.Z()/lappdBaryDist);
	}


	double pmt_varAngle=0.;
	double pmt_skewAngle=0.;
	double pmt_kurtAngle=0.;
	double lappd_varAngle=0.;
	double lappd_skewAngle=0.;
	double lappd_kurtAngle=0.;
	double pmt_rmsAngleBary=0.;
	double pmt_varAngleBary=0.;
	double pmt_skewAngleBary=0.;
	double pmt_kurtAngleBary=0.;
	double lappd_rmsAngleBary=0.;
	double lappd_varAngleBary=0.;
	double lappd_skewAngleBary=0.;
	double lappd_kurtAngleBary=0.;

	double pmt_rmsPhiBary = 0.;
	double pmt_varPhiBary = 0.;
	double pmt_fracLargeAnglePhi=0.;
	double pmt_fracLargeAngle=0.;
	int pmt_hits_largeangle=0;
	int pmt_hits_largeangle_phi=0;	

	//std::cout <<"ParticleIDPDF tool: Calculate variance/skewness/etc for csv file: "<<std::endl;

	for (unsigned int i_pmt=0; i_pmt < pmtQ.size(); i_pmt++){
		pmt_varAngle+=(pow(pmtAngle.at(i_pmt),2)*pmtQ.at(i_pmt)/pmt_totalQ);
		pmt_skewAngle+=(pow(pmtAngle.at(i_pmt),3)*pmtQ.at(i_pmt)/pmt_totalQ);
		pmt_kurtAngle+=(pow(pmtAngle.at(i_pmt),4)*pmtQ.at(i_pmt)/pmt_totalQ);
		hist_pmtAngleBary_all->Fill(pmtAngle.at(i_pmt)-pmtBaryAngle);
		hist_pmtAngleBary_all_ChargeWeighted->Fill((pmtAngle.at(i_pmt)-pmtBaryAngle),pmtQ.at(i_pmt)/pmt_totalQ);
		double diff_angle2 = (pmtAngle2.at(i_pmt)-pmtBaryAngle2);
		if (diff_angle2 > TMath::Pi()) diff_angle2 = -(2*TMath::Pi()-diff_angle2);
		else if (diff_angle2 < -TMath::Pi()) diff_angle2 = 2*TMath::Pi()+diff_angle2;
		hist_pmtAngle2Bary_all->Fill(diff_angle2);
		hist_pmtAngle2Bary_all_ChargeWeighted->Fill(diff_angle2,pmtQ.at(i_pmt)/pmt_totalQ);
		double diff_phi = (pmtPhi.at(i_pmt)-pmtBaryPhi);
		if (diff_phi > TMath::Pi()) diff_phi = -(2*TMath::Pi()-diff_phi);
		else if (diff_phi < -TMath::Pi()) diff_phi = 2*TMath::Pi()+diff_phi;
		hist_pmtPhiBary_all->Fill(diff_phi);
		hist_pmtPhiBary_all_ChargeWeighted->Fill(diff_phi,pmtQ.at(i_pmt)/pmt_totalQ);
		hist_pmtYBary_all->Fill(pmtY.at(i_pmt)-pmtBaryY);
		hist_pmtYBary_all_ChargeWeighted->Fill(pmtY.at(i_pmt)-pmtBaryY,pmtQ.at(i_pmt));
		pmt_rmsAngleBary+=pow(pmtAngle.at(i_pmt)-pmtBaryAngle,2);
		pmt_varAngleBary+=(pow(pmtAngle.at(i_pmt)-pmtBaryAngle,2)*pmtQ.at(i_pmt)/pmt_totalQ);
		pmt_skewAngleBary+=((pow(pmtAngle.at(i_pmt)-pmtBaryAngle,3)*pmtQ.at(i_pmt)/pmt_totalQ));
		pmt_kurtAngleBary+=((pow(pmtAngle.at(i_pmt)-pmtBaryAngle,4)*pmtQ.at(i_pmt)/pmt_totalQ));
		if (fabs(pmtAngle.at(i_pmt)-pmtBaryAngle) > 0.9) pmt_hits_largeangle++;
		if (fabs(diff_phi) > 1.) pmt_hits_largeangle_phi++;
		pmt_rmsPhiBary+=(pow(diff_phi,2));
		pmt_varPhiBary+=(pow(diff_phi,2)*pmtQ.at(i_pmt)/pmt_totalQ);
	}
	pmt_kurtAngle-=(3*pmt_varAngle*pmt_varAngle);
      	pmt_kurtAngleBary-=(3*pmt_varAngleBary*pmt_varAngleBary);
	if (pmtQ.size()>0) {
		pmt_rmsAngleBary = sqrt(pmt_rmsAngleBary/pmtQ.size());
		pmt_rmsPhiBary = sqrt(pmt_rmsPhiBary/pmtQ.size());
		pmt_fracLargeAngle = double(pmt_hits_largeangle)/pmtQ.size();
		pmt_fracLargeAnglePhi = double(pmt_hits_largeangle_phi)/pmtQ.size();	
	}	

	for (unsigned int i_lappd=0; i_lappd< lappdQ.size(); i_lappd++){
		lappd_varAngle+=(pow(lappdAngle.at(i_lappd),2)*lappdQ.at(i_lappd)/lappd_totalQ);
		lappd_skewAngle+=(pow(lappdAngle.at(i_lappd),3)*lappdQ.at(i_lappd)/lappd_totalQ);
		lappd_kurtAngle+=(pow(lappdAngle.at(i_lappd),4)*lappdQ.at(i_lappd)/lappd_totalQ);
		lappd_rmsAngleBary+=pow(lappdAngle.at(i_lappd)-lappdBaryAngle,2);
		lappd_varAngleBary+=(pow(lappdAngle.at(i_lappd)-lappdBaryAngle,2)*lappdQ.at(i_lappd)/lappd_totalQ);
		lappd_skewAngleBary+=((pow(lappdAngle.at(i_lappd)-lappdBaryAngle,3)*lappdQ.at(i_lappd)/lappd_totalQ));
		lappd_kurtAngleBary+=((pow(lappdAngle.at(i_lappd)-lappdBaryAngle,4)*lappdQ.at(i_lappd)/lappd_totalQ));
	}
	lappd_kurtAngle-=(3*lappd_varAngle*lappd_varAngle);
	lappd_kurtAngleBary-=(3*lappd_varAngleBary*lappd_varAngleBary);
	if (lappdQ.size()>0) lappd_rmsAngleBary = sqrt(lappd_rmsAngleBary/lappdQ.size());

	/*std::cout <<"pmt_QDownstream: "<<pmt_QDownstream<<std::endl;
        std::cout <<"pmt_avgT: "<<pmt_avgT<<std::endl;
        std::cout <<"pmtBaryQ: ("<< pmtBaryQ.X()<<","<<pmtBaryQ.Y()<<","<<pmtBaryQ.Z()<<")"<<std::endl;
        std::cout <<"pmt_totalQ: "<<pmt_totalQ<<std::endl;
        std::cout <<"pmt_highestQ: "<<pmt_highestQ<<std::endl;
        std::cout <<"pmt_hits: "<<pmt_hits<<std::endl;
        std::cout <<"pmt_rmsAngle: "<<pmt_rmsAngle<<std::endl;
        std::cout <<"pmt_avgDist: "<<pmt_avgDist<<std::endl;
        std::cout <<"pmt_QRing: "<<pmt_QRing<<std::endl;
        std::cout <<"pmt_hitsRing: "<<pmt_hitsRing<<std::endl;
	std::cout <<"pmt_varAngleBary: "<<pmt_varAngleBary<<std::endl;
	std::cout <<"pmt_skewAngleBary: "<<pmt_skewAngleBary<<std::endl;
	std::cout <<"pmt_kurtAngleBary: "<<pmt_kurtAngleBary<<std::endl;
        std::cout <<"lappdBaryQ: ("<<lappdBaryQ.X()<<", "<<lappdBaryQ.Y()<<","<<lappdBaryQ.Z()<<")"<<std::endl;
        std::cout <<"lappd_totalQ: "<<lappd_totalQ<<std::endl;
        std::cout <<"lappd_avgT: "<<lappd_avgT<<std::endl;
        std::cout <<"lappd_hits: "<<lappd_hits<<std::endl;
        std::cout <<"lappd_rmsAngle: "<<lappd_rmsAngle<<std::endl;
        std::cout <<"lappd_avgDist: "<<lappd_avgDist<<std::endl;
        std::cout <<"lappd_QRing: "<< lappd_QRing<<std::endl;
        std::cout <<"lappd_hitsRing: "<<lappd_hitsRing<<std::endl;
*/
	//-----------------------------------------------------------------------
	//--------------- Filling properties into histograms---------------------
	//-----------------------------------------------------------------------
//	std::cout <<"ParticleIDPDF tool: Filling properties into histograms: "<<std::endl;

  	hist_pmtHits->Fill(pmt_hits);
	hist_pmtPEtotal->Fill(pmt_totalQ);
	hist_pmtAvgTime->Fill(pmt_avgT);
	hist_pmtAvgDist->Fill(pmt_avgDist);
	hist_pmtAngleBary->Fill(pmtBaryAngle);
	hist_pmtAngleRMS->Fill(pmt_rmsAngle);
	hist_pmtAngleVar->Fill(pmt_varAngle);
	hist_pmtAngleSkew->Fill(pmt_skewAngle);
	hist_pmtAngleKurt->Fill(pmt_kurtAngle);
	hist_pmtBaryRMS->Fill(pmt_rmsAngleBary);
	hist_pmtBaryVar->Fill(pmt_varAngleBary);
	hist_pmtBarySkew->Fill(pmt_skewAngleBary);
	hist_pmtBaryKurt->Fill(pmt_kurtAngleBary);
	hist_pmtPhiBaryRMS->Fill(pmt_rmsPhiBary);
	hist_pmtPhiBaryVar->Fill(pmt_varPhiBary);
	hist_pmtFracRing->Fill(pmt_fracRing);
	hist_pmtFracDownstream->Fill(pmt_fracQDownstream);
	hist_pmtFracRingNoWeight->Fill(pmt_fracRingNoWeight);
	hist_pmtFracHighestQ->Fill(pmt_frachighestQ);
	hist_pmtFracClustered->Fill(pmt_fracClustered);
	hist_pmtFracLateTime->Fill(pmt_fracLate);
	hist_pmtFracLowCharge->Fill(pmt_fracLowQ);
	hist_pmtFracEarlyTime->Fill(pmt_fracEarly);
	hist_pmtBaryFracLargeAngle->Fill(pmt_fracLargeAngle);
	hist_pmtPhiBaryFracLargeAngle->Fill(pmt_fracLargeAnglePhi);


	hist_lappdHits->Fill(lappd_hits);
	hist_lappdPEtotal->Fill(lappd_totalQ);
	hist_lappdAvgTime->Fill(lappd_avgT);
	hist_lappdAvgDist->Fill(lappd_avgDist);
	hist_lappdAngleBary->Fill(lappdBaryAngle);
	hist_lappdAngleRMS->Fill(lappd_rmsAngle);
	hist_lappdAngleVar->Fill(lappd_varAngle);
	hist_lappdAngleSkew->Fill(lappd_skewAngle);
	hist_lappdAngleKurt->Fill(lappd_kurtAngle);
	hist_lappdBaryRMS->Fill(lappd_rmsAngleBary);
	hist_lappdBaryVar->Fill(lappd_varAngleBary);
	hist_lappdBarySkew->Fill(lappd_skewAngleBary);
	hist_lappdBaryKurt->Fill(lappd_kurtAngleBary);
	hist_lappdFracRing->Fill(lappd_fracRing);
	hist_lappdFracRingNoWeight->Fill(lappd_fracRingNoWeight);

	hist_energy->Fill(TrueMuonEnergy);
	hist_nrings->Fill(nrings);
	hist_multiplerings->Fill(!no_pik);

	//-----------------------------------------------------------------------
	//--------------- Store properties in vectors ---------------------------
	//-----------------------------------------------------------------------
//	std::cout <<"ParticleIDPDF tool: Storing properties in vectors: "<<std::endl;

	if (use_mctruth){
		vec_distHor.push_back(distWallHor);
		vec_distVert.push_back(distWallVert);
		vec_distInnerHor.push_back(distInnerStrHor);
		vec_distInnerVert.push_back(distInnerStrVert);
		vec_pmt_fracRing.push_back(pmt_fracRing);
		vec_pmt_fracRingNoWeight.push_back(pmt_fracRingNoWeight);
		vec_lappd_fracRing.push_back(lappd_fracRing);
		vec_lappd_fracRingNoWeight.push_back(lappd_fracRingNoWeight);
	}           

  	vec_pmt_avgDist.push_back(pmt_avgDist);
	vec_pmt_hits.push_back(pmt_hits);
	vec_pmt_totalQ.push_back(pmt_totalQ);
	vec_pmt_avgT.push_back(pmt_avgT);
	vec_pmt_baryAngle.push_back(pmtBaryAngle);
	vec_pmt_rmsAngle.push_back(pmt_rmsAngle);
	vec_pmt_varAngle.push_back(pmt_varAngle);
	vec_pmt_skewAngle.push_back(pmt_skewAngle);
	vec_pmt_kurtAngle.push_back(pmt_kurtAngle);
	vec_pmt_rmsBary.push_back(pmt_rmsAngleBary);
	vec_pmt_varBary.push_back(pmt_varAngleBary);
	vec_pmt_skewBary.push_back(pmt_skewAngleBary);
	vec_pmt_kurtBary.push_back(pmt_kurtAngleBary);
	vec_pmt_rmsPhi.push_back(pmt_rmsPhiBary);
	vec_pmt_varPhi.push_back(pmt_varPhiBary);
	vec_pmt_fracDownstream.push_back(pmt_fracQDownstream);
	vec_pmt_fracHighestQ.push_back(pmt_frachighestQ);
	vec_pmt_highestQ.push_back(pmt_highestQ);
	vec_pmt_fracClustered.push_back(pmt_fracClustered);
	vec_pmt_fracLowQ.push_back(pmt_fracLowQ);
	vec_pmt_fracLateT.push_back(pmt_fracLate);
	vec_pmt_fracEarlyT.push_back(pmt_fracEarly);
	vec_pmt_largeangleAngle.push_back(pmt_fracLargeAngle);
	vec_pmt_largeanglePhi.push_back(pmt_fracLargeAnglePhi);

  	vec_lappd_avgDist.push_back(lappd_avgDist);
	vec_lappd_hits.push_back(lappd_hits);
	vec_lappd_totalQ.push_back(lappd_totalQ);
	vec_lappd_avgT.push_back(lappd_avgT);
	vec_lappd_baryAngle.push_back(lappdBaryAngle);
	vec_lappd_rmsAngle.push_back(lappd_rmsAngle);
	vec_lappd_varAngle.push_back(lappd_varAngle);
	vec_lappd_skewAngle.push_back(lappd_skewAngle);
	vec_lappd_kurtAngle.push_back(lappd_kurtAngle);
	vec_lappd_rmsBary.push_back(lappd_rmsAngleBary);
	vec_lappd_varBary.push_back(lappd_varAngleBary);
	vec_lappd_skewBary.push_back(lappd_skewAngleBary);
	vec_lappd_kurtBary.push_back(lappd_kurtAngleBary);

	vec_energy.push_back(TrueMuonEnergy);
	vec_evnum.push_back(evnum);
	vec_mcevnum.push_back(mcevnum);
	vec_nrings.push_back(nrings);
	vec_multiplerings.push_back(!no_pik);
	//vec_pid.push_back(pid);


        //-----------------------------------------------------------------------
        //----------------------- Read out MC - MRD hits ------------------------
        //-----------------------------------------------------------------------	


	//std::cout <<"ParticleIDPDF tool: Reading out MRD data: "<<std::endl;

  	if(!TDCData){
    		Log("ParticleIDPDF tool: No TDC data to process!",v_warning,verbosity);
  	} else {
		TH1F *x_layer[11], *y_layer[11];
		for (int i_layer=0;i_layer<11;i_layer++){
			std::stringstream ss_x_layer, ss_y_layer;
			ss_x_layer <<"hist_x_layer"<<i_layer;
			ss_y_layer <<"hist_y_layer"<<i_layer;
			x_layer[i_layer] = new TH1F(ss_x_layer.str().c_str(),ss_x_layer.str().c_str(),100,1,0);
			y_layer[i_layer] = new TH1F(ss_y_layer.str().c_str(),ss_y_layer.str().c_str(),100,1,0);
		}
    		if(TDCData->size()==0){
      			Log("ParticleIDPDF tool: No TDC hits.",v_message,verbosity);
			hist_mrdPaddles->Fill(0);
			hist_mrdLayers->Fill(0);
			hist_mrdconsLayers->Fill(0);
			vec_mrd_paddles.push_back(0);
			vec_mrd_layers.push_back(0);
			vec_mrd_conslayers.push_back(0);
			vec_mrd_adjacenthits.push_back(0); //layer-wise: how many hits do we have directly next to each other? 
			vec_mrd_xspread.push_back(0.);  //layer-wise: how much do the hits in the MRD spread out? (x-direction)
			vec_mrd_yspread.push_back(0.);	//layer-wise: how much do the hits in the MRD spread out? (y-direction)
    			vec_mrd_paddlesperlayer.push_back(0.);  //layer-wise: how many paddles are hit on average?	
	} else {
			int num_mrd_paddles=0;
			int num_mrd_layers=0;
			int num_mrd_conslayers=0;
			int num_mrd_adjacent=0;
			bool layer_occupied[11] = {0};
			double mrd_paddlesize[11];
			std::vector<std::vector<double>> mrd_hits;
			for (int i_layer=0; i_layer<11; i_layer++){
				std::vector<double> empty_hits;
				mrd_hits.push_back(empty_hits);
			}
			std::vector<int> temp_cons_layers;
			//std::cout <<"ParticleIDPDF tool: Start loop over MRD hits: "<<std::endl;
      			for(auto&& anmrdpmt : (*TDCData)){
        			unsigned long chankey = anmrdpmt.first;
        			Detector *thedetector = geom->ChannelToDetector(chankey);
				if(thedetector->GetDetectorElement()!="MRD") {
            				continue;                 // this is a veto hit, not an MRD hit.
        			}
      				num_mrd_paddles++;
				int detkey = thedetector->GetDetectorID();
				Paddle *apaddle = geom->GetDetectorPaddle(detkey);
				int layer = apaddle->GetLayer();
				layer_occupied[layer-1]=true;
				if (apaddle->GetOrientation()==1) {
					x_layer[layer-2]->Fill(0.5*(apaddle->GetXmin()+apaddle->GetXmax()));
					mrd_hits.at(layer-2).push_back(0.5*(apaddle->GetXmin()+apaddle->GetXmax()));
					mrd_paddlesize[layer-2]=apaddle->GetPaddleWidth();
				}
				else if (apaddle->GetOrientation()==0) {
					y_layer[layer-2]->Fill(0.5*(apaddle->GetYmin()+apaddle->GetYmax()));
					mrd_hits.at(layer-2).push_back(0.5*(apaddle->GetYmin()+apaddle->GetYmax()));
					mrd_paddlesize[layer-2]=apaddle->GetPaddleWidth();
				}
			}
    			hist_mrdPaddles->Fill(num_mrd_paddles);
			vec_mrd_paddles.push_back(num_mrd_paddles);
			if (num_mrd_paddles > 0) {
				for (int i_layer=0;i_layer<11;i_layer++){
					//std::cout <<"ParticleIDPDF: i_layer "<<i_layer<<std::endl;
					if (layer_occupied[i_layer]==true) {
						//std::cout <<"Layer is occupied!"<<std::endl;
						num_mrd_layers++;
						if (num_mrd_conslayers==0) num_mrd_conslayers++;
						else {
							if (layer_occupied[i_layer-1]==true) num_mrd_conslayers++;
							else {
								temp_cons_layers.push_back(num_mrd_conslayers);
								num_mrd_conslayers=0;
							}
						}
					}
					for (unsigned int i_hitpaddle=0; i_hitpaddle<mrd_hits.at(i_layer).size(); i_hitpaddle++){
						for (unsigned int j_hitpaddle= i_hitpaddle+1; j_hitpaddle < mrd_hits.at(i_layer).size(); j_hitpaddle++){
							if (fabs(mrd_hits.at(i_layer).at(i_hitpaddle)-mrd_hits.at(i_layer).at(j_hitpaddle))-mrd_paddlesize[i_layer] < 0.001) num_mrd_adjacent++;
						}
					} 		
				}
			}
			std::vector<int>::iterator it = std::max_element(temp_cons_layers.begin(),temp_cons_layers.end());
			if (it != temp_cons_layers.end()) num_mrd_conslayers = *it;
			else num_mrd_conslayers=0;
			hist_mrdLayers->Fill(num_mrd_layers);
			hist_mrdconsLayers->Fill(num_mrd_conslayers);
			vec_mrd_layers.push_back(num_mrd_layers);
			vec_mrd_conslayers.push_back(num_mrd_conslayers);
    			double padperlayer = double(num_mrd_paddles)/num_mrd_layers;
			vec_mrd_paddlesperlayer.push_back(padperlayer);
			hist_mrdPaddlesPerLayer->Fill(padperlayer);
			vec_mrd_adjacenthits.push_back(num_mrd_adjacent);
			hist_mrdAdjacentHits->Fill(num_mrd_adjacent);
			double mean_xspread=0.;
			double mean_yspread=0.;
			int num_xspread=0;
			int num_yspread=0;
			for (int i_layer=0; i_layer < 11; i_layer++){
				if (x_layer[i_layer]->GetEntries()>0){
					mean_xspread+=(x_layer[i_layer]->GetRMS());
					num_xspread++;
				}
				if (y_layer[i_layer]->GetEntries()>0){
					mean_yspread+=(y_layer[i_layer]->GetRMS());
					num_yspread++;	
				}
			}
			if (num_xspread>0) mean_xspread/=num_xspread;
			if (num_yspread>0) mean_yspread/=num_yspread;
			vec_mrd_xspread.push_back(mean_xspread);
			vec_mrd_yspread.push_back(mean_yspread);
			hist_mrdXSpread->Fill(mean_xspread);
			hist_mrdYSpread->Fill(mean_yspread);
		}
		for (int i_layer = 0; i_layer < 11; i_layer++){
			delete x_layer[i_layer];
			delete y_layer[i_layer];
		}
  	}
	hist_mrdClusters->Fill(NumMrdTimeClusters);
	vec_mrd_clusters.push_back(NumMrdTimeClusters);
	}

	Log("ParticleIDPDF tool: Execution step finished",v_message,verbosity);

	return true;

}


bool ParticleIDPDF::Finalise(){

  	Log("ParticleIDPDF tool: Finalisation started",v_message,verbosity);
  	Log("ParticleIDPDF tool: Writing overview plots to root-file.",v_message,verbosity);
  	
        //-----------------------------------------------------------------------
        //-------------- Write overview histograms to file ----------------------
        //-----------------------------------------------------------------------	

	if (draw_overview){
		std::cout <<"ParticleIDPDF tool: Writing histograms to root-file: "<<std::endl;
		std::string filename_root = filename + ".root";
		TFile *file_overview = new TFile(filename_root.c_str(),"RECREATE");
  		file_overview->cd();
  		        
        	hist_pmtPE->Write();
        	hist_pmtTime->Write();
		hist_pmtAngle->Write();
		hist_pmtAngle2->Write();
		hist_pmtPhi->Write();
		hist_pmtY->Write();
		hist_pmtDist->Write();
        	hist_pmtHits->Write();
        	hist_pmtPEtotal->Write();
        	hist_pmtAvgTime->Write();
        	hist_pmtAngleBary->Write();
        	hist_pmtAngleRMS->Write();
        	hist_pmtAngleVar->Write();
        	hist_pmtAngleSkew->Write();
        	hist_pmtAngleKurt->Write();
       		hist_pmtBaryRMS->Write();
        	hist_pmtBaryVar->Write();
        	hist_pmtBarySkew->Write();
        	hist_pmtBaryKurt->Write();
		hist_pmtPhiBaryRMS->Write();
		hist_pmtPhiBaryVar->Write();
        	hist_pmtFracRing->Write();
        	hist_pmtFracDownstream->Write();
        	hist_pmtFracRingNoWeight->Write();
        	hist_pmtFracHighestQ->Write();
		hist_pmtFracClustered->Write();
		hist_pmtFracLowCharge->Write();
		hist_pmtFracLateTime->Write();
		hist_pmtFracEarlyTime->Write();
		hist_pmtAngleBary_all->Write();
		hist_pmtAngleBary_all_ChargeWeighted->Write();
		hist_pmtAngle2Bary_all->Write();
		hist_pmtAngle2Bary_all_ChargeWeighted->Write();
		hist_pmtPhiBary_all->Write();
		hist_pmtPhiBary_all_ChargeWeighted->Write();
		hist_pmtYBary_all->Write();
		hist_pmtYBary_all_ChargeWeighted->Write();
		hist_pmtBaryFracLargeAngle->Write();
		hist_pmtPhiBaryFracLargeAngle->Write();

	
        	hist_lappdPE->Write();
        	hist_lappdTime->Write();
        	hist_lappdAngle->Write();
        	hist_lappdDist->Write();
        	hist_lappdHits->Write();
        	hist_lappdPEtotal->Write();
       		hist_lappdAvgTime->Write();
        	hist_lappdAngleBary->Write();
        	hist_lappdAngleRMS->Write();
        	hist_lappdAngleVar->Write();
        	hist_lappdAngleSkew->Write();
        	hist_lappdAngleKurt->Write();
        	hist_lappdBaryRMS->Write();
        	hist_lappdBaryVar->Write();
        	hist_lappdBarySkew->Write();
        	hist_lappdBaryKurt->Write();
        	hist_lappdFracRing->Write();
        	hist_lappdFracRingNoWeight->Write();
	
		hist_mrdPaddles->Write();
		hist_mrdLayers->Write();
		hist_mrdconsLayers->Write();
		hist_mrdClusters->Write();	

		hist_mrdPaddlesPerLayer->Write();
		hist_mrdXSpread->Write();
		hist_mrdYSpread->Write();
		hist_mrdAdjacentHits->Write();
	
		hist_energy->Write();
		hist_nrings->Write();
		hist_multiplerings->Write();

		file_overview->Close();
		delete file_overview;
	}

        //-----------------------------------------------------------------------
        //-------------- Create csv files ---------------------------------------
        //-----------------------------------------------------------------------	
	
	std::string filename_csv = filename + ".csv";
	
	if (use_mctruth) {				//use truth information like event vertex and primary energy in the calculation of variables
	
		Log("ParticleIDPDF Tool: Start writing to csv-file (with MC truth values).",v_message,verbosity);

  		ofstream csv_file(filename_csv.c_str());
		csv_file << "pmtHits,pmtTotalQ,pmtAvgT,baryAngle,rmsAngle,varAngle,skewAngle,kurtAngle,rmsBary,varBary,skewBary,kurtBary,rmsPhi,varPhi,fracHighestQ,fracDownstream,fracClustered,fracLowQ,fracLateT,fracEarlyT,largeAngleAngle,largeAnglePhi,fracRing,distVert,distHor,distInnerVert,distInnerHor,lappdHits,lappdAvgT,lappdBaryAngle,lappdRmsAngle,lappdVarAngle,lappdSkewAngle,lappdKurtAngle,lappdRMSBary,lappdVarBary,lappdSkewBary,lappdKurtBary,lappdFracRing,mrdPaddles,mrdLayers,mrdConslayers,mrdClusters,mrdPaddlesPerLayer,mrdXSpread,mrdYSpread,mrdAdjacentHits,nrings,multiplerings,energy,evnum,mcevnum"<<std::endl;
  		for (unsigned int i_vector = 0; i_vector < vec_pmt_hits.size(); i_vector++){
//			std::cout <<"i_vector (mc truth): "<<i_vector<<std::endl;
			csv_file << vec_pmt_hits.at(i_vector) 
				<<","<<vec_pmt_totalQ.at(i_vector)
				<<","<<vec_pmt_avgT.at(i_vector)
				<<","<<vec_pmt_baryAngle.at(i_vector)
				<<","<<vec_pmt_rmsAngle.at(i_vector)
				<<","<<vec_pmt_varAngle.at(i_vector)
				<<","<<vec_pmt_skewAngle.at(i_vector)
				<<","<<vec_pmt_kurtAngle.at(i_vector)
				<<","<<vec_pmt_rmsBary.at(i_vector)
				<<","<<vec_pmt_varBary.at(i_vector)
				<<","<<vec_pmt_skewBary.at(i_vector)
				<<","<<vec_pmt_kurtBary.at(i_vector)
				<<","<<vec_pmt_rmsPhi.at(i_vector)
                                <<","<<vec_pmt_varPhi.at(i_vector)
				<<","<<vec_pmt_fracHighestQ.at(i_vector)
				<<","<<vec_pmt_fracDownstream.at(i_vector)
				<<","<<vec_pmt_fracClustered.at(i_vector)
				<<","<<vec_pmt_fracLowQ.at(i_vector)
				<<","<<vec_pmt_fracLateT.at(i_vector)
				<<","<<vec_pmt_fracEarlyT.at(i_vector)
				<<","<<vec_pmt_largeangleAngle.at(i_vector)
				<<","<<vec_pmt_largeanglePhi.at(i_vector)
				<<","<<vec_pmt_fracRing.at(i_vector)
				<<","<<vec_distVert.at(i_vector)
				<<","<<vec_distHor.at(i_vector)
				<<","<<vec_distInnerVert.at(i_vector)
				<<","<<vec_distInnerHor.at(i_vector)
				<<","<<vec_lappd_hits.at(i_vector)
				<<","<<vec_lappd_avgT.at(i_vector)
				<<","<<vec_lappd_baryAngle.at(i_vector)
				<<","<<vec_lappd_rmsAngle.at(i_vector)
				<<","<<vec_lappd_varAngle.at(i_vector)
				<<","<<vec_lappd_skewAngle.at(i_vector)
				<<","<<vec_lappd_kurtAngle.at(i_vector)
				<<","<<vec_lappd_rmsBary.at(i_vector)
				<<","<<vec_lappd_varBary.at(i_vector)
				<<","<<vec_lappd_skewBary.at(i_vector)
				<<","<<vec_lappd_kurtBary.at(i_vector)
				<<","<<vec_lappd_fracRing.at(i_vector)
				<<","<<vec_mrd_paddles.at(i_vector)
				<<","<<vec_mrd_layers.at(i_vector)
				<<","<<vec_mrd_conslayers.at(i_vector)
				<<","<<vec_mrd_clusters.at(i_vector)
				<<","<<vec_mrd_paddlesperlayer.at(i_vector)
				<<","<<vec_mrd_xspread.at(i_vector)
				<<","<<vec_mrd_yspread.at(i_vector)
				<<","<<vec_mrd_adjacenthits.at(i_vector)
				<<","<<vec_nrings.at(i_vector)
				<<","<<vec_multiplerings.at(i_vector)
				<<","<<vec_energy.at(i_vector)
				<<","<<vec_evnum.at(i_vector)
				<<","<<vec_mcevnum.at(i_vector)
				<<std::endl;
  		}
  		std::string logmessage = "Cross check data sizes: PMT mctruth vectors: "+std::to_string(vec_pmt_hits.size())+", LAPPD mctruth vectors: "+std::to_string(vec_lappd_hits.size());
		Log(logmessage,v_message,verbosity);
  		
		csv_file.close();
  	} else {					//don't use truth information in the calculation of input variables

		Log("ParticleIDPDF Tool: Start writing information to csv-file.",v_message,verbosity);

  		ofstream csv_file(filename_csv.c_str());
		csv_file << "pmtHits,pmtTotalQ,pmtAvgT,pmtBaryAngle,pmtRmsAngle,pmtVarAngle,pmtSkewAngle,pmtKurtAngle,pmtRmsBary,pmtVarBary,pmtSkewBary,pmtKurtBary,pmtRmsPhi,pmtVarPhi,pmtFracHighestQ,pmtFracDownstream,pmtFracClustered,pmtFracLowQ,pmtFracLateT,pmtFracEarlyT,pmtLargeAngleAngle,pmtLargeAnglePhi,lappdHits,lappdAvgT,lappdBaryAngle,lappdRmsAngle,lappdVarAngle,lappdSkewAngle,lappdKurtAngle,lappdRMSBary,lappdVarBary,lappdSkewBary,lappdKurtBary,mrdPaddles,mrdLayers,mrdConslayers,mrdClusters,mrdPaddlesPerLayer,mrdXSpread,mrdYSpread,mrdAdjacentHits,nrings,multiplerings,energy,evnum,mcevnum"<<std::endl;
  		for (unsigned int i_vector = 0; i_vector < vec_pmt_totalQ.size(); i_vector++){
//			std::cout <<"i_vector (no mc truth): "<<i_vector<<std::endl;
       			csv_file << vec_pmt_hits.at(i_vector) 
				<<","<<vec_pmt_totalQ.at(i_vector)
				<<","<<vec_pmt_avgT.at(i_vector)
				<<","<<vec_pmt_baryAngle.at(i_vector)
				<<","<<vec_pmt_rmsAngle.at(i_vector)
				<<","<<vec_pmt_varAngle.at(i_vector)
				<<","<<vec_pmt_skewAngle.at(i_vector)
				<<","<<vec_pmt_kurtAngle.at(i_vector)
				<<","<<vec_pmt_rmsBary.at(i_vector)
				<<","<<vec_pmt_varBary.at(i_vector)
				<<","<<vec_pmt_skewBary.at(i_vector)
				<<","<<vec_pmt_kurtBary.at(i_vector)
				<<","<<vec_pmt_rmsPhi.at(i_vector)
				<<","<<vec_pmt_varPhi.at(i_vector)
				<<","<<vec_pmt_fracHighestQ.at(i_vector)
				<<","<<vec_pmt_fracDownstream.at(i_vector)
				<<","<<vec_pmt_fracClustered.at(i_vector)
				<<","<<vec_pmt_fracLowQ.at(i_vector)
				<<","<<vec_pmt_fracLateT.at(i_vector)
                                <<","<<vec_pmt_fracEarlyT.at(i_vector)
                                <<","<<vec_pmt_largeangleAngle.at(i_vector)
                                <<","<<vec_pmt_largeanglePhi.at(i_vector)
				<<","<<vec_lappd_hits.at(i_vector)
				<<","<<vec_lappd_avgT.at(i_vector)
				<<","<<vec_lappd_baryAngle.at(i_vector)
				<<","<<vec_lappd_rmsAngle.at(i_vector)
				<<","<<vec_lappd_varAngle.at(i_vector)
				<<","<<vec_lappd_skewAngle.at(i_vector)
				<<","<<vec_lappd_kurtAngle.at(i_vector)
				<<","<<vec_lappd_rmsBary.at(i_vector)
				<<","<<vec_lappd_varBary.at(i_vector)
				<<","<<vec_lappd_skewBary.at(i_vector)
				<<","<<vec_lappd_kurtBary.at(i_vector)
				<<","<<vec_mrd_paddles.at(i_vector)
				<<","<<vec_mrd_layers.at(i_vector)
				<<","<<vec_mrd_conslayers.at(i_vector)
				<<","<<vec_mrd_clusters.at(i_vector)
                                <<","<<vec_mrd_paddlesperlayer.at(i_vector)
                                <<","<<vec_mrd_xspread.at(i_vector)
                                <<","<<vec_mrd_yspread.at(i_vector)
                                <<","<<vec_mrd_adjacenthits.at(i_vector)
				<<","<<vec_nrings.at(i_vector)
				<<","<<vec_multiplerings.at(i_vector)
				<<","<<vec_energy.at(i_vector)
				<<","<<vec_evnum.at(i_vector)
				<<","<<vec_mcevnum.at(i_vector)
				<<std::endl;
  		}
  		csv_file.close();
	}
  
	Log("ParticleIDPDF Tool: Information written to csv-file.",v_message,verbosity);
	Log("ParticleIDPDF Tool: Finalisation complete",v_message,verbosity);

	return true;

}

