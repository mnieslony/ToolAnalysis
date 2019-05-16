#include "FindTrackLengthInWater.h"

FindTrackLengthInWater::FindTrackLengthInWater():Tool(){}


bool FindTrackLengthInWater::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////
  // get configuration variables for this tool 
  m_variables.Get("verbosity",verbosity);
  m_variables.Get("WriteTrackLengthToFile",writefile);
  m_variables.Get("Outputfile",myfile);
  csvfile.open(myfile);   

      std::cout<<" open file.. max number of hits: "<<maxhits0<<std::endl;  
      if(maxhits0>1100){ 
        std::cerr<<" Please change the dim of double lambda_vec[1100]={0.}; double digitt[1100]={0.}; from 1100 to max number of hits"<<std::endl; 
      }
      //--- write to file: ---//
      //if(first==1 && deny_access==0){
      //    deny_access=1;
          for (int i=0; i<maxhits0;++i){
             stringstream strs;
             strs << i;
             string temp_str = strs.str();
             string X_name= "l_";
             X_name.append(temp_str);
             const char * xname = X_name.c_str();
             csvfile<<xname<<",";
          }
          for (int ii=0; ii<maxhits0;++ii){
             stringstream strs4;
             strs4 << ii;
             string temp_str4 = strs4.str();
             string T_name= "T_";
             T_name.append(temp_str4);
             const char * tname = T_name.c_str();
             csvfile<<tname<<",";
          }
        csvfile<<"lambda_max"<<","; //first estimation of track length(using photons projection on track)
        csvfile<<"totalPMTs"<<",";
        csvfile<<"totalLAPPDs"<<",";
        csvfile<<"lambda_max"<<",";
        csvfile<<"TrueTrackLengthInWater"<<",";
        csvfile<<"neutrinoE"<<",";
        csvfile<<"trueKE"<<",";
        csvfile<<"diffDirAbs"<<",";
        csvfile<<"TrueTrackLengthInMrd"<<",";
        csvfile<<"recoDWallR"<<",";
        csvfile<<"recoDWallZ"<<",";
        csvfile<<"dirX"<<",";
        csvfile<<"dirY"<<",";
        csvfile<<"dirZ"<<",";
        csvfile<<"vtxX"<<",";
        csvfile<<"vtxY"<<",";
        csvfile<<"vtxZ";
        csvfile<<'\n';
      // }
  
  get_ok = m_data->Stores.at("ANNIEEvent")->Header->Get("AnnieGeometry",anniegeom);
  if(not get_ok){
    Log("FindTrackLengthInWater Tool: No Geometry in ANNIEEvent!",v_error,verbosity);
    return false;
  }
  totalPMTs   = anniegeom->GetNumTankPMTs();
  totalLAPPDs = anniegeom->GetNumLAPPDs();
  tank_radius = anniegeom->GetTankRadius()*100.;
  tank_halfheight = anniegeom->GetTankHalfheight()*100.;
  
  return true;
}

bool FindTrackLengthInWater::Execute(){
   
   // read the input hit and reconstruction info
   // ==========================================
   Int_t recoStatus;
   double vtxX,vtxY,vtxZ,dirX,dirY,dirZ,TrueTrackLengthInMrd,TrueTrackLengthInWater,TrueNeutrinoEnergy,trueEnergy;
   std::string TrueInteractionType;
   std::vector<double> digitX; std::vector<double> digitY;  std::vector<double> digitZ;
   std::vector<double> digitT;
   std::map<unsigned long,std::vector<MCHit>>* MCHits=nullptr;
   
   // Get hits from the ANNIEEvent
   get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCHits", MCHits);   // FIXME change to 'Hits' to support data
   if(not get_ok){
      Log("FindTrackLengthInWater Tool: Failed to retrieve the MCHits!",v_error,verbosity);
      return false;
   }
   
   // Get reconstructed vertex from the RecoEvent
   RecoVertex theExtendedVertex;
   get_ok = m_data->Stores.at("RecoEvent")->Get("ExtendedVertex", theExtendedVertex);
   if(not get_ok){
   	Log("FindTrackLengthInWater Tool: Failed to retrieve the ExtendedVertex from RecoEvent Store!",v_error,verbosity);
   	return false;
   }
   // get the reconstructed vertex and direction
   recoStatus = theExtendedVertex.GetStatus();
   vtxX = theExtendedVertex.GetPosition().X();
   vtxY = theExtendedVertex.GetPosition().Y();
   vtxZ = theExtendedVertex.GetPosition().Z();
   dirX = theExtendedVertex.GetDirection().X();
   dirY = theExtendedVertex.GetDirection().Y();
   dirZ = theExtendedVertex.GetDirection().Z();
   
   // Get neutrino info from GenieEvent TODO FIXME
   get_ok = m_data->Stores.at("GenieEvent")->Get("TrueNeutrinoEnergy", TrueNeutrinoEnergy);
   if(not get_ok){
   	Log("FindTrackLengthInWater Tool: Failed to retrieve TrueNeutrinoEnergy!",v_error,verbosity);
   	return false;
   }
   get_ok = m_data->Stores.at("GenieEvent")->Get("TrueInteractionType", TrueInteractionType);
   if(not get_ok){
   	Log("FindTrackLengthInWater Tool: Failed to retrieve the TrueInteractionType!",v_error,verbosity);
   	return false;
   }
   
   // Get the primary muon information
   // ================================
   // find muon
   int PrimaryMuonIndex;
   get_ok = m_data->Stores.at("ANNIEEvent")->Get("PrimaryMuonIndex",PrimaryMuonIndex);
   std::vector<MCParticle>* MCParticles=nullptr;
   get_ok &= m_data->Stores.at("ANNIEEvent")->Get("MCParticles", MCParticles);
   MCParticle* primarymuon=nullptr;
   if((get_ok==0) || (PrimaryMuonIndex<0)){
     Log("FindTrackLengthInWater Tool: No PrimaryMuonIndex in ANNIEEvent",v_error,verbosity); // FIXME for data?
     return false;
   } else {
     primarymuon = &(MCParticles->at(PrimaryMuonIndex));
   }
   
   // Get info
   TrueTrackLengthInWater = primarymuon->GetTrackLengthInTank();
   TrueTrackLengthInMrd = primarymuon->GetTrackLengthInMrd();
   trueEnergy = primarymuon->GetStartEnergy();
   
   // Get the hit information
   // =======================
   int nhits = MCHits->size();
	Log("TotalLightMap Tool: Looping over PMTs with a hit",v_debug,verbosity);
	for(std::pair<const unsigned long,std::vector<MCHit>>& nextpmt : *MCHits ){
		// if it's not a tank PMT, ignore it
		Detector* thepmt = anniegeom->ChannelToDetector(nextpmt.first);
		if(thepmt->GetDetectorElement()!="Tank") continue;
		if(thepmt->GetTankLocation()=="OD") continue;  // don't utilize OD pmts?
		Position PMT_position = thepmt->GetPositionInTank();
		// loop over hits on this PMT
		// ==========================
		for(MCHit& nexthit : nextpmt.second){
			double hit_time = nexthit.GetTime();
			digitT.push_back(hit_time);
		}
		digitX.clear(); digitY.clear(); digitZ.clear();
		digitX.resize(digitT.size(), PMT_position.X());
		digitY.resize(digitT.size(), PMT_position.Y());
		digitZ.resize(digitT.size(), PMT_position.Z());
	}
   
   // Declare the variables involved in track length estimation
   // =========================================================
   float recoDWallR2=0; float recoDWallZ2=0;
   float diffDirAbs2=0;
   double lambda_min = 10000000;  double lambda_max = -99999999.9; double lambda = 0;
   double lambda_vec[1100]={0.}; double digitt[1100]={0.};
   
   // Estimate the track length in the tank
   // =====================================
   uint32_t EventNumber;
   get_ok = m_data->Stores.at("ANNIEEvent")->Get("EventNumber", EventNumber);
   std::cout<<"EventNumber: "<<EventNumber<<endl;
   if(recoStatus == 0){ count1++;
     if((TrueInteractionType == "QES - Weak[CC]") && TrueTrackLengthInMrd>0.){
   	//std::cout<<"current entry: "<<EventNumber<<" with nhits: "<<nhits<<std::endl;

        //calculate diff dir with (0,0,1)  
        double diffDirAbs0 = TMath::ACos(dirZ)*TMath::RadToDeg();
        //cout<<"diffDirAbs0: "<<diffDirAbs0<<endl;    
        diffDirAbs2=diffDirAbs0/90.;
        double recoVtxR2 = vtxX*vtxX + vtxZ*vtxZ;//vtxY*vtxY;
        double recoDWallR = tank_radius-TMath::Sqrt(recoVtxR2);   // FIXME is this coordinate-system
        double recoDWallZ = tank_halfheight-TMath::Abs(vtxY);     // dependant? Is it subtracting the tank
        recoDWallR2      = recoDWallR/tank_radius;                // origin??
        recoDWallZ2      = recoDWallZ/tank_halfheight;
        


	// Estimate the track length
	// =========================
	for(int k=0; k<nhits; k++){
          //std::cout<<"k: "<<k<<", "<<digitT.at(k)<<std::endl;
          digitt[k]=digitT.at(k);

	  // Estimate length as the distance between the reconstructed vertex last photon emission point
          lambda = find_lambda(vtxX,vtxY,vtxZ,dirX,dirY,dirZ,digitX.at(k),digitY.at(k),digitZ.at(k),42.);
          if( lambda <= lambda_min ){
	      lambda_min = lambda;
	  }
	  if( lambda >= lambda_max ){
	      lambda_max = lambda;
	  }
          lambda_vec[k]=lambda;
         //m_data->Stores["ANNIEEvent"]->Set("WaterRecoTrackLength",lambda_max);
  	}
        //std::cout<<"the track length in the water tank (1st approx) is: "<<lambda_max<<std::endl;
       
       // Post-processing of variables to store
       // =====================================
       recoDWallR2      = recoDWallR/tank_radius;
       recoDWallZ2      = recoDWallZ/tank_halfheight;
       float lambda_max_2     = TMath::Abs(lambda_max)/500;
       float totalPMTs2=double(totalPMTs)/1000.;
       float totalLAPPDs2=double(totalLAPPDs)/1000.;
       float TrueTrackLengthInWater2 = TrueTrackLengthInWater/500.;
       float TrueTrackLengthInMrd2 = TrueTrackLengthInMrd/200.;
       
       // Write to .csv file - including variables for track length & energy reconstruction
       // ==================
        for(int i=0; i<maxhits0;++i){
           csvfile<<lambda_vec[i]<<",";
        }
        for(int i=0; i<maxhits0;++i){
           csvfile<<digitt[i]<<",";
        }
        csvfile<<lambda_max<<",";
        csvfile<<totalPMTs<<",";
        csvfile<<totalLAPPDs<<",";
        csvfile<<lambda_max<<",";
        csvfile<<TrueTrackLengthInWater<<",";
        csvfile<<TrueNeutrinoEnergy<<",";
        csvfile<<trueEnergy<<",";
        csvfile<<diffDirAbs2<<",";
        csvfile<<TrueTrackLengthInMrd2<<",";
        csvfile<<recoDWallR2<<",";
        csvfile<<recoDWallZ2<<",";
        csvfile<<dirX<<",";
        csvfile<<dirY<<",";
        csvfile<<dirZ<<",";
        csvfile<<vtxX<<",";
        csvfile<<vtxY<<",";
        csvfile<<vtxZ;
        csvfile<<'\n';
        
        // Put these variables in the ANNIEEvent
        // =====================================
        m_data->Stores.at("EnergyReco")->Set("lambda_vec",lambda_vec);
        m_data->Stores.at("EnergyReco")->Set("digitt",digitt);
        m_data->Stores.at("EnergyReco")->Set("lambda_max",lambda_max);
        m_data->Stores.at("EnergyReco")->Set("TrueTrackLengthInWater",TrueTrackLengthInWater);
        m_data->Stores.at("EnergyReco")->Set("trueNeuE",TrueNeutrinoEnergy);
        m_data->Stores.at("EnergyReco")->Set("trueE",trueEnergy);
        m_data->Stores.at("EnergyReco")->Set("diffDirAbs2",diffDirAbs2);
        m_data->Stores.at("EnergyReco")->Set("TrueTrackLengthInMrd2",TrueTrackLengthInMrd2);
        m_data->Stores.at("EnergyReco")->Set("recoDWallR2",recoDWallR2);
        m_data->Stores.at("EnergyReco")->Set("recoDWallZ2",recoDWallZ2);
        m_data->Stores.at("EnergyReco")->Set("dirVec",primarymuon->GetStartDirection());
        m_data->Stores.at("EnergyReco")->Set("vtxVec",primarymuon->GetStartVertex());
        
     }
   }

  return true;
}

//---------------------------------------------------------------------------------------------------------//
// Find Distance between Muon Vertex and Photon Production Point (lambda) 
//---------------------------------------------------------------------------------------------------------//
double FindTrackLengthInWater::find_lambda(double xmu_rec,double ymu_rec,double zmu_rec,double xrecDir,double yrecDir,double zrecDir,double x_pmtpos,double y_pmtpos,double z_pmtpos,double theta_cher)
{
     double lambda1 = 0.0;    double lambda2 = 0.0;    double length = 0.0 ;
     double xmupos_t1 = 0.0;  double ymupos_t1 = 0.0;  double zmupos_t1 = 0.0;
     double xmupos_t2 = 0.0;  double ymupos_t2 = 0.0;  double zmupos_t2 = 0.0;
     double xmupos_t = 0.0;   double ymupos_t = 0.0;   double zmupos_t = 0.0;

     double theta_muDir_track = 0.0;
     double theta_muDir_track1 = 0.0;  double theta_muDir_track2 = 0.0;
     double cos_thetacher = cos(theta_cher*TMath::DegToRad());
     double xmupos_tmin = 0.0; double ymupos_tmin = 0.0; double zmupos_tmin = 0.0;
     double xmupos_tmax = 0.0; double ymupos_tmax = 0.0; double zmupos_tmax = 0.0;
     double lambda_min = 10000000;  double lambda_max = -99999999.9;  double lambda = 0.0;

     double alpha = (xrecDir*xrecDir + yrecDir*yrecDir + zrecDir*zrecDir) * ( (xrecDir*xrecDir + yrecDir*yrecDir + zrecDir*zrecDir) - (cos_thetacher*cos_thetacher) );
     double beta = ( (-2)*(xrecDir*(x_pmtpos - xmu_rec) + yrecDir*(y_pmtpos - ymu_rec) + zrecDir*(z_pmtpos - zmu_rec) )*((xrecDir*xrecDir + yrecDir*yrecDir + zrecDir*zrecDir) - (cos_thetacher*cos_thetacher)) );
     double gamma = ( ( (xrecDir*(x_pmtpos - xmu_rec) + yrecDir*(y_pmtpos - ymu_rec) + zrecDir*(z_pmtpos - zmu_rec))*(xrecDir*(x_pmtpos - xmu_rec) + yrecDir*(y_pmtpos - ymu_rec) + zrecDir*(z_pmtpos - zmu_rec)) ) - (((x_pmtpos - xmu_rec)*(x_pmtpos - xmu_rec) + (y_pmtpos - ymu_rec)*(y_pmtpos - ymu_rec) + (z_pmtpos - zmu_rec)*(z_pmtpos - zmu_rec))*(cos_thetacher*cos_thetacher)) );


     double discriminant = ( (beta*beta) - (4*alpha*gamma) );

     lambda1 = ( (-beta + sqrt(discriminant))/(2*alpha));
     lambda2 = ( (-beta - sqrt(discriminant))/(2*alpha));

     xmupos_t1 = xmu_rec + xrecDir*lambda1;  xmupos_t2 = xmu_rec + xrecDir*lambda2;
     ymupos_t1 = ymu_rec + yrecDir*lambda1;  ymupos_t2 = ymu_rec + yrecDir*lambda2;
     zmupos_t1 = zmu_rec + zrecDir*lambda1;  zmupos_t2 = zmu_rec + zrecDir*lambda2;

     double tr1 = sqrt((x_pmtpos - xmupos_t1)*(x_pmtpos - xmupos_t1) + (y_pmtpos - ymupos_t1)*(y_pmtpos - ymupos_t1) + (z_pmtpos - zmupos_t1)*(z_pmtpos - zmupos_t1));
     double tr2 = sqrt((x_pmtpos - xmupos_t2)*(x_pmtpos - xmupos_t2) + (y_pmtpos - ymupos_t2)*(y_pmtpos - ymupos_t2) + (z_pmtpos - zmupos_t2)*(z_pmtpos - zmupos_t2));
     theta_muDir_track1 = (acos( (xrecDir*(x_pmtpos - xmupos_t1) + yrecDir*(y_pmtpos - ymupos_t1) + zrecDir*(z_pmtpos - zmupos_t1))/(tr1) )*TMath::RadToDeg());
     theta_muDir_track2 = (acos( (xrecDir*(x_pmtpos - xmupos_t2) + yrecDir*(y_pmtpos - ymupos_t2) + zrecDir*(z_pmtpos - zmupos_t2))/(tr2) )*TMath::RadToDeg());

     //---------------------------- choose lambda!! ---------------------------------
     if( theta_muDir_track1 < theta_muDir_track2 ){
       lambda = lambda1;
       xmupos_t = xmupos_t1;
       ymupos_t = ymupos_t1;
       zmupos_t = zmupos_t1;
       theta_muDir_track=theta_muDir_track1;
     }else if( theta_muDir_track2 < theta_muDir_track1 ){
       lambda = lambda2;
       xmupos_t = xmupos_t2;
       ymupos_t = ymupos_t2;
       zmupos_t = zmupos_t2;
       theta_muDir_track=theta_muDir_track2;
     }

     return lambda;
}

bool FindTrackLengthInWater::Finalise(){
  csvfile.close();
  return true;
}
