#include "HighEReco.h"
//#include <TMinuit.h>

//Included in ANNIE ToolAnalysis framework, taken from TITUS software package (https://github.com/TITUSHK/ts-WChRecoSandBox) and shortened

using namespace std;

//ClassImp(HighEReco)

HighEReco::HighEReco()
{

}

HighEReco::~HighEReco()
{

}

const double HighEReco::C_VAC = 29.9792458;
const double HighEReco::N_REF = 1.34;
const double HighEReco::C_WAT = C_VAC/N_REF;

double HighEReco::ExpectedPMTPhotoelectrons(double vtxX, double vtxY, double vtxZ, double dirX, double dirY, double dirZ, int pmt, bool isHit, double kineticEnergy, int ipnu, double *expectedTime){
      
	TH3D * photonsLookup = (ipnu == 11) ? (TH3D*)electronPhotons : (TH3D*)muonPhotons;

	TH3D * timesLookup = (ipnu == 11) ? (TH3D*)electronTimes : (TH3D*)muonTimes;
	double vtxToPMTx, vtxToPMTy, vtxToPMTz;
	if (isHit) {
		vtxToPMTx = hitPMTx[pmt] - vtxX;
		vtxToPMTy = hitPMTy[pmt] - vtxY;
		vtxToPMTz = hitPMTz[pmt] - vtxZ;
	}
	else {
		vtxToPMTx = unhitPMTx[pmt] - vtxX;
		vtxToPMTy = unhitPMTy[pmt] - vtxY;
		vtxToPMTz = unhitPMTz[pmt] - vtxZ;
	}

    double pmtDist = sqrt(pow(vtxToPMTx,2)+pow(vtxToPMTy,2)+pow(vtxToPMTz,2));
	double pmtTrackAngle = TMath::ACos((vtxToPMTx*dirX+vtxToPMTy*dirY+vtxToPMTz*dirZ)/pmtDist);
	double minDist = photonsLookup->GetXaxis()->GetBinCenter(1);
	if(pmtDist <= minDist) pmtDist = minDist+0.0001;
	double maxDist = photonsLookup->GetXaxis()->GetBinCenter(photonsLookup->GetNbinsX());
	if(pmtDist >= maxDist) pmtDist = maxDist-0.0001;
	double minAngle = photonsLookup->GetYaxis()->GetBinCenter(1);
	if(pmtTrackAngle <= minAngle) pmtTrackAngle = minAngle+0.0001;
	double maxAngle = photonsLookup->GetYaxis()->GetBinCenter(photonsLookup->GetNbinsY());
	if(pmtTrackAngle >= maxAngle) pmtTrackAngle = maxAngle-0.0001;
	double minKE = photonsLookup->GetZaxis()->GetBinCenter(1);
	if(kineticEnergy <= minKE) kineticEnergy= minKE+0.0001;
	double maxKE = photonsLookup->GetZaxis()->GetBinCenter(photonsLookup->GetNbinsZ());
	if(kineticEnergy >= maxKE) kineticEnergy = maxKE-0.0001;


	double expectedPhotons = photonsLookup->GetBinContent(photonsLookup->FindBin(pmtDist, pmtTrackAngle, kineticEnergy));//*100000000;//Remove this factor when likelihood tables normalisation is fixed
	if(expectedTime){
		*expectedTime = timesLookup->Interpolate(pmtDist, pmtTrackAngle, kineticEnergy);
	}	
	
	/*
	//////////   static const double QE = 2.*0.22; 											//m. nieslony: pmt Q.E. already included in WCSim simulation
	double QE = 2.*0.22; //WChSandBox has half as many photons so need so double QE
	//Account for solid angle covered by PMT at given distance
	////////static double PMTsize = 12*2.54/2.; //12 inch diameter 
	double PMTsize = 12*2.54/2.; //12 inch diameter
	double solidAngleFactor = (1.-pmtDist/ TMath::Sqrt(pmtDist*pmtDist+PMTsize*PMTsize)); // solid angle covered by PMT
 	/////////////   const double halfAngleBinWidth = photonsLookup->GetYaxis()->GetBinWidth(photonsLookup->GetYaxis()->FindFixBin(pmtTrackAngle))/2.;
	double halfAngleBinWidth = photonsLookup->GetYaxis()->GetBinWidth(photonsLookup->GetYaxis()->FindFixBin(pmtTrackAngle))/2.;
	solidAngleFactor /= 2*TMath::Sin(pmtTrackAngle)*TMath::Sin(halfAngleBinWidth); //solid angle covered by bin in theta
	//Account for reduction of solid angle covered by PMT not directly facing towards vertex
  	/////////////const double cosPMTangle = isHit ?
               //          -(vtxToPMTx* hitPMTDirX[pmt] + vtxToPMTy* hitPMTDirY[pmt] + vtxToPMTz* hitPMTDirZ[pmt])/pmtDist :
               ///          -(vtxToPMTx* unhitPMTDirX[pmt] + vtxToPMTy* unhitPMTDirY[pmt] + vtxToPMTz* unhitPMTDirZ[pmt])/pmtDist;
 	////////////   const double expectedPhotoelectrons = TMath::Max(expectedPhotons*QE*solidAngleFactor*cosPMTangle,1e-5);
	double cosPMTangle = isHit ?
		-(vtxToPMTx* hitPMTDirX[pmt] + vtxToPMTy* hitPMTDirY[pmt] + vtxToPMTz* hitPMTDirZ[pmt])/pmtDist :
		-(vtxToPMTx* unhitPMTDirX[pmt] + vtxToPMTy* unhitPMTDirY[pmt] + vtxToPMTz* unhitPMTDirZ[pmt])/pmtDist;
	double expectedPhotoelectrons = expectedPhotons*QE*solidAngleFactor*cosPMTangle;*/

	double expectedPhotoelectrons = expectedPhotons;
	if(expectedPhotoelectrons<0.00001) expectedPhotoelectrons = 0.00001;
	return expectedPhotoelectrons;
}


double HighEReco::LnLikelihood(const double *par, int ipnu, bool total) {
	//cerr<<"Likelihood"<<endl;        
	double vtxZ = par[1];
	double vtxY = par[2];
	double vtxX = par[3];
	double time = par[4];
	double dirCosTheta = par[5];
	double dirSinTheta = TMath::Sqrt(1-dirCosTheta*dirCosTheta);
	double dirPhi = par[6];
	double dirX = dirSinTheta* TMath::Cos(dirPhi);
	double dirY = dirSinTheta* TMath::Sin(dirPhi);
	double dirZ = dirCosTheta;
	double kineticEnergy = par[7];
	double trackOffset = par[0];
	vtxZ -= dirZ*trackOffset;
	vtxY -= dirY*trackOffset;
	vtxX -= dirX*trackOffset;
	if(vtxX*vtxX+vtxZ*vtxZ>1.534*1.534 || (vtxY)>1.9 || (vtxY)<-2.1
			|| vtxX!=vtxX || vtxY!=vtxY || vtxZ!=vtxZ || time!=time || dirCosTheta!=dirCosTheta
			|| dirPhi!=dirPhi || kineticEnergy!=kineticEnergy || trackOffset!=trackOffset
	  ){
		return -9999999999;
	}
	double lnLikelihood = 0;
	double *hitExpectedPE = new double[nHitPMT];
	double *hitExpectedT = new double[nHitPMT];
	double totalExpectedPE = 0;
	double totalObservedPE = 0;
	//Hit PMTs
	for(int iPMT = 0; iPMT< nHitPMT; iPMT++){
		hitExpectedPE[iPMT] = ExpectedPMTPhotoelectrons(vtxX,vtxY,vtxZ,dirX,dirY,dirZ,iPMT,true,kineticEnergy,ipnu,hitExpectedT + iPMT);
		totalExpectedPE += hitExpectedPE[iPMT];
		totalObservedPE += hitPMTPEs[iPMT];
	}
	//Unhit PMTs
	double *unhitExpectedPE = new double[nUnhitPMT];
	for(int iPMT = 0; iPMT< nUnhitPMT; iPMT++){
		unhitExpectedPE[iPMT] = ExpectedPMTPhotoelectrons(vtxX,vtxY,vtxZ,dirX,dirY,dirZ,iPMT,false,kineticEnergy,ipnu);		//added time variable because it would not compile otherwise
	}

	int n;
	n=nUnhitPMT;

	for(int iPMT = nUnhitPMT-n; iPMT<nUnhitPMT; iPMT++) //use PMTs with larger expected PE than nth element
		totalExpectedPE += unhitExpectedPE[iPMT];
	//if(print) cout << totalExpectedPE << " " << totalObservedPE << endl;
	if(!total){
		double normalisation = 1./totalExpectedPE;
		for(int iPMT = 0; iPMT< nHitPMT; iPMT++){
			//Poisson:
			//lnLikelihood += hitPMTPEs[iPMT]*TMath::Log(hitExpectedPE[iPMT]);

			//Multinomial:
			lnLikelihood += hitPMTPEs[iPMT]*TMath::Log(hitExpectedPE[iPMT]*normalisation);

			//Negative Binomial, r=2:
			//const double r = 1000;
			//double p = 1/(1+r/(normalisation*hitExpectedPE[iPMT]));
			//lnLikelihood += hitPMTPEs[iPMT]*TMath::Log(p)+r*TMath::Log(1-p);
		}
		//Poisson:
		//lnLikelihood -= totalExpectedPE;
		for(int iPE = 0; iPE<nPEs; iPE++){
			//Gaussian for time:
			const double sigma = 4.;
			lnLikelihood -= 0.5*TMath::Power((hitT[iPE]-hitExpectedT[hitPMT[iPE]])/(sigma+hitPMTTimeRes[hitPMT[iPE]]),2);
		}
	} else {				//m. nieslony: not used anyway, so just leave it in? --> absolute numbers are optimised for TITUS, so will inevitably be wrong for ANNIE
		//Total PE term
		//if(ipnu == 13 && totalExpectedPE>600) totalExpectedPE = 2.*totalExpectedPE - 600.;				//correction terms, need to figure out on my own
		//else if(ipnu==11 && totalExpectedPE>500/0.6) totalExpectedPE = 1.6*totalExpectedPE - 500.;
		lnLikelihood += totalObservedPE* TMath::Log(totalExpectedPE) - totalExpectedPE;
	}
	delete[] hitExpectedPE;
	delete[] hitExpectedT;
	delete[] unhitExpectedPE;
	return lnLikelihood;
}

double HighEReco::ElectronLnLikelihood(const double *par){
	//cerr<<"ElectronLnLikelihood"<<endl;        
	return -LnLikelihood(par, 11, false);
}

double HighEReco::MuonLnLikelihood(const double *par){
	//    for(int i =0; i<8; i++) cout << par[i] << " ";
	//cerr<<"MuonLnLikelihood"<<endl;        
	double result = -LnLikelihood(par, 13, false);
	//    cout << result << endl;
	return result;
}

double HighEReco::ElectronLnLikelihoodTotal(const double *par){
	//cerr<<"ElectronLnLikelihood"<<endl;        
	return -LnLikelihood(par, 11, true);
}

double HighEReco::MuonLnLikelihoodTotal(const double *par){
	return -LnLikelihood(par, 13, true);
}


void HighEReco::LikelihoodFit(double &trackCorrection, double &recoVtxX, double &recoVtxY,
		double &recoVtxZ, double &recoT, double &recoDirPhi,
		double &recoDirTheta, double &recoKE, double &recoLnL, int ipnu) {

	ROOT::Math::Minimizer*minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
	minimizer->SetMaxFunctionCalls(1000000);
	//    minimizer->SetTolerance(0.0001);
	minimizer->SetPrintLevel(0);

	ROOT::Math::Functor func;
	if(ipnu==11) {
		func = ROOT::Math::Functor(this,&HighEReco::ElectronLnLikelihood,8);
	}
	else{
		func = ROOT::Math::Functor(this,&HighEReco::MuonLnLikelihood,8);
	}
	minimizer->SetFunction(func);
	double recoDirCosTheta = TMath::Cos(recoDirTheta);

	if(recoKE < 100) recoKE = 101;			//bring variable values back into suitable range
	if(recoKE > 2000) recoKE = 1999;		
	minimizer->SetLimitedVariable(0,"trackOffset",trackCorrection,10,-1000,1000);
	minimizer->SetLimitedVariable(1,"vtxZ",recoVtxZ,1.,-1.55,1.55);
	minimizer->SetLimitedVariable(2,"vtxY",recoVtxY,1.,-2.0,2.0);
	minimizer->SetLimitedVariable(3,"vtxX",recoVtxX,1.,-1.55,1.55);
	minimizer->SetLimitedVariable(4,"T",recoT,1.,-100,100);
	minimizer->SetLimitedVariable(5,"dirCosTheta", recoDirCosTheta,0.1,-0.99999,0.99999);
	minimizer->SetLimitedVariable(6,"dirPhi",recoDirPhi,0.1,-TMath::Pi(), TMath::Pi());
	minimizer->SetLimitedVariable(7,"kineticEnergy",recoKE,10,50,3000);

	minimizer->FixVariable(1);
	minimizer->FixVariable(2);
	minimizer->FixVariable(3);
	minimizer->FixVariable(4);
	minimizer->FixVariable(5);
	minimizer->FixVariable(6);
	minimizer->FixVariable(7);

	minimizer->ReleaseVariable(5);
	minimizer->ReleaseVariable(6);
	minimizer->Minimize();
	minimizer->ReleaseVariable(1);
	minimizer->ReleaseVariable(2);
	minimizer->ReleaseVariable(3);
	minimizer->FixVariable(0);
	minimizer->Minimize();

	recoVtxZ = minimizer->X()[1];
	recoVtxY = minimizer->X()[2];
	recoVtxX = minimizer->X()[3];
	recoT = minimizer->X()[4]- minimizer->X()[0]/C_VAC;
	recoDirCosTheta = minimizer->X()[5];
	recoDirPhi = minimizer->X()[6];
	recoDirTheta= TMath::ACos(recoDirCosTheta);
	recoKE = minimizer->X()[7];
	recoVtxZ -= minimizer->X()[0]*recoDirCosTheta;
	recoVtxY -= minimizer->X()[0]* TMath::Sin(recoDirTheta)* TMath::Sin(recoDirPhi);
	recoVtxX -= minimizer->X()[0]* TMath::Sin(recoDirTheta)* TMath::Cos(recoDirPhi);
	
	double recoPar[8]={0,recoVtxZ,recoVtxY,recoVtxX,recoT,recoDirCosTheta,recoDirPhi, recoKE};
	recoLnL = LnLikelihood(recoPar, ipnu, false);

}

