#ifndef HighEReco_H
#define HighEReco_H

//Included in ANNIE ToolAnalysis framework, taken from TITUS software package (https://github.com/TITUSHK/ts-WChRecoSandBox) and shortened

#include "TObject.h"
#include "TH3.h"
#include "TMath.h"
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

class HighEReco : public TObject {

	public:

		static const double C_VAC;
		static const double N_REF;
		static const double C_WAT;

		double *hitPMTx;
		double *hitPMTy;
		double *hitPMTz;
		double *hitPMTDirX;
		double *hitPMTDirY;
		double *hitPMTDirZ;
		double *hitPMTTimeRes;
		int *hitPMTPEs;
		int nHitPMT;
		double *unhitPMTx;
		double *unhitPMTy;
		double *unhitPMTz;
		double *unhitPMTDirX;
		double *unhitPMTDirY;
		double *unhitPMTDirZ;
		double *unhitPMTTimeRes;
		int nUnhitPMT;
		double *hitT;
		int *hitPMT;
		int nPEs;
		TH3D * electronPhotons;
		TH3D * muonPhotons;
		TH3D * electronTimes;
		TH3D * muonTimes;

		HighEReco();
		~HighEReco();

		void LikelihoodFit(double &trackCorrection, double &recoVtxX, double &recoVtxY,
				double &recoVtxZ, double &recoT, double &recoDirPhi,
				double &recoDirTheta, double &recoKE, double &recoLnL, int ipnu);
		double ElectronLnLikelihood(const double *par);
		double MuonLnLikelihood(const double *par);
		double ExpectedPMTPhotoelectrons(double vtxX, double vtxY, double vtxZ, double dirX, double dirY, double dirZ, int pmt, bool isHit, double kineticEnergy, int ipnu, double * expectedTime = 0);

	private:

		double PMTlnLikelihood(const double& vtxX,const double& vtxY,const double& vtxZ,const double& dirX,const double& dirY,const double& dirZ,const int& pmt,const bool& isHit,const double& kineticEnergy,const int& ipnu);
		double LnLikelihood(const double *par, int ipnu, bool total);
		double ElectronLnLikelihoodTotal(const double *par);
		double MuonLnLikelihoodTotal(const double *par);

//		ClassDef(HighEReco,0)

};


inline double HighEReco::PMTlnLikelihood(const double& vtxX,const double& vtxY,const double& vtxZ,const double& dirX,const double& dirY,const double& dirZ,const int& pmt,const bool& isHit,const double& kineticEnergy,const int& ipnu){
double expectedPhotoelectrons = ExpectedPMTPhotoelectrons(vtxX, vtxY, vtxZ, dirX, dirY, dirZ, pmt, isHit, kineticEnergy, ipnu);
double observedPhotoelectrons = hitPMTPEs[pmt];
double l = observedPhotoelectrons* TMath::Log(expectedPhotoelectrons)-expectedPhotoelectrons;
return l;
}




#endif