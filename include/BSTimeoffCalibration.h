#ifndef __BSTimeoffCal__
#define __BSTimeoffCal__

#include "G4KM2A_Geometry.h"
#include "TClonesArray.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "LHEvent.h"
#include "slamac.h"
#include "math.h"

#define MNED 5500
#define Source_Ra 83.6332
#define Source_Dec 22.0145

class BSTimeoffCalibration
{
public:
	static BSTimeoffCalibration* GetInstance()
	{
		if (m_myself == 0) { m_myself = new BSTimeoffCalibration(); }
		return m_myself;
	}
	~BSTimeoffCalibration();

private:
	BSTimeoffCalibration();
	static BSTimeoffCalibration* m_myself;
	static G4KM2A_Geometry* geom;
	static LHEvent* lhevent;
	static int ITER;
	static TH2D* h2;

public:
	int Init();
	int SetEvent(LHEvent* lhe);
	int SetA(double alpha);
	int SetC0(double C0);
	int SetITER(int iter);
	int SetPreToffset(const char* meanCalfile, const char* Pre_Toffsetfile);
	int SetToffset(const char* file);
	int SetSource(double mjd = -1, double ra = Source_Ra, double dec = Source_Dec, int flag = 0);
	int SetCore(double core_x = 0, double core_y = 0);
	int Setnparticle(TClonesArray& tHits, int np, double pe);
	int SetBadHit(TClonesArray& tHits, int np, const char* style);
	double GetTheta();
	double GetPhi();
	double GetA();
	double GetC0();
	double GetSigma();
	double GetC0_set();
	double GetMPValue(int id, double delta = 10);
	double GetRefTime(int IsSet = 0);
	TH2D* GetHist();
	double* GetstatED();
	double* GetToffsetED();
	int GetNfilter();
	int ApplyDaz();
	int ApplyPreToffset();
	int Addresidual(TClonesArray& tHits, int np, const char* style);
	int GenerateToffset(const char* style, int MPmode = 0, double delta = 10, double step = 1);
	double NormToffset();
	int spacetimefilter(TClonesArray& tHits, int np, int twind, int rwind, const char* style);
	int noisefilter(TClonesArray& tHits, int np, int twindl, int twindh, int rwind, const char* style);
	int planarfit(TClonesArray& tHits, int np, const char* style);
	int conicalfit(TClonesArray& tHits, int np, float alpha, const char* style, int IsSet = 0);
	int conicalfit_dir(TClonesArray& tHits, int np, float alpha, const char* style, int nmode = 0, double mjd = 0, double ra0 = 0, double dec0 = 0);
	int eventcallineED();
	int eventcallineED_dir();

private:
	int Nevents;
	int NfiltED;
	double theta, phi;
	double rchi;
	double rchi_sum;
	double conical_residualsED[MNED];
	double Pre_ToffsetED[MNED];
	double ToffsetED[MNED];
	double cal_statED[MNED];
	double dir[3];
	double dir_set[3];
	double core[2];
	double a;
	double c0;
	double c0_set;
	double Sigma;
	double Sigma_mean;
	double Dazcos[3];


};

#endif