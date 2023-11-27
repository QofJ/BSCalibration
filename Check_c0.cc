#include "BSTimeoffCalibration.h"
#include "LHRecEvent.h"
#include "G4KM2A_Reconstruction.h"
#include "EOSopen.h"
#include "TNtupleD.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TH1D.h"
#include "TH2D.h"
#include <stdio.h>
#include <stdlib.h>    
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>

#define DRa 1.
#define DDec 1.
double precession2Old(double mjd, double* hra, double* ddec);
int h2e_argo(double mjd, double zen, double az, double* hra, double* ddec);
int Calibration(BSTimeoffCalibration* BS, double rin, double rout, double El, double dCx, double dCy, int Nfit);

int main(int argc, char* argv[])
{
	char day[100];
	int n;
	char TDCcal[200];
	int flag;
	strcpy(day,argv[1]);
	n = atoi(argv[2]);
	strcpy(TDCcal,argv[3]);
	flag = atoi(argv[4]);	

	std::shared_ptr<LHEvent> lhevent = std::make_shared<LHEvent>();
	LHEvent* ptr = lhevent.get();
	int NTotalEntry;
	double mjd, the, phi, ra, dec, corex, corey, E, npe3, num4;
	BSTimeoffCalibration* bscal = BSTimeoffCalibration::GetInstance();
	char rootfile[200] = "root://eos01.ihep.ac.cn//eos/user/q/qijincan/brightsource/result/LHE";
	strcat(rootfile, day);
	TFile* fin = TFile::Open(rootfile);
	TTree* tree = (TTree*)fin->Get("Event");
	TNtupleD* nt = (TNtupleD*)fin->Get("nt1");
	tree->SetBranchAddress("event", &ptr);
	nt->SetBranchAddress("r_Corex", &corex);
	nt->SetBranchAddress("r_Corey", &corey);
	nt->SetBranchAddress("NpE3", &npe3);
	nt->SetBranchAddress("NuM4", &num4);
	nt->SetBranchAddress("Ra", &ra);
	nt->SetBranchAddress("Dec", &dec);
	nt->SetBranchAddress("Energy", &E);
	nt->SetBranchAddress("mjd", &mjd);
	nt->SetBranchAddress("r_Theta", &the);
	nt->SetBranchAddress("r_Phi", &phi);

	NTotalEntry = tree->GetEntriesFast();
	if (n<0 || n>NTotalEntry)
	{
		printf("the index is not in scope\n");
		return -1;
	}

	tree->GetEntry(n);
	nt->GetEntry(n);

	printf("Theta0 = %lf, Phi0 = %lf, Ra = %lf, Dec = %lf, E = %lf\n", the, phi, ra, dec, E);
	printf("Core_x = %lf, Core_y = %lf\n", corex, corey);

	int ned = lhevent->GetNhitE();
	TClonesArray* HitsE;
	HitsE = lhevent->GetHitsE();

	bscal->Init();
	bscal->SetEvent(lhevent.get());
	bscal->SetCore(corex, corey);
	if (strcmp(TDCcal, "Zero") != 0)bscal->SetToffset(TDCcal);
	
	if(flag == 0){

	bscal->SetBadHit(*HitsE, ned, "ED");
	bscal->Setnparticle(*HitsE, ned, 40.);
	int tWind, rWind;
	int NfiltE;
	tWind = 400; rWind = 100;
	NfiltE = bscal->spacetimefilter(*HitsE, ned, tWind, rWind, "ED");
	printf("NfiltE_st = %d\n", NfiltE);
	NfiltE = bscal->planarfit(*HitsE, ned, "ED");
	printf("planarfit_return = %d\n", NfiltE);
	//if (bscal->planarfit(*HitsE, ned, "ED") != 0)printf("Planarfit error!\n");
	NfiltE = bscal->conicalfit_dir(*HitsE, ned, 0.035, "ED");
	printf("conicalfit_return = %d\n", NfiltE);
	NfiltE = bscal->noisefilter(*HitsE, ned, -50, 100, rWind + 100, "ED");
	printf("NED = %d, NfiltE = %d\n", ned, NfiltE);
	if (NfiltE > 200)
	{
		bscal->noisefilter(*HitsE, ned, -50, 100, rWind + 300, "ED");
		NfiltE = bscal->conicalfit_dir(*HitsE, ned, -1, "ED");
	}
	else NfiltE = bscal->conicalfit_dir(*HitsE, ned, 0.035, "ED");
	printf("conicalfit_return = %d\n", NfiltE);
	the = bscal->GetTheta(); phi = bscal->GetPhi();
	double a = bscal->GetA();
	printf("c0_dif = %lf\n", bscal->GetC0());
	}

	//compare
	else{

	bscal->SetSource(mjd, ra, dec);
	if(bscal->eventcallineED()!=0)printf("error\n");
	printf("c0 = %lf\n", bscal->GetC0());	
	}

	printf("Theta = %lf, Phi = %lf, alpha = %lf\n", the, phi, bscal->GetA());
	return 0;
							
}
