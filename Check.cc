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
#include "slamac.h"

//#include "G4KM2A_Reconstruction.h"

#define DRa 1.
#define DDec 1.
double precession2Old(double mjd, double* hra, double* ddec);
int h2e_argo(double mjd, double zen, double az, double* hra, double* ddec);
int Calibration(BSTimeoffCalibration* BS, double rin, double rout, double El, double dCx, double dCy, int Nfit);

//.I ./include
//.L ./src/G4KM2A_Geometry.cc
//.L ./src/LHEvent.cc
//.L ./src/LHRecEvent.cc
//.L ./src/EOSopen.cc
//.L ./src/G4KM2A_Reconstruction.cc
//.L ./src/BSTimeoffCalibration.cc
//.L Check.cc

void GenerateOffTime()
{
	G4KM2A_Geometry* geom = G4KM2A_Geometry::GetInstance(7);
	double x, y, z;
	FILE* f = fopen("./config/pre_TimeOffset_abs.txt", "w");
	fprintf(f, "Type ID Pre_TimeOffset\n");
	for (int i = 1; i < 5500; i++)
	{
		if (geom->GetEDxyz(i, x, y, z, 1) > 0)
		{
			fprintf(f, "ED %d %lf\n", i, fabs(y) / 50.);
			//fprintf(f, "ED %d %lf\n", i, gRandom->Gaus(0, 3));
		}
	}
	fclose(f);
}

TGraph2D* ShowTimeOffset(const char* calfile)
{
	int id;
	double x, y, z, time;
	int p = 0;
	FILE* f = fopen(calfile, "r");
	G4KM2A_Geometry* geom = G4KM2A_Geometry::GetInstance(7);
	//TCanvas* c = new TCanvas("c1");
	TGraph2D* graph = new TGraph2D;

	fscanf(f, "%*[^\n]%*c");
	while (!feof(f))
	{
		fscanf(f, "ED %d %lf\n", &id, &time);
		if (geom->GetEDxyz(id, x, y, z, 1))
		{
			graph->SetPoint(p++, x, y, time);
		}
	}
	for (int i = -12; i < 15; i++)
	{
		for (int j = -15; j < 16; j++)
		{
			graph->SetPoint(p++, cos(30*DD2R) * double(i * 10) + cos(120*DD2R) * double(j * 10), sin(30*DD2R) * double(i * 10) + sin(120*DD2R) * double(j * 10), 0.);
			//if (abs(i) == 15 || abs(j) == 13)graph->SetPoint(p++, double(i * 10), double(j * 10), double(j * 10) / 50);
		}
	}

	gStyle->SetPalette(kRainBow);
	gStyle->SetOptTitle(1);
	graph->SetTitle(";x;y;z");
	//printf("%p", (graph->GetXaxis()));
	graph->Draw("surf1");
	//printf("%s\n", graph->GetXaxis()->GetTitle());
	graph->GetXaxis()->SetTitle("x (m)");
	graph->GetYaxis()->SetTitle("y (m)");
	graph->GetZaxis()->SetTitle("Preset time offset (ns)");
	graph->GetXaxis()->SetTitleSize(0.04);
	graph->GetXaxis()->SetTitleFont(42);
	graph->GetXaxis()->SetRangeUser(-200,300);
	graph->GetXaxis()->CenterTitle();
	

	//TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);

	//graph->Draw("surf1");
	//c->Update();

	return graph;
}

void ShowCol(double Cx, double Cy, BSTimeoffCalibration* BS, TClonesArray& tHits, int np, int mode = 0)
{
	LHHit* tHit;
	G4KM2A_Geometry* geom = G4KM2A_Geometry::GetInstance(7);
	TGraph* g_fit = new TGraph();
	TGraph* g_sur = new TGraph();
	double corex, corey;
	double x, y, z, dx, dy, dz, dr, dt;
	double dirl, dirm, dirn, alpha, C0, C0_set;
	dirl = sin(BS->GetTheta()) * cos(BS->GetPhi());
	dirm = sin(BS->GetTheta()) * sin(BS->GetPhi());
	dirn = cos(BS->GetTheta());
	alpha = BS->GetA();
	C0 = BS->GetC0();
	C0_set = BS->GetC0_set();
	corex = Cx; corey = Cy;

	g_fit->SetMarkerStyle(97);
	g_fit->SetMarkerColor(4);
	g_fit->GetXaxis()->SetRangeUser(0, 300);
	g_fit->GetYaxis()->SetRangeUser(-20, 60);
	g_fit->GetXaxis()->SetTitle("X (m)");
	g_fit->GetXaxis()->CenterTitle();
	g_fit->GetXaxis()->SetLabelSize(0.05);
	g_fit->GetXaxis()->SetTitleSize(0.05);
	g_fit->GetYaxis()->SetTitle("Y (m)");
	g_fit->GetYaxis()->CenterTitle();
	g_fit->GetYaxis()->SetLabelSize(0.05);
	g_fit->GetYaxis()->SetTitleSize(0.05);

	g_sur->SetMarkerStyle(29);
	g_sur->SetMarkerColor(2);
	g_sur->SetMarkerSize(2);

	int nup = 0;
	for (int i = 0; i < np; i++)
	{
		tHit = (LHHit*)((tHits)[i]);
		if (tHit->GetStatus() < 2.5 || tHit->GetPe() < 1.e-3)continue;
		// printf("\n%d %d ",i,tHit->GetId());
		if (geom->Getxyz(tHit->GetId(), x, y, z, 1, "ED") < 0)continue;
		
		dx = x - corex;
		dy = y - corey;
		dz = z;
		dr = sqrt(dx * dx + dy * dy + dz * dz - pow(dx * dirl + dy * dirm - dz * dirn, 2.));
		if (mode == 0) {
			dt = tHit->GetTime() - (BS->GetToffsetED())[tHit->GetId()] - (dirl * x + dirm * y - dirn * z + alpha * dr + C0_set) / 0.2998;
		}
		else
		{
			dt = tHit->GetTime() - (BS->GetToffsetED())[tHit->GetId()] - (dirl * x + dirm * y - dirn * z + alpha * dr + C0) / 0.2998;
		}
		//if (dt > 5 || dt < -5)continue;
		g_fit->SetPoint(nup++, dr, alpha * dr + dt * 0.2998);
	}

	g_sur->SetPoint(0, 0, 0);
	g_sur->SetPoint(1, 100, 100 * alpha);
	g_sur->SetPoint(2, 200, 200 * alpha);
	printf("alpha = %lf\n", alpha);
	
	g_fit->Draw("ap");
	g_sur->Draw("csame");
	g_sur->GetXaxis()->SetRangeUser(0, 300);
	g_sur->GetYaxis()->SetRangeUser(-20, 60);

}

//TGraph2D* DrawEDArray(int id_l, int id_h)
//{
//	TGraph2D* g = new TGraph2D;
//	auto geom = G4KM2A_Geometry::GetInstance(7);
//	int p = 0;
//	double x, y, z;
//	for (int i = id_l; i < id_h; i++)
//	{
//		if (geom->GetEDId2(i) < 0)continue;
//		geom->GetEDxyz(i, x, y, z, 1);
//		g->SetPoint(p++, x, y);
//	}
//	return g;
//}

//.I ./include
//.L ./src/G4KM2A_Geometry.cc
//.L ./src/LHEvent.cc
//.L ./src/LHRecEvent.cc
//.L ./src/EOSopen.cc
//.L ./src/G4KM2A_Reconstruction.cc
//.L ./src/BSTimeoffCalibration.cc
//.L Check.cc
void REC_LHE(const char* infile, const char* outfile)
{
	int i, j, k, n, m, id, nentries, nTotalHits, mode, tag, nflag, nFile, Nout, arrayFlag, Num[4];
	double t, pe, x, y, z, dt, mjd, pd, peda, pedd, Tresolutioni, ntu[30], Dnum[10];
	float dtt, daa, drr;
	char sTemp[80], name[50];
	arrayFlag = 7;
	printf("Flag=%d\n", arrayFlag);
	printf("***********************************************\n");
	printf("***  LHAASO reconstruction v2.1 @2021-06-16 ***\n");
	printf("***  Please email to chensz@ihep.ac.cn for  ***\n");
	printf("***  any problem or bug !                   ***\n");
	printf("***********************************************\n");
	//init the detector position firstly !!!
	G4KM2A_Geometry* geom = G4KM2A_Geometry::GetInstance(arrayFlag);
	//Calibration class 
	
	//output file and tree
	LHRecEvent* rec = new LHRecEvent();

	G4KM2A_Reconstruction* reconstruct = G4KM2A_Reconstruction::GetInstance(arrayFlag);

	//decode data event
	
	//rec event
	LHEvent* lhevent = new LHEvent();
	//output root file
	TFile* newfile = EOSopen(outfile, 1); //0 open, 1: created
	TNtupleD* nt1 = new TNtupleD("nt1", "LHAASO Reconstructed Data",
		"Nevent:mjd:NhitE:NfiltE:NhitM:NfiltM:NpE1:NpE2:NpE3:NuM1:NuM2:NuM3:NuM4:NuM5:r_Corex:r_Corey:r_Theta:r_Phi:r_C0:a:chi:NtrigE:NtrigE2:size:age:dr:nMD50");
	int nED, nMD, nEDgood, nMDgood;

	//import the files and reconstructed them

	TFile* hfile = EOSopen(infile, 0);
	printf("inputing file -> %s\n", infile);
	if (!hfile || hfile->IsZombie() || hfile->GetEND() < 50) {
		printf("file error %s!!\n", infile);
		return;
	}
	//input file and tree
	TTree* EventTree = (TTree*)gDirectory->Get("Event");
	EventTree->SetBranchAddress("event", &lhevent);
	nentries = Int_t(EventTree->GetEntriesFast());
	printf("Total event %d\n", nentries);

	Nout = int(nentries / 10.);
	//Nout=10000;
	for (i = 0; i < 1000; i++) {
		//        for(i=0; i<nentries; i++) {
		EventTree->GetEntry(i);
		reconstruct->eventrecline(lhevent, rec, 10);
		ntu[0] = lhevent->GetEvN();
		ntu[1] = lhevent->GetMjd();
		ntu[2] = lhevent->GetNhitE();
		ntu[3] = rec->NfiltE;
		ntu[4] = lhevent->GetNhitM();
		ntu[5] = rec->NfiltM;
		ntu[6] = rec->NpE1;
		ntu[7] = rec->NpE2;
		ntu[8] = rec->NpE3;
		ntu[9] = rec->NuM1;
		ntu[10] = rec->NuM2;
		ntu[11] = rec->NuM3;
		ntu[12] = rec->NuM4;
		ntu[13] = rec->NuW1;
		ntu[14] = rec->rec_x;
		ntu[15] = rec->rec_y;
		ntu[16] = rec->rec_theta;
		ntu[17] = rec->rec_phi;
		ntu[18] = rec->rec_c0;
		ntu[19] = rec->rec_a;
		ntu[20] = rec->rec_sigma;
		ntu[21] = rec->NtrigE;
		ntu[22] = rec->NtrigW;
		ntu[23] = rec->rec_Esize;
		ntu[24] = rec->rec_Eage;
		ntu[25] = rec->NpW;
		ntu[26] = rec->NhitW;
		nt1->Fill(ntu);
		if (i % Nout == 0 || i == nentries - 1)printf("---------->Reconstructed event %d %d\n", i, k);

	}
	
	hfile->Close();
	
	printf("----------End of Globle Reconstruction--------------\n");
	newfile->Write();
	newfile->Close();
}

//.I ./include
//.L ./src/G4KM2A_Geometry.cc
//.L ./src/LHEvent.cc
//.L ./src/LHRecEvent.cc
//.L ./src/EOSopen.cc
//.L ./src/G4KM2A_Reconstruction.cc
//.L ./src/BSTimeoffCalibration.cc
//.L Check.cc
//gSystem->AddLinkedLibs("-L./lib/ -l:slalib64.a")


//TH1D* Reconstruct(const char* day, int n, const char* TDCcal = "./config/TrueCal_mean_211001_221001.txt")
void Reconstruct(const char* day, int n, const char* TDCcal = "./config/TrueCal_mean_211001_221001.txt")
{
	std::shared_ptr<LHEvent> lhevent = std::make_shared<LHEvent>();
	LHEvent* ptr = lhevent.get();
	int NTotalEntry;
	double mjd, the, phi, ra, dec, corex, corey, E, npe3, num4, c0;
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
	nt->SetBranchAddress("r_C0", &c0);

	NTotalEntry = tree->GetEntriesFast();
	/*if (n<0 || n>NTotalEntry)
	{
		printf("the index is not in scope\n");
		return;
	}*/

	double sigma = 0;
	//TH1D* h1 = new TH1D("h1", "c0 - c0_rec", 600, -30, 30);
	for (int i = 0; i < n; i++) {
		tree->GetEntry(i);
		nt->GetEntry(i);

		//tree->GetEntry(n);
		//nt->GetEntry(n);

		//printf("Theta0 = %lf, Phi0 = %lf, Ra = %lf, Dec = %lf, E = %lf\n", the, phi, ra, dec, E);
		//printf("Core_x = %lf, Core_y = %lf\n", corex, corey);

		int ned = lhevent->GetNhitE();
		TClonesArray* HitsE;
		HitsE = lhevent->GetHitsE();

		bscal->Init();
		bscal->SetEvent(lhevent.get());
		bscal->SetCore(corex, corey);
		if (strcmp(TDCcal, "Zero") != 0)bscal->SetToffset(TDCcal);
		bscal->SetBadHit(*HitsE, ned, "ED");
		bscal->Setnparticle(*HitsE, ned, 40.);
		int tWind, rWind;
		int NfiltE;
		tWind = 400; rWind = 100;
		NfiltE = bscal->spacetimefilter(*HitsE, ned, tWind, rWind, "ED");
		//printf("NfiltE_st = %d\n", NfiltE);
		NfiltE = bscal->planarfit(*HitsE, ned, "ED");
		//printf("planarfit_return = %d\n", NfiltE);
		//if (bscal->planarfit(*HitsE, ned, "ED") != 0)printf("Planarfit error!\n");
		NfiltE = bscal->conicalfit_dir(*HitsE, ned, 0.035, "ED");
		//printf("conicalfit_return = %d\n", NfiltE);
		NfiltE = bscal->noisefilter(*HitsE, ned, -50, 100, rWind + 100, "ED");
		//printf("NED = %d, NfiltE = %d\n", ned, NfiltE);
		if (NfiltE > 200)
		{
			bscal->noisefilter(*HitsE, ned, -50, 100, rWind + 300, "ED");
			NfiltE = bscal->conicalfit_dir(*HitsE, ned, -1, "ED");
		}
		else NfiltE = bscal->conicalfit_dir(*HitsE, ned, 0.035, "ED");
		sigma += pow(bscal->GetSigma(), 2);
		//printf("conicalfit_return = %d\n", NfiltE);
		the = bscal->GetTheta(); phi = bscal->GetPhi();
		double a = bscal->GetA();
		if (fabs(c0 - bscal->GetC0()) > 50)printf("n(delta_c0 > 50) = %d\n", n);
		if (i == n - 1)ShowCol(corex, corey, bscal, *HitsE, ned, 1);
		//bscal->SetSource(mjd, ra + 0.5, dec + 0.5, 1); //bscal->SetC0(c0);
		//auto c2 = new TCanvas();
		//ShowCol(corex, corey, bscal, *HitsE, ned);
		//auto c3 = new TCanvas();
		//bscal->SetSource(mjd, ra - 0.5, dec - 0.5, 1);
		//ShowCol(corex, corey, bscal, *HitsE, ned);
		//h1->Fill(c0 - bscal->GetC0());
		//printf("c0_dir = %lf\n", bscal->GetC0());
		//printf("r_C0 = %lf\n", c0);
	}

	printf("Sigma = %lf\n", sigma);
	
	//return h1;
	//compare
	/*bscal->SetSource(mjd, ra, dec);
	bscal->eventcallineED();
	printf("c0 = %lf\n", bscal->GetC0());*/
	
	//printf("Theta = %lf, Phi = %lf, alpha = %lf\n", the, phi, bscal->GetA());
	
}

//int main(int argc, char* argv[])
//{
//	char day[100], TDCcal[100];
//	strcpy(day, argv[1]);
//	strcpy(TDCcal, argv[2]);
//
//	std::shared_ptr<LHEvent> lhevent = std::make_shared<LHEvent>();
//	LHEvent* ptr = lhevent.get();
//	int NTotalEntry;
//	double mjd, the, phi, ra, dec, corex, corey, E, npe3, num4;
//	BSTimeoffCalibration* bscal = BSTimeoffCalibration::GetInstance();
//	char rootfile[200] = "root://eos01.ihep.ac.cn//eos/user/q/qijincan/brightsource/result/LHE";
//	strcat(rootfile, day);
//	TFile* fin = TFile::Open(rootfile);
//	TTree* tree = (TTree*)fin->Get("Event");
//	TNtupleD* nt = (TNtupleD*)fin->Get("nt1");
//	tree->SetBranchAddress("event", &ptr);
//	nt->SetBranchAddress("r_Corex", &corex);
//	nt->SetBranchAddress("r_Corey", &corey);
//	nt->SetBranchAddress("NpE3", &npe3);
//	nt->SetBranchAddress("NuM4", &num4);
//	nt->SetBranchAddress("Ra", &ra);
//	nt->SetBranchAddress("Dec", &dec);
//	nt->SetBranchAddress("Energy", &E);
//	nt->SetBranchAddress("mjd", &mjd);
//	nt->SetBranchAddress("r_Theta", &the);
//	nt->SetBranchAddress("r_Phi", &phi);
//
//	NTotalEntry = tree->GetEntriesFast();
//	/*if (n<0 || n>NTotalEntry)
//	{
//		printf("the index is not in scope\n");
//		return;
//	}
//	
//	tree->GetEntry(n);
//	nt->GetEntry(n);*/
//
//	for (int i = 0; i < NTotalEntry; i++)
//	{
//		tree->GetEntry(i);
//		nt->GetEntry(i);
//		printf("Ra0 = %lf, Dec0 = %lf\n", ra, dec);
//		double ra0 = ra;
//		double dec0 = dec;
//		h2e_argo(mjd, the * DR2D, phi * DR2D, &ra, &dec);
//		precession2Old(mjd, &ra, &dec);
//		ra0 = ra; dec0 = dec;
//		printf("MJD = %lf\n", mjd);
//		printf("Theta0 = %lf, Phi0 = %lf, Ra = %lf, Dec = %lf, E = %lf\n", the, phi, ra, dec, E);
//		printf("Core_x = %lf, Core_y = %lf\n", corex, corey);
//
//		int ned = lhevent->GetNhitE();
//		TClonesArray* HitsE;
//		HitsE = lhevent->GetHitsE();
//
//		bscal->Init();
//		bscal->SetEvent(lhevent.get());
//		bscal->SetCore(corex, corey);
//		if (strcmp(TDCcal, "Zero") != 0)bscal->SetToffset(TDCcal);
//		bscal->SetBadHit(*HitsE, ned, "ED");
//		bscal->Setnparticle(*HitsE, ned, 40.);
//		int tWind, rWind;
//		int NfiltE;
//		tWind = 400; rWind = 100;
//		bscal->spacetimefilter(*HitsE, ned, tWind, rWind, "ED");
//		if (bscal->planarfit(*HitsE, ned, "ED") != 0)continue;
//		if (bscal->conicalfit_dir(*HitsE, ned, 0.035, "ED") != 0)continue;
//		NfiltE = bscal->noisefilter(*HitsE, ned, -50, 100, rWind + 100, "ED");
//		if (NfiltE > 200)
//		{
//			bscal->noisefilter(*HitsE, ned, -50, 100, rWind + 300, "ED");
//			bscal->conicalfit_dir(*HitsE, ned, -1, "ED");
//		}
//		else bscal->conicalfit_dir(*HitsE, ned, 0.035, "ED");
//		the = bscal->GetTheta(); phi = bscal->GetPhi();
//		double a = bscal->GetA();
//		if (the < 0)continue;
//		printf("NED = %d, NfiltE = %d\n", ned, NfiltE);
//		h2e_argo(mjd, the * DR2D, phi * DR2D, &ra, &dec);
//		precession2Old(mjd, &ra, &dec);
//		printf("Theta = %lf, Phi = %lf, Ra = %lf, Dec = %lf\na = %lf\n", the, phi, ra, dec, a);
//
//		if (pow(ra - ra0, 2.) + pow(dec - dec0, 2.) > 5)
//		{
//			printf("Index = %d\n", i);
//			printf("Theta_f = %lf, Phi_f = %lf\n", ra, dec);
//			break;
//		}
//	}
//	return 1;
//}

//.I ./include
//.L ./src/G4KM2A_Geometry.cc
//.L ./src/LHEvent.cc
//.L ./src/BSTimeoffCalibration.cc
//.L Check.cc

void FitOneSignal(double Energy, int iter = 20, int Nfit = 1)
{
	int i, j, k, ITER = iter;
	int ID[5500] = { 0 };
	double Toff[5500] = { 0 };
	int id_temp;
	double toff_temp;
	FILE* f = fopen("/eos/user/q/qijincan/brightsource/result/TrueCal_mean_211001_221001.txt", "r");
	fscanf(f, "%*[^\n]%*c");
	i = 0;
	while (!feof(f))
	{
		fscanf(f, "ED %d %lf\n", &id_temp, &toff_temp);
		if (id_temp < 0 || id_temp >= 5500)continue;
		ID[id_temp] = 1;
		Toff[id_temp] = toff_temp;
	}
	fclose(f);

	TMatrixD Result(MNED, 5);
	Result.Zero();

	double* count, * toff;
	double rin, rout, El, dCenterx, dCentery;
	rin = 0;
	rout = 0.1;
	El = Energy;
	dCenterx = 0;
	dCentery = 0;
	BSTimeoffCalibration* bscal = BSTimeoffCalibration::GetInstance();

	for (i = 0; i < ITER; i++)
	{
		bscal->SetITER(i);
		bscal->Init();
		if (i == 0)bscal->SetToffset("./config/TrueCal_mean_211001_221001.txt");
		Calibration(bscal, rin, rout, El, dCenterx, dCentery, Nfit);
		bscal->GenerateToffset("ED");
		bscal->NormToffset();
		printf("ITER:%d\n", i + 1);

	}

	count = bscal->GetstatED();
	toff = bscal->GetToffsetED();

	double mean = 0;
	j = 0; k = 0;
	for (i = 0; i < MNED; i++)
	{
		if (count[i] == 0 )continue;
		Result(j, 0) = (double)i;
		Result(j, 1) = toff[i];
		if (ID[i] == 1)
		{
			Result(j, 3) = Toff[i];
			mean += Toff[i];
			k++;
		}

		j++;
	}
	mean /= k;
	for (i = 0; i < j; i++)
	{
		Result(i, 2) = Result(i, 1) - mean;
		Result(i, 4) = Result(i, 2) - Result(i, 3);
	}

	Result.ResizeTo(j, 5);
	Result.Print();

	return;
}

//int main()
//{
//	int i, j, NTotalEntry, n_event = 0;
//	double mjd, ra, dec, E, corex, corey, npe3, num4;
//	double the_temp, phi_temp, ra_temp, dec_temp;
//	LHEvent* lhevent = new LHEvent();
//	char filehead[200] = "root://eos01.ihep.ac.cn//eos/user/q/qijincan/brightsource/result/LHE";
//	char file[200];
//	char suffix[20] = ".root";
//	const char* filedate;
//	BSTimeoffCalibration* bscal = BSTimeoffCalibration::GetInstance();
//
//	TFile* histfile = TFile::Open("root://eos01.ihep.ac.cn//eos/user/q/qijincan/brightsource/result/LHEhist2.root", "recreate");
//	TH2D* h2 = new TH2D("h2", "Ra & Dec Map", 250., 70., 95., 250., 10., 35.);
//	TH2D* h2_m = new TH2D("h2_MeanCal", "Ra & Dec Map", 250., 70., 95., 250., 10., 35.);
//	TH2D* h2_z = new TH2D("h2_ZeroCal", "Ra & Dec Map", 250., 70., 95., 250., 10., 35.);
//	h2->SetDirectory(histfile); h2_m->SetDirectory(histfile); h2_z->SetDirectory(histfile);
//
//	//for (i = 2021202; i < 2022250; i++)
//	for (i = 2021202; i < 2021210; i++)
//	{
//		if (i % 1000 == 366)
//		{
//			i = i + 1000 - 366;
//			continue;
//		}
//
//		filedate = std::to_string(i).c_str();
//		strcpy(file, filehead);
//		strcat(file, filedate);
//		strcat(file, suffix);
//
//		TFile* infile = TFile::Open(file);
//		if (!infile) {
//			std::cout << "Can't open file:" << file << std::endl;
//			continue;
//		}
//		if (infile->IsZombie())
//		{
//			std::cout << file << " ";
//			printf("file error!!\n");
//			infile->Close();
//			continue;
//		}
//
//		TTree* tree = (TTree*)infile->Get("Event");
//		if (!tree)
//		{
//			printf("get event error!!\n");
//			infile->Close();
//			continue;
//		}
//
//		TNtupleD* nt = (TNtupleD*)infile->Get("nt1");
//		if (!nt)
//		{
//			printf("get nt1 error!!\n");
//			infile->Close();
//			continue;
//		}
//
//		tree->SetBranchAddress("event", &lhevent);
//		nt->SetBranchAddress("r_Corex", &corex);
//		nt->SetBranchAddress("r_Corey", &corey);
//		nt->SetBranchAddress("NpE3", &npe3);
//		nt->SetBranchAddress("NuM4", &num4);
//		nt->SetBranchAddress("Ra", &ra);
//		nt->SetBranchAddress("Dec", &dec);
//		nt->SetBranchAddress("Energy", &E);
//		nt->SetBranchAddress("mjd", &mjd);
//
//		NTotalEntry = tree->GetEntriesFast();
//
//		for (j = 0; j < NTotalEntry; j++)
//		{
//			nt->GetEntry(j);
//			if (pow(10., E) < 10)continue;
//			if (log10((num4 + 0.0001) / npe3) > -2.3)continue;
//			tree->GetEntry(j);
//
//			h2->Fill(ra, dec);
//			//printf("Ra = %lf, Dec = %lf\n", ra, dec);
//
//			int ned = lhevent->GetNhitE();
//			TClonesArray* HitsE;
//			HitsE = lhevent->GetHitsE();
//
//			bscal->Init();
//			bscal->SetEvent(lhevent);
//			bscal->SetCore(corex, corey);
//			bscal->SetBadHit(*HitsE, ned, "ED");
//			bscal->Setnparticle(*HitsE, ned, 40.);
//			int tWind, rWind;
//			int NfiltE;
//			tWind = 400; rWind = 100;
//			bscal->spacetimefilter(*HitsE, ned, tWind, rWind, "ED");
//
//			if (bscal->planarfit(*HitsE, ned, "ED") == 0)
//			{
//				if (bscal->conicalfit_dir(*HitsE, ned, 0.035, "ED") != 0)continue;
//				NfiltE = bscal->noisefilter(*HitsE, ned, -50, 100, rWind + 100, "ED");
//				if (NfiltE > 200)
//				{
//					bscal->noisefilter(*HitsE, ned, -50, 100, rWind + 300, "ED");
//					if (bscal->conicalfit_dir(*HitsE, ned, -1, "ED") != 0)continue;
//				}
//				else
//				{
//					if (bscal->conicalfit_dir(*HitsE, ned, 0.035, "ED") != 0)continue;
//				}
//				the_temp = bscal->GetTheta(); phi_temp = bscal->GetPhi();
//				if (the_temp > -0.5)
//				{
//					h2e_argo(mjd, the_temp * DR2D, phi_temp * DR2D, &ra_temp, &dec_temp);
//					precession2Old(mjd, &ra_temp, &dec_temp);
//					h2_z->Fill(ra_temp, dec_temp);
//					//printf("Ra_z = %lf, Dec_z = %lf\n", ra_temp, dec_temp);
//				}
//			}
//
//			bscal->SetToffset("./config/TrueCal_mean_211001_221001.txt");
//			if (bscal->planarfit(*HitsE, ned, "ED") == 0)
//			{
//				bscal->conicalfit_dir(*HitsE, ned, 0.035, "ED");
//				NfiltE = bscal->noisefilter(*HitsE, ned, -50, 100, rWind + 100, "ED");
//				if (NfiltE > 200)
//				{
//					bscal->noisefilter(*HitsE, ned, -50, 100, rWind + 300, "ED");
//					if (bscal->conicalfit_dir(*HitsE, ned, -1, "ED") != 0)continue;
//				}
//				else
//				{
//					if (bscal->conicalfit_dir(*HitsE, ned, 0.035, "ED") != 0)continue;
//				}
//				the_temp = bscal->GetTheta(); phi_temp = bscal->GetPhi();
//				if (the_temp > -0.5)
//				{
//					h2e_argo(mjd, the_temp * DR2D, phi_temp * DR2D, &ra_temp, &dec_temp);
//					precession2Old(mjd, &ra_temp, &dec_temp);
//					h2_m->Fill(ra_temp, dec_temp);
//					//printf("Ra_m = %lf, Dec_m = %lf\n", ra_temp, dec_temp);  
//				}
//			}
//
//			n_event++;
//
//		}
//		printf("NTotalEntry = %d\n", n_event);
//		infile->Close();
//
//	}
//
//	histfile->Write(0, TObject::kOverwrite);
//	histfile->Close();
//	printf("The num of events for calibration : %d\n", n_event);
//	delete(lhevent);
//}

int Calibration(BSTimeoffCalibration* BS, double rin, double rout, double El, double dCx = 0, double dCy = 0, int Nfit = 0)
{
	int i, j, NTotalEntry, n_event = 0;
	double mjd, ra, dec, E, corex, corey, npe3, num4;
	LHEvent* lhevent = new LHEvent();
	char filehead[200] = "root://eos01.ihep.ac.cn//eos/user/q/qijincan/brightsource/result/LHE";
	char file[200];
	char suffix[20] = ".root";
	const char* filedate;

	for (i = 2021202; i < 2022250; i++)
	{
		if (i % 1000 == 366)
		{
			i = i + 1000 - 366;
			continue;
		}

		filedate = std::to_string(i).c_str();
		strcpy(file, filehead);
		strcat(file, filedate);
		strcat(file, suffix);

		TFile* infile = TFile::Open(file);
		if (!infile) {
			std::cout << "Can't open file:" << file << std::endl;
			continue;
		}
		if (infile->IsZombie())
		{
			std::cout << file << " ";
			printf("file error!!\n");
			infile->Close();
			continue;
		}

		TTree* tree = (TTree*)infile->Get("Event");
		if (!tree)
		{
			printf("get event error!!\n");
			infile->Close();
			continue;
		}

		TNtupleD* nt = (TNtupleD*)infile->Get("nt1");
		if (!nt)
		{
			printf("get nt1 error!!\n");
			infile->Close();
			continue;
		}

		tree->SetBranchAddress("event", &lhevent);
		nt->SetBranchAddress("r_Corex", &corex);
		nt->SetBranchAddress("r_Corey", &corey);
		nt->SetBranchAddress("NpE3", &npe3);
		nt->SetBranchAddress("NuM4", &num4);
		nt->SetBranchAddress("Ra", &ra);
		nt->SetBranchAddress("Dec", &dec);
		nt->SetBranchAddress("Energy", &E);
		nt->SetBranchAddress("mjd", &mjd);

		NTotalEntry = tree->GetEntriesFast();

		for (j = 0; j < NTotalEntry; j++)
		{
			nt->GetEntry(j);;
			//if (ra<Source_Ra - DRa || ra>Source_Ra + DRa || dec<Source_Dec - DDec || dec>Source_Dec + DDec)continue;
			double dr = pow(ra - (Source_Ra + dCx), 2) + pow(dec - (Source_Dec + dCy), 2);
			if (dr<rin * rin || dr>rout * rout)continue;
			if (pow(10., E) < El)continue;
			if ((log10((num4 + 0.0001) / npe3) > -2.3))continue;
			tree->GetEntry(j);

			//BS->SetSource(mjd, Source_Ra + dCx, Source_Dec + dCy);
			BS->SetSource(mjd, ra, dec);
			BS->SetEvent(lhevent);
			BS->SetCore(corex, corey);
			if (BS->eventcallineED() < 0)continue;
			n_event++;

			if (n_event == Nfit)break;

		}
		printf("NTotalEntry = %d\n", n_event);
		infile->Close();

		if (n_event == Nfit)break;

	}

	printf("The num of events for calibration : %d\n", n_event);
	delete(lhevent);
	return n_event;
}

/*
int main(int argc, char* argv[])
{
	int i, j, ITER = 20;
	double* count, * toff;
	double rin, rout, El, dCenterx, dCentery;
	rin = atof(argv[3]);
	rout = atof(argv[4]);
	El = atof(argv[5]);
	dCenterx = atof(argv[6]);
	dCentery = atof(argv[7]);
	char outpath[200] = "root://eos01.ihep.ac.cn//eos/user/q/qijincan/brightsource/result/";
	char temp[200];
	BSTimeoffCalibration* bscal = BSTimeoffCalibration::GetInstance();

	strcat(outpath, argv[1]);
	TFile* outfile = TFile::Open(outpath, "recreate");
	TNtupleD* nt_out;

	ITER = atoi(argv[2]);
	for (i = 0; i < ITER; i++)
	{
		bscal->SetITER(i);
		bscal->Init();
		if (i == 0)bscal->SetToffset("./config/TrueCal_mean_211001_221001.txt");
		Calibration(bscal, rin, rout, El, dCenterx, dCentery);
		bscal->GenerateToffset("ED");
		bscal->NormToffset();
		printf("ITER:%d\n", i + 1);

		const char* i_str = std::to_string(i + 1).c_str();
		strcpy(temp, "EDs_ITER_");
		strcat(temp, i_str);

		count = bscal->GetstatED();
		toff = bscal->GetToffsetED();
		nt_out = new TNtupleD(temp, "timeoff for EDs of KM2A", "ID:cal_count:timeoffset");
		nt_out->SetDirectory(outfile);
		for (j = 0; j < MNED; j++)
		{
			if (count[j] == 0)continue;
			nt_out->Fill((double)j, count[j], toff[j]);
		}
		outfile->Write(0, TObject::kOverwrite);
		delete(nt_out);

	}

	outfile->Close();

	printf("outpath:%s ,temp:%s\n", outpath, temp);

	//WritetxtFile(outpath, temp);
	return 0;
}
*/