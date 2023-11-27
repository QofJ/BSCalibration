#include "BSTimeoffCalibration.h"
#include "TNtupleD.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include <stdio.h>
#include <stdlib.h>    
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>

#define DRa 1.
#define DDec 1.

int Calibration(BSTimeoffCalibration* BS, double rin, double rout, double El, double dCx = 0, double dCy = 0)
{
	int i, j, NTotalEntry, n_event = 0;
	double mjd, ra, dec, E, corex, corey, npe3, num4, c0, alpha;
	LHEvent* lhevent = new LHEvent();
	char filehead[200] = "root://eos01.ihep.ac.cn//eos/user/q/qijincan/brightsource/result/LHE";
	char file[200];
	char suffix[20] = ".root";
	const char* filedate;

	//for (i = 2021202; i < 2021250; i++)
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
		nt->SetBranchAddress("a", &alpha);
		nt->SetBranchAddress("r_C0", &c0);

		NTotalEntry = tree->GetEntriesFast();

		for (j = 0; j < NTotalEntry; j++)
		{
			nt->GetEntry(j);
			TClonesArray* HitsE;
			int ned = lhevent->GetNhitE();
			HitsE = lhevent->GetHitsE();
			
			//if (ra<Source_Ra - DRa || ra>Source_Ra + DRa || dec<Source_Dec - DDec || dec>Source_Dec + DDec)continue;
			double dr = pow(ra - (Source_Ra + dCx), 2) + pow(dec - (Source_Dec + dCy), 2);
			if (dr<rin * rin || dr>rout * rout)continue;
			if (pow(10., E) < El)continue;
			if ((log10((num4 + 0.0001) / npe3) > -2.3))continue;
			tree->GetEntry(j);

			BS->SetSource(mjd, ra, dec);
			BS->SetSource(mjd, ra, dec);
			BS->SetEvent(lhevent);
			BS->SetCore(corex, corey);
			BS->SetA(alpha);
			BS->SetC0(c0);
		//	BS->SetC0(c0 + gRandom->Gaus(0,1.79));
			BS->SetBadHit(*HitsE, ned, "ED");
			BS->Setnparticle(*HitsE, ned, 40.);
			BS->noisefilter(*HitsE, ned, -50, 100, 100 + 100, "ED");
			BS->Addresidual(*HitsE, ned, "ED");
			//if (BS->eventcallineED() < 0)continue;
			n_event++;

		}
		printf("NTotalEntry = %d\n", n_event);
		infile->Close();

	}

	printf("The num of events for calibration : %d\n", n_event);
	delete(lhevent);
	return n_event;
}

int WritetxtFile(const char* outroot, const char* outnt)
{
	int i, ID, NED;
	double id, toff;
	TFile* outfile = TFile::Open(outroot);
	TNtupleD* nt_out = (TNtupleD*)outfile->Get(outnt);
	nt_out->SetBranchAddress("ID", &id);
	nt_out->SetBranchAddress("timeoffset", &toff);

	char txtpath[200] = "/eos/user/q/qijincan/brightsource/result/calEDs.txt";
	NED = nt_out->GetEntriesFast();
	FILE* f = fopen(txtpath, "w");

	fprintf(f, "Type ID Timeoffset\n");
	for (i = 0; i < NED; i++)
	{
		nt_out->GetEntry(i);
		ID = (int)id;
		fprintf(f, "ED %d %lf\n", ID, toff);
	}
	fclose(f);
	outfile->Close();

	return 0;

}


//argv[1]:outfile argv[2]:ITER
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
		if(i==0)bscal->SetToffset("./config/TrueCal_mean_211001_221001.txt");
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
