#include <stdio.h>
#include <stdlib.h>    
#include <math.h>
#include <string.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
#include <iostream>
#include "TFile.h"
#include "TNtupleD.h"


int Fitpol3()
{
	char filelist[200] = "/eos/user/q/qijincan/brightsource/result/filelistforFit.txt";
	char file[200];
	const int Nrow = 9;
	int ID_C[5500] = { 0 };
	int i, j, k;
	int ID;
	double toff;
	double temp;
	double Toff[5500][Nrow + 1] = { 0 };
	double r[Nrow] = { 1,1.4,1.7,2.0,2.2,2.4,2.6,2.8,3 };
	double sup[3] = { 0,2,3 };
	TMatrixD A(3, 3);
	TVectorD aTb(3), result(3);
	aTb.Zero();
	result.Zero();
	
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			temp = 0;
			for (k = 0; k < Nrow; k++)
			{
				temp += pow(r[k], double(sup[i] + sup[j]));
			}
			A(i, j) = temp;
		}
	}

	i = 0;
	FILE* flist = fopen(filelist, "r");
	while (!feof(flist))
	{
		fscanf(flist, "%s\n", file);
		FILE* f = fopen(file, "r");
		if (!f)
		{
			printf("can't open file:%s\n", file);
		}

		printf("%s\n", file);

		fscanf(f, "%*[^\n]%*c");
		while (!feof(f))
		{
			fscanf(f, "ED %d %lf\n", &ID, &toff);
			if (ID < 0 || ID >= 5500)continue;
			ID_C[ID] += 1;
			Toff[ID][i] = toff;
		}
		fclose(f);
		i++;
	}
	fclose(flist);

	TDecompSVD svd(A);
	for (i = 0; i < 5500; i++)
	{
		printf("Fit ID:%d\n", i);

		if (ID_C[i] < Nrow)continue;
		for (j = 0; j < Nrow; j++)
		{
			aTb(0) += Toff[i][j];
			aTb(1) += pow(r[j], 2) * Toff[i][j];
			aTb(2) += pow(r[j], 3) * Toff[i][j];
		}
		Bool_t ok;
		result = svd.Solve(aTb, ok);
		if (!ok)
		{
			ID_C[i] = -1;
			continue;
		}
		Toff[i][Nrow] = result(0);
		aTb.Zero();
	}

	FILE* fout = fopen("/eos/user/q/qijincan/brightsource/result/calEDs_pol3.txt", "w");
	fprintf(fout, "r: 1 1.4 1.7 2.0 2.2 2.4 2.6 2.8 3.0\n");
	for (i = 0; i < 5500; i++)
	{
		printf("i = %d\n", i);

		if (ID_C[i] < Nrow)continue;
		fprintf(fout, "ED %d %lf\n", i, Toff[i][Nrow]);
	}
	fclose(fout);
	return 0;
}

int WritetxtFile(const char* outroot, const char* outnt, const char* outtxt = "/eos/user/q/qijincan/brightsource/result/calEDs.txt")
{
	int i, ID, NED;
	double id, toff;
	TFile* outfile = TFile::Open(outroot);
	TNtupleD* nt_out = (TNtupleD*)outfile->Get(outnt);
	nt_out->SetBranchAddress("ID", &id);
	nt_out->SetBranchAddress("timeoffset", &toff);
	//nt_out->SetBranchAddress("Toffset", &toff);

	//char txtpath[200] = "/eos/user/q/qijincan/brightsource/result/calEDs.txt";
	char txtpath[200];
	strcpy(txtpath, outtxt);
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

//char outroot[200] = "root://eos01.ihep.ac.cn//eos/user/q/qijincan/brightsource/result/Cal_E10_1p5_2.root";
//char outnt[20] = "EDs_ITER_30";
//WritetxtFile(outroot, outnt);



int Normlization(const char* calfile)
{
	char outfile[200];
	strcpy(outfile, calfile);
	strcat(outfile, "_norm");
	int id[5500], n(0);
	double time[5500], mean(0);
	FILE* f1 = fopen(calfile, "r");
	FILE* f2 = fopen(outfile, "w");
	
	fscanf(f1, "%*[^\n]%*c");
	while (!feof(f1))
	{
		fscanf(f1, "%*s %d %lf\n", &id[n], &time[n]);
		mean += time[n];
		n++;
	}
	if (n != 0)mean /= n;

	fprintf(f2, "ED ID NormlizedTimeoffset\n");
	for (int i = 0; i < n; i++)
	{
		fprintf(f2, "ED %d %lf\n", id[i], time[i] - mean);
	}

	fclose(f1);
	fclose(f2);

	return 0;
}

//char calfile[200] = "/eos/user/q/qijincan/brightsource/result/calEDs.txt";
//Normlization(calfile);

TGraph* DrawIter(const char* rootfile, int ITER, int index, int Norm = 0)
{
	int i;
	TFile* f = TFile::Open(rootfile);

	double id, toff;
	char Leaf[30];
	char Prefix[20] = "EDs_ITER_";

	TGraph* g = new TGraph;
	g->SetMarkerStyle(8);

	for (i = 0; i < ITER; i++)
	{
		strcpy(Leaf, Prefix);
		strcat(Leaf, std::to_string(i + 1).c_str());

		TNtupleD* nt = (TNtupleD*)f->Get(Leaf);
		nt->SetBranchAddress("ID", &id);
		nt->SetBranchAddress("timeoffset", &toff);

		double mean = 0;
		if (Norm == 1)
		{
			int NEntry = nt->GetEntriesFast();
			for (int j = 0; j < NEntry; j++)
			{
				nt->GetEntry(j);
				mean += toff;
			}
			mean /= NEntry;
		}

		nt->GetEntry(index);

		g->SetPoint(i, (double)(i + 1), toff - mean);

	}
	printf("ID = %d\n", (int)id);
	g->Draw("ACP");
	return g;
}


//char outroot[200] = "root://eos01.ihep.ac.cn//eos/user/q/qijincan/brightsource/result/CalVf_Daz_E20_0_1.root";
//auto g = DrawIter(outroot, 20, 0);

TGraph* DrawIterStd(const char* rootfile, int ITER)
{
	int i;
	TFile* f = TFile::Open(rootfile);

	double id, toff;
	char Leaf[30];
	char Prefix[20] = "EDs_ITER_";

	TGraph* g = new TGraph;
	g->SetMarkerStyle(8);

	for (i = 0; i < ITER; i++)
	{
		strcpy(Leaf, Prefix);
		strcat(Leaf, std::to_string(i + 1).c_str());

		TNtupleD* nt = (TNtupleD*)f->Get(Leaf);
		nt->Draw("timeoffset>>hist_name(500,-50,50)", "", "goff");
		TH1D* hist = (TH1D*)gDirectory->Get("hist_name");
		double std = hist->GetStdDev();

		g->SetPoint(i, (double)(i + 1), std);

	}
	printf("ID = %d\n", (int)id);
	g->GetXaxis()->SetTitle("Iterations");
	g->GetYaxis()->SetTitle("Time offset StdDev");
	g->Draw("ACP");
	return g;
}

TGraph* AddIterStd(const char* rootfile, int ITER, TGraph* g)
{
	int i;
	TFile* f = TFile::Open(rootfile);

	double id, toff, num = g->GetN();
	char Leaf[30];
	char Prefix[20] = "EDs_ITER_";

	g->SetMarkerStyle(8);

	for (i = 0; i < ITER; i++)
	{
		strcpy(Leaf, Prefix);
		strcat(Leaf, std::to_string(i + 1).c_str());

		TNtupleD* nt = (TNtupleD*)f->Get(Leaf);
		nt->Draw("timeoffset>>hist_name(500,-50,50)", "", "goff");
		TH1D* hist = (TH1D*)gDirectory->Get("hist_name");
		double std = hist->GetStdDev();

		g->SetPoint(num + i, (double)(num + i + 1), std);

	}
	printf("ID = %d\n", (int)id);
	//g->Draw("ACP");
	return g;
}

TGraph* DrawIterDelta(const char* rootfile, int ITER, string mode = "std")
{
	int i;
	TFile* f = TFile::Open(rootfile);

	double id1, toff1, id2, toff2;
	double dif_ave = 0;
	char Leaf[30];
	char Prefix[20] = "EDs_ITER_";

	TGraph* g = new TGraph;
	g->SetMarkerStyle(8);

	for (i = 0; i < ITER - 1; i++)
	{
		dif_ave = 0;
		strcpy(Leaf, Prefix);
		//strcat(Leaf, std::to_string(i + 1).c_str());
		char temp[30];
		strcpy(temp, Leaf);
		strcat(temp, std::to_string(i + 1).c_str());

		TNtupleD* nt1 = (TNtupleD*)f->Get(temp);
		nt1->SetBranchAddress("ID", &id1);
		nt1->SetBranchAddress("timeoffset", &toff1);

		strcpy(temp, Leaf);
		strcat(temp, std::to_string(i + 2).c_str());
		TNtupleD* nt2 = (TNtupleD*)f->Get(temp);
		nt2->SetBranchAddress("ID", &id2);
		nt2->SetBranchAddress("timeoffset", &toff2);

		TH1D h("h1", "", 1000, -25, 25);


		int NEntry = nt1->GetEntriesFast();
		double nt1_ave = 0, nt2_ave = 0;
		for (int j = 0; j < NEntry; j++)
		{
			nt1->GetEntry(j);
			nt2->GetEntry(j);
			nt1_ave += toff1;
			nt2_ave += toff2;
		}	
		nt1_ave /= NEntry;
		nt2_ave /= NEntry;
		for (int j = 0; j < NEntry; j++)
		{
			nt1->GetEntry(j);
			nt2->GetEntry(j);
			h.Fill(toff1 - toff2);
			dif_ave += fabs(toff2 - toff1 - nt2_ave + nt1_ave);

		}	
		dif_ave /=NEntry;

        if(mode == "std")
		g->SetPoint(i, (double)(i + 1), h.GetStdDev());
		else
		g->SetPoint(i, (double)(i + 1), dif_ave);

	}
	g->Draw("ACP");
	return g;
}

int Scatter(const char* calfile1, const char* calfile2)
{
	int n1(0), n2(0);
	double id1[5500], id2[5500];
	double toff1[5500], toff2[5500];
	double mean1 = 0, mean2 = 0;
	char temp[20];
	FILE* f1 = fopen(calfile1, "r");
	FILE* f2 = fopen(calfile2, "r");

	fscanf(f1, "%*[^\n]%*c");
	fscanf(f2, "%*[^\n]%*c");

	while (!feof(f1))
	{
		fscanf(f1, "%s %lf %lf", temp, id1 + n1, toff1 + n1);
		if (!(strcmp(temp, "ED") == 0))break;
		mean1 += toff1[n1];
		n1++;
	}
	if (n1 != 0)mean1 /= n1;
	while (!feof(f2))
	{
		fscanf(f2, "%s %lf %lf", temp, &id2[n2], &toff2[n2]);
		if (!(strcmp(temp, "ED") == 0))break;
		mean2 += toff2[n2];
		n2++;
	}
	if (n2 != 0)mean2 /= n2;
	for (int i = 0; i < n1 || i < n2; i++)
	{
		if (i < n1)toff1[i] -= mean1;
		if (i < n2)toff2[i] -= mean2;
	}

	TGraph* g1 = new TGraph(n1, id1, toff1);
	TGraph* g2 = new TGraph(n2, id2, toff2);
	g1->SetMarkerStyle(34);
	g2->SetMarkerStyle(26);
	g1->Draw("AP");
	g2->Draw("Psame");
	fclose(f1);
	fclose(f2);
	
	return 0;
}

//char calfile1[200] = "/lhaasofs/user/lhaasorec/cal/km2a/TDC/TDCCalPublish/cal20211212-gt20-2.txt";
//char calfile2[200] = "/eos/user/q/qijincan/brightsource/result/calEDs1.txt";
//Scatter(calfile1, calfile2);

TH1D* Hist(const char* calfile1, const char* calfile2)
{
	int n1(0), n2(0);
	int id;
	double time, dtime;
	int num = 0;
	double Exy = 0;
	double id1[5500] = { 0 }, id2[5500] = { 0 };
	double toff1[5500], toff2[5500];
	double mean1 = 0, mean2 = 0;
	char temp[20];

	TH1D* h1 = new TH1D("h1", "cal minus true cal", 100, -10., 10.);
	h1->SetTitle(";Timeoffset_BS method - Timeoffset_CP method (ns);Number of EDs");

	FILE* f1 = fopen(calfile1, "r");
	FILE* f2 = fopen(calfile2, "r");

	fscanf(f1, "%*[^\n]%*c");
	fscanf(f2, "%*[^\n]%*c");

	while (!feof(f1))
	{
		fscanf(f1, "%s %d %lf", temp, &id, &time);
		if (!(strcmp(temp, "ED") == 0))break;
		if (id >= 0 && id < 5500)
		{
			id1[id] = 1;
			toff1[id] = time;
		}
		mean1 += time;
		n1++;
	}
	if (n1 != 0)mean1 /= n1;
	while (!feof(f2))
	{
		fscanf(f2, "%s %d %lf", temp, &id, &time);
		if (!(strcmp(temp, "ED") == 0))break;
		if (id >= 0 && id < 5500)
		{
			id2[id] = 1;
			toff2[id] = time;
		}
		mean2 += time;
		n2++;
	}
	if (n2 != 0)mean2 /= n2;
	double dif_ave = 0;
	for (int i = 0; i < 5500; i++)
	{
		if (id1[i] == 1)toff1[i] -= mean1;
		if (id2[i] == 1)toff2[i] -= mean2;
		if (id1[i] == 1 && id2[i] == 1)
		{
			dtime = toff1[i] - toff2[i];
			h1->Fill(dtime);
			Exy += toff1[i] * toff2[i];
			dif_ave += fabs(dtime);
			num++;
		}
	}

	Exy /= num;
	dif_ave /= num;
	printf("Exy = %lf\n", Exy);
	printf("Average Difference = %lf\n", dif_ave);
	fclose(f1);
	fclose(f2);

	return h1;
}

//char calfile1[200] = "/eos/user/q/qijincan/brightsource/result/calEDs2122E10I30C1p5_2_norm.txt";
//char calfile2[200] = "/lhaasofs/user/lhaasorec/cal/km2a/TDC/TDCCalPublish/cal20211212-gt20-2.txt";
//auto h1 = Hist(calfile1, calfile2);
//h1->Draw();

TH1D* HistDif(const char* calfile1, const char* calfile2)
{
	int n1(0), n2(0);
	int id;
	double time, dtime;
	int num = 0;
	double Exy = 0;
	double id1[5500] = { 0 }, id2[5500] = { 0 };
	double toff1[5500], toff2[5500];
	double mean1 = 0, mean2 = 0;
	char temp[20];

	TH1D* h1 = new TH1D("h1", "cal minus true cal", 5500, 0, 5500);

	FILE* f1 = fopen(calfile1, "r");
	FILE* f2 = fopen(calfile2, "r");

	fscanf(f1, "%*[^\n]%*c");
	fscanf(f2, "%*[^\n]%*c");

	while (!feof(f1))
	{
		fscanf(f1, "%s %d %lf", temp, &id, &time);
		if (!(strcmp(temp, "ED") == 0))break;
		if (id >= 0 && id < 5500)
		{
			id1[id] = 1;
			toff1[id] = time;
		}
		mean1 += time;
		n1++;
	}
	if (n1 != 0)mean1 /= n1;
	while (!feof(f2))
	{
		fscanf(f2, "%s %d %lf", temp, &id, &time);
		if (!(strcmp(temp, "ED") == 0))break;
		if (id >= 0 && id < 5500)
		{
			id2[id] = 1;
			toff2[id] = time;
		}
		mean2 += time;
		n2++;
	}
	if (n2 != 0)mean2 /= n2;
	for (int i = 0; i < 5500; i++)
	{
		if (id1[i] == 1)toff1[i] -= mean1;
		if (id2[i] == 1)toff2[i] -= mean2;
		if (id1[i] == 1 && id2[i] == 1)
		{
			dtime = toff1[i] - toff2[i];
			h1->Fill(i + 0.5, dtime * 100);
			Exy += toff1[i] * toff2[i];
			num++;
		}
	}

	Exy /= num;
	printf("Exy = %lf\n", Exy);

	fclose(f1);
	fclose(f2);

	return h1;
}

//char calfile1[200] = "/eos/user/q/qijincan/brightsource/result/test.txt";
//char calfile2[200] = "/eos/user/q/qijincan/brightsource/result/TrueCal_mean_211001_221001.txt";
//auto h1 = HistDif(calfile1, calfile2);
//h1->Draw();

void Filter(TH1D* h1)
{
	int nbin = h1->GetNbinsX();
	for (int i = 1; i <= nbin; i++)
	{
		double content = h1->GetBinContent(i);
		if (content < 2)
		{
			h1->SetBinContent(i, 0);
		}
	}
}



TH1D* Hist1(const char* calfile)
{
	int n1(0);
	int id;
	double time;
	double id1[5500] = { 0 };
	double toff1[5500];
	double mean1 = 0;
	char temp[20];

	TH1D* h1 = new TH1D("h1", "BS cal", 100, -10., 10.);

	FILE* f = fopen(calfile, "r");

	fscanf(f, "%*[^\n]%*c");

	while (!feof(f))
	{
		fscanf(f, "%s %d %lf", temp, &id, &time);
		if (!(strcmp(temp, "ED") == 0))break;
		if (id >= 0 && id < 5500)
		{
			id1[id] = 1;
			toff1[id] = time;
		}
		mean1 += time;
		n1++;
	}
	if (n1 != 0)mean1 /= n1;
	printf("mean = %lf\n", mean1);
	
	for (int i = 0; i < 5500; i++)
	{
		if (id1[i] == 1)
		{
			toff1[i] -= mean1;
			h1->Fill(toff1[i]);
		}
	}
	fclose(f);

	return h1;
}

//char calfile[200] = "/lhaasofs/user/lhaasorec/cal/km2a/TDC/TDCCalPublish/cal20211212-gt20-2.txt";
//auto h2 = Hist1(calfile);
//h2->Draw();

int DrawCompare(const char* calfile1, const char* calfile2)
{
	int id, n1(0), n2(0), id1[5500] = { 0 }, id2[5500] = { 0 };
	double time, toff1[5500], toff2[5500];
	char temp[20];

	FILE* f1 = fopen(calfile1, "r");
	FILE* f2 = fopen(calfile2, "r");;

	fscanf(f1, "%*[^\n]%*c");
	while (!feof(f1))
	{
		fscanf(f1, "%s %d %lf", temp, &id, &time);
		if (strcmp(temp, "ED") != 0)break;
		if (id >= 0 && id < 5500)
		{
			id1[id] = 1;
			toff1[id] = time;
		}
	}

	fscanf(f2, "%*[^\n]%*c");
	while (!feof(f2))
	{
		fscanf(f2, "%s %d %lf", temp, &id, &time);
		if (!(strcmp(temp, "ED") == 0))break;
		if (id >= 0 && id < 5500)
		{
			id2[id] = 1;
			toff2[id] = time;
		}
	}

	fclose(f1);
	fclose(f2);

	int i, j(0);
	TGraph* gp = new TGraph();
	TGraph* gl = new TGraph();
	gp->SetTitle(";Time offset_CP method (ns);Time offset_BS method (ns)");
	gp->SetMarkerStyle(8);
	gl->SetPoint(0, -5, -5);
	gl->SetPoint(1, 12, 12);
	gl->SetLineWidth(3);
	gl->SetLineColor(kBlue);
	for (i = 0; i < 5500; i++)
	{
		if (id1[i] != 1)continue;
		if (id2[i] != 1)continue;
		gp->SetPoint(j++, toff1[i], toff2[i]);
	}

	gp->Draw("AP");
	gl->Draw("Csame");
	return 0;
	
}

//char calfile1[200] = "/eos/user/q/qijincan/brightsource/result/calEDs_2122_E10_I50_1p5_norm.txt";
//char calfile2[200] = "/lhaasofs/user/lhaasorec/cal/km2a/TDC/TDCCalPublish/cal20211212-gt20-2.txt";
//DrawCompare(calfile2, calfile1);

void HistDDegDC0(TH1D& h1, TH1D& h2, TH1D& h3, int n, const char* f1, const char* f2)
{
	TFile* froot1 = TFile::Open(f1);
	TFile* froot2 = TFile::Open(f2);
	TNtupleD* nt1 = (TNtupleD*)froot1->Get("nt1");
	TNtupleD* nt2 = (TNtupleD*)froot2->Get("nt1");
	double the_1, the_2, phi_1, phi_2, c0_1, c0_2, ddeg, dc0;
	double cx_1, cx_2, cy_1, cy_2, dcore;
	int i;

	nt1->SetBranchAddress("r_Theta", &the_1);
	nt1->SetBranchAddress("r_Phi", &phi_1);
	nt1->SetBranchAddress("r_C0", &c0_1);
	nt1->SetBranchAddress("r_Corex", &cx_1);
	nt1->SetBranchAddress("r_Corey", &cy_1);
	nt2->SetBranchAddress("r_Theta", &the_2);
	nt2->SetBranchAddress("r_Phi", &phi_2);
	nt2->SetBranchAddress("r_C0", &c0_2);
	nt2->SetBranchAddress("r_Corex", &cx_2);
	nt2->SetBranchAddress("r_Corey", &cy_2);
	
	for (i = 0; i < n; i++)
	{
		nt1->GetEntry(i);
		nt2->GetEntry(i);
		ddeg = sqrt(pow(57.29578 * (the_1 - the_2), 2) + pow(57.29578 * (phi_1 - phi_2), 2));
		if (ddeg > 2)printf("index(DDeg>2) = %d\n", i);
		dc0 = c0_1 - c0_2;
		dcore = sqrt(pow(cx_1 - cx_2, 2) + pow(cy_1 - cy_2, 2));
		h1.Fill(ddeg);
		h2.Fill(dc0);
		h3.Fill(dcore);
	}
}

//int n = 1000;
//char f1[200] = "/home/lhaaso/qijincan/KM2ADatarec_V2/rectest20210819.root";
// // char f2[200] = "/eos/lhaaso/rec/km2a/full_array/2021_V2/0819/ES.118545.KM2A_EVENT.TEST_5_TELESCOPE_EDMD_OVERLAP.es-11.20210819235924.284.dat.root";
//char f2[200] = "/home/lhaaso/qijincan/KM2ADatarec_V2/rectest20210819_onlyTDC.root";
//TH1D h1("h1", "DDeg", 1000, -25, 25);
//TH1D h2("h2", "DC0", 1000, -25, 25);
//TH1D h3("h3", "DCore", 1000, -25, 25);
//HistDDegDC0(h1, h2, h3, n, f1, f2);

TH1D* ShowEDCal(int ID, int Iter, TFile* fin)
{
	const char* iter_str = std::to_string(Iter).c_str();
	char temp[20];
	TH1D* h1;

	strcpy(temp, "h2_");
	strcat(temp, iter_str);

	h1 = ((TH2D*)fin->Get(temp))->ProjectionX("h1", ID, ID);
	return h1;
}

//auto fin = TFile::Open("root://eos01.ihep.ac.cn//eos/user/q/qijincan/brightsource/result/Calrec_2_2122_E20_0_1.root");
//int Iter = 1;
//int ID = 1;
//TH1D* h1 = ShowEDCal(ID, Iter, fin);
//h1->Draw();

double GetMPValue(TH2D* h2, int id, double delta = 10)
{
	TH1D* h1 = h2->ProjectionX("h1", id, id);
	double mean = h1->GetMean();
	double temp = 0;
	int ite = 30;
	for (int i = 0; i < ite; i++)
	{
		temp = mean;
		h1->SetAxisRange(mean - delta, mean + delta);
		mean = h1->GetMean();
		if (fabs(mean - temp) < 0.01)break;
	}
	return mean;
}

//auto f1 = TFile::Open("root://eos01.ihep.ac.cn//eos/user/q/qijincan/brightsource/result/Calrec_MPV_2122_E20_0_1.root");
//auto h2 = (TH2D*)f1->Get("h2_1");
//double MPV;
//MPV = GetM	PValue(h2, 1591);

TH1D* FindEntries50(TH2D* h2)
{
	TH1D* h3 = new TH1D("h1_e", "h1", 500, 0, 500);
	for (int i = 1; i < 5500; i++)
	{
		TH1D* h1 = h2->ProjectionX("h1", i, i);
		if (h1->GetEntries() != 0)h3->Fill(double(h1->GetEntries()));
		if (h1->GetEntries() < 50 && h1->GetEntries() != 0)
		{
			printf("ID = %d, N_Entries = %lf\n", i, h1->GetEntries());
		}
	}

	h3->Draw();

	return h3;
}

TH1D* TrueCalPlusGaus(double sigma, const char* calfile1 = "/eos/user/q/qijincan/brightsource/result/TrueCal_mean_211001_221001.txt")
{
	int n1(0);
	int id;
	double time, dtime;
	double id1[5500] = { 0 };
	double toff1[5500], toff2[5500];
	double mean1 = 0, mean2 = 0;
	char temp[20];

	TH1D* h1 = new TH1D("h1", "cal minus true cal", 100, -10., 10.);

	FILE* f1 = fopen(calfile1, "r");

	fscanf(f1, "%*[^\n]%*c");

	while (!feof(f1))
	{
		fscanf(f1, "%s %d %lf", temp, &id, &time);
		if (!(strcmp(temp, "ED") == 0))break;
		if (id >= 0 && id < 5500)
		{
			id1[id] = 1;
			toff1[id] = time;
			toff2[id] = time + 2 * gRandom->Gaus(0, sigma);
		}
		mean1 += time;
		mean2 += toff2[id];
		n1++;
	}
	if (n1 != 0)
	{
		mean1 /= n1;
		mean2 /= n1;
	}
	
	for (int i = 0; i < 5500; i++)
	{
		if (id1[i] == 1)
		{
			toff1[i] -= mean1;
			toff2[i] -= mean2;
		}
		
		if (id1[i] == 1)
		{
			dtime = toff1[i] - toff2[i];
			h1->Fill(dtime);
		}
	}
	fclose(f1);

	return h1;
}

//double sigma = 0.8;
//auto h1 = TrueCalPlusGaus(sigma);
//h1->Draw();

//读取文件中某一探测器的timeoffset
void readTimeOffset(const char* filename, int ITER, const int n) {
    // 打开文件
	char fileName[200] = "root://eos01.ihep.ac.cn//eos/user/q/qijincan/brightsource/result/";
	char ntupleName[200] = "EDs_ITER_";
	strcat(fileName,filename);
	strcat(ntupleName,std::to_string(ITER).c_str());
	const char* leafName = "timeoffset";
    TFile* file = TFile::Open(fileName);
    if (!file) {
        std::cout << "Failed to open file: " << fileName << std::endl;
        return;
    }
    // 获取指定TNtupleD
    TNtupleD* ntuple = dynamic_cast<TNtupleD*>(file->Get(ntupleName));
    if (!ntuple) {
        std::cout << "Failed to get TNtupleD: " << ntupleName << std::endl;
        file->Close();
        return;
    }
    // 获取指定TLeafD
    TLeafD* leaf = dynamic_cast<TLeafD*>(ntuple->GetLeaf(leafName));
    if (!leaf) {
        std::cout << "Failed to get TLeafD: " << leafName << std::endl;
        file->Close();
        return;
    }
    // 读取指定元素
    double value = 0;
    if (ntuple->GetEntries() > n) {
		ntuple->GetEntry(n);
        value = leaf->GetValue(0);
        std::cout << "The " << n << "th element of " << leafName << " is: " << value << std::endl;
    }
    else {
        std::cout << "The " << n << "th element of " << leafName << " doesn't exist!" << std::endl;
    }
    // 关闭文件
    file->Close();
}




//char outroot[200] = "root://eos01.ihep.ac.cn//eos/user/q/qijincan/brightsource/result/CalVf_MPV_PreToff_confit0p0001_E20_0_1.root";
//char outnt[20] = "EDs_ITER_20";
//WritetxtFile(outroot, outnt, "/eos/user/q/qijincan/brightsource/result/test.txt");

//char calfile[200] = "/eos/user/q/qijincan/brightsource/result/calEDs.txt";
//Normlization(calfile);

//char calfile1[200] = "/eos/user/q/qijincan/brightsource/result/TrueCal_mean_211001_221001.txt";
//char calfile2[200] = "/eos/user/q/qijincan/brightsource/result/calEDs1.txt";
//Scatter(calfile1, calfile2);

//char calfile1[200] = "/eos/user/q/qijincan/brightsource/result/calEDs2122E10I30C0_0p25_norm.txt";
//char calfile2[200] = "/eos/user/q/qijincan/brightsource/result/TrueCal_mean_211001_221001.txt";
//auto h1 = Hist(calfile1, calfile2);
//h1->Draw();

//char calfile[200] = "/eos/user/q/qijincan/brightsource/result/TrueCal_mean_211001_221001.txt";
//auto h2 = Hist1(calfile);
//h2->Draw();

//char outroot[200] = "root://eos01.ihep.ac.cn//eos/user/q/qijincan/brightsource/result/CalVf_MPV_PreToffGaus_confit0p0001_E10_0_1_Ra113.root";
//char outnt[20] = "EDs_ITER_6";
//WritetxtFile(outroot, outnt, "/eos/user/q/qijincan/brightsource/result/test.txt");
//char calfile1[200] = "/eos/user/q/qijincan/brightsource/result/test.txt";
////char calfile2[200] = "/eos/user/q/qijincan/brightsource/result/TrueCal_mean_211001_221001_norm.txt";
//char calfile2[200] = "/home/lhaaso/qijincan/KM2ADatarec_V2/config/pre_TimeOffset_Gaus.txt";
////DrawCompare(calfile2, calfile1);
//auto h1 = Hist(calfile1, calfile2);h1->Draw()
//auto h2 = Hist1(calfile1);

//TGraph* g = new TGraph;
//AddIterStd("/eos/user/q/qijincan/brightsource/result/CalVf_MPV_PreToff_confit0p0001_E20_0_1.root", 20, g);
//AddIterStd("/eos/user/q/qijincan/brightsource/result/CalVf_MPV_PreToff_confit0p0001_E20_0_1_FromI20.root", 20, g);
//AddIterStd("/eos/user/q/qijincan/brightsource/result/CalVf_MPV_PreToff_confit0p0001_E20_0_1_FromI40.root", 20, g);
//AddIterStd("/eos/user/q/qijincan/brightsource/result/CalVf_MPV_PreToff_confit0p0001_E20_0_1_FromI60.root", 30, g);
//AddIterStd("/eos/user/q/qijincan/brightsource/result/CalVf_MPV_PreToff_confit0p0001_E20_0_1_FromI90.root", 30, g);
//AddIterStd("/eos/user/q/qijincan/brightsource/result/CalVf_MPV_PreToff_confit0p0001_E20_0_1_FromI120.root", 30, g);
//AddIterStd("/eos/user/q/qijincan/brightsource/result/CalVf_MPV_PreToff_confit0p0001_E20_0_1_FromI150.root", 30, g);
//AddIterStd("/eos/user/q/qijincan/brightsource/result/CalVf_MPV_PreToff_confit0p0001_E20_0_1_FromI180.root", 30, g);
//AddIterStd("/eos/user/q/qijincan/brightsource/result/CalVf_MPV_PreToff_confit0p0001_E20_0_1_FromI210.root", 30, g);
//AddIterStd("/eos/user/q/qijincan/brightsource/result/CalVf_MPV_PreToff_confit0p0001_E20_0_1_FromI240.root", 30, g);
//AddIterStd("/eos/user/q/qijincan/brightsource/result/CalVf_MPV_PreToff_confit0p0001_E20_0_1_FromI270.root", 30, g);
//g->Draw("ACP")