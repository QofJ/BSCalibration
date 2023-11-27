#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "astro-csz.h"
#include "rece.h"
#include "KM2AEvent.h"
#include "LHEvent.h"
#include "G4KM2A_Geometry.h"
#include "LHCalibration.h"

//int KM2Arecdec(int argc,char* argv[])
int main(int argc, char* argv[])
{
	int i, j, nflag;
	double Theta, Phi, Corex, Corey, MJD, T0, Ra, Dec, a, age, size, E, NuM4, NpE3;
	Double_t* ntptr;
	Double_t out_array[30];
	int out_index[] = { 0,2,3,4,5,6,7,9,10,11,13,20,21,22,25,26 };
	int lenindex = sizeof(out_index) / sizeof(out_index[0]);
	int NTotalEntry;
	KM2AEvent* e = new KM2AEvent();
	LHEvent* ge = new LHEvent();

	G4KM2A_Geometry* geom = G4KM2A_Geometry::GetInstance(7);
	LHCalibration* cal = LHCalibration::GetInstance();
	FILE* filelist;
	char headname[100] = "root://eos01.ihep.ac.cn/";
	char decodepath[100] = "/eos/lhaaso/decode/km2a/2021/";
	char filepath[200];
	char decodefile[200];
	char outhead[100] = "root://eos01.ihep.ac.cn//eos/user/q/qijincan/brightsource/result/LHE_BK";
	char outpath[200];
	char stemp[200];

	filelist = fopen(argv[1], "r");

	strcpy(outpath, outhead);
	strcat(outpath, argv[2]);
	TFile* outfile = TFile::Open(outpath, "RECREATE");
	TTree* t_out = new TTree("Event", "Reduced Data for Calibration");
	TNtupleD* nt_out = new TNtupleD("nt1", "KM2A Reconstructed Data", "Nevent:mjd:NhitE:NfiltE:NhitM:NfiltM:NpE1:NpE2:NpE3:"
		"NuM1:NuM2:NuM3:NuM4:NuM5:r_Corex:r_Corey:r_Theta:r_Phi:r_C0:a:chi:NtrigE:NtrigE2:size:age:dr:nMD50:Ra:Dec:Energy");

	t_out->Branch("event", "LHEvent", &ge);
	//t_out->SetAutoSave(100);

	while (!feof(filelist))
	{
		fscanf(filelist, "%s\n", stemp);
		printf("The present file:%s\n", stemp);
		strcpy(filepath, headname);
		strcat(filepath, stemp);

		strcpy(decodefile, headname);
		strcat(decodefile, decodepath);
		strcat(decodefile, strstr(stemp, "/ES") - 4);
		printf("The decode file:%s\n", decodefile);

		//strcpy(outpath, outhead);
		//strcat(outpath, strchr(stemp, '-') + 1);
		TFile* infile = TFile::Open(filepath);
		if (!infile)continue;
		if (infile->IsZombie() || infile->GetEND() < 10000) {
			printf("file error!!\n");
			infile->Close();
			continue;
		}
		TNtupleD* nt = (TNtupleD*)infile->Get("nt1");
		if (!nt) {
			printf("get nt error!!\n");
			infile->Close();
			continue;
		}

		TFile* decfile = TFile::Open(decodefile);
		if (!decfile)continue;
		if (decfile->IsZombie() || decfile->GetEND() < 10000)
		{
			printf("decfile error!!\n");
			infile->Close();
			decfile->Close();
			continue;
		}

		ntptr = nt->GetArgs();
		
		nt->SetBranchAddress("r_Theta", &Theta);
		nt->SetBranchAddress("r_Phi", &Phi);
		nt->SetBranchAddress("r_Corex", &Corex);
		nt->SetBranchAddress("r_Corey", &Corey);
		nt->SetBranchAddress("mjd", &MJD);
		nt->SetBranchAddress("r_C0", &T0);
		nt->SetBranchAddress("NpE3", &NpE3);
		nt->SetBranchAddress("NuM4", &NuM4);
		nt->SetBranchAddress("a", &a);
		nt->SetBranchAddress("age", &age);
		nt->SetBranchAddress("size", &size);
		

		TTree* eve = (TTree*)decfile->Get("event");
		if (!eve)
		{
			printf("get event error!!\n");
			decfile->Close();
			continue;
		}
		eve->SetBranchAddress("Event", &e);


		NTotalEntry = nt->GetEntriesFast();

		for (i = 0; i < NTotalEntry; ++i)
		{
			nt->GetEntry(i);
			//h2e_argo(MJD, Theta * DR2D, Phi * DR2D, &Ra, &Dec);
			//precession2Old(MJD, &Ra, &Dec);
			h2e_argo(MJD, Theta * DR2D, Phi * DR2D, &Ra, &Dec);
			precession2Old(MJD, &Ra, &Dec);
			if ((Ra < 100 || Dec < 22.0145-2 || Dec > 22.0145+2))continue;

			E = recer50new3(age, size, Theta);

			if ((log10((NuM4 + 0.0001) / NpE3) <= -2.3) && (E != -1))
			{

				eve->GetEntry(i);
				nflag = cal->ApplyKM2AEvent(e, ge);
				ge->SetEvN(e->EvN());
				ge->SetMjd(e->Mjd());
				ge->SetDt(e->Dt());

				if (nflag < 2)continue;
				for (j = 0; j < lenindex; j++)
				{
					out_array[out_index[j]] = ntptr[out_index[j]];
				}
				out_array[1] = MJD;
				out_array[8] = NpE3;
				out_array[12] = NuM4;
				out_array[14] = Corex;
				out_array[15] = Corey;
				out_array[16] = Theta;
				out_array[17] = Phi;
				out_array[18] = T0;
				out_array[19] = a;
				out_array[23] = size;
				out_array[24] = age;
				out_array[27] = Ra;
				out_array[28] = Dec;
				out_array[29] = E;

				/*
				ge->SetNhit(0);
				ge->SetNmd(0);
				ge->SetTheta(Theta);
				ge->SetPhi(Phi);
				ge->SetCorex(Corex);
				ge->SetCorey(Corey);
				ge->SetMjd(MJD);
				ge->SetC0(T0);
				ge->SetA(a);
				ge->SetRa(Ra);
				ge->SetDec(Dec);
				ge->SetE(E);
				ge->SetNpE3(NpE3);
				ge->SetNuM4(NuM4);
				
				//nt_out->Fill(Theta, Phi, Corex, Corey, MJD, T0, a, Ra, Dec, E);

				for (j = 0; j < (e->NHit() + e->Nmd()); ++j)
				{

					ge->AddHit(((KM2AHit*)(*e->Hitlist())[j])->ID(), ((KM2AHit*)(*e->Hitlist())[j])->Mode(), ((KM2AHit*)(*e->Hitlist())[j])->Time());
				}*/


				t_out->Fill();
				nt_out->Fill(out_array);
				ge->Initcsz();
			}
		}
		//nt_out->SetDirectory(outfile);
		//outfile->Write();
		//delete(nt_out);nt_out = NULL;
		//outfile->Close();
		infile->Close();
		decfile->Close();
	}
	outfile->Write(0, TObject::kOverwrite);
	outfile->Close();
	fclose(filelist);
	return 0;
}
