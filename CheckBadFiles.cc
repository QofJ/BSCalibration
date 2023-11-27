//#include "BSTimeoffCalibration.h"
#include "TNtupleD.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TFile.h"
#include <stdio.h>
#include <stdlib.h>    
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>

#define DRa 1.
#define DDec 1.

void CheckBadFiles()
{
	int i, j, NTotalEntry, n_event = 0;
	double mjd, ra, dec, E, corex, corey, npe3, num4;
//	LHEvent* lhevent = new LHEvent();
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
		infile->Close();
	}
}
