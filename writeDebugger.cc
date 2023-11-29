#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TNtupleD.h>

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

int main(int argc, char* argv[])
{
    const char* outroot = argv[1];
    const char* outnt = argv[2];
    const char* outtxt = argv[3];
    WritetxtFile(outroot, outnt, outtxt);
    return 0;
}