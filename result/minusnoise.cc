#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int minusnoise()
{
	int n0 = 75381, n1(74577), n2(101778);
	double r0(1), r1(1.5), r2(2);
	double N_noise, ratio;
	double a[2];
	int ID[5500] = { 0 };
	double Toff[5500];
	int id0, id1, id2;
	double t0, t1, t2;
	double toff_noise, temp1, temp2;

	N_noise = (n1 + n2) * pow(r0, 2) / (pow(r2, 2) - pow(r0, 2));
	ratio = (n0 - N_noise) / N_noise;
	a[0] = (pow(r1, 3) - pow(r0, 3)) / (pow(r1, 2) - pow(r0, 2));
	a[1] = (pow(r2, 3) - pow(r1, 3)) / (pow(r2, 2) - pow(r1, 2));
	FILE* f0 = fopen("/home/lhaaso/qijincan/KM2ADatarec_V2/result/calEDs2122E10I30C0_1_norm.txt", "r");
	FILE* f1 = fopen("/home/lhaaso/qijincan/KM2ADatarec_V2/result/calEDs2122E10I30C1_1p5_norm.txt", "r");
	FILE* f2 = fopen("/home/lhaaso/qijincan/KM2ADatarec_V2/result/calEDs2122E10I30C1p5_2_norm.txt", "r");
	FILE* fout = fopen("/home/lhaaso/qijincan/KM2ADatarec_V2/result/calEDs_minusnoise.txt", "w");
	fscanf(f0, "%*[^\n]%*c");
	fscanf(f1, "%*[^\n]%*c");
	fscanf(f2, "%*[^\n]%*c");

	while (!feof(f0) && !feof(f1) && !feof(f2))
	{
		fscanf(f0, "%*s %d %lf", &id0, &t0);
		fscanf(f1, "%*s %d %lf", &id1, &t1);
		fscanf(f2, "%*s %d %lf", &id2, &t2);
		if (id0 == id1 && id0 == id2)
		{
			ID[id0] = 1;
			temp2 = (t2 - t1) / (a[1] - a[0]);
			temp1 = t1 - temp2 * a[0];
			toff_noise = temp1 + temp2 * r0;
			Toff[id0] = (t0 - toff_noise + t0 * ratio) / ratio;
		}
	}

	fprintf(fout, "TYPE ID Toffset_noisereduction\n");
	for (int i = 0; i < 5500; i++)
	{
		if (ID[i] != 1)continue;
		fprintf(fout, "ED %d %lf\n", i, Toff[i]);
	}

	fclose(f0);
	fclose(f1);
	fclose(f2);
	fclose(fout);
	return 0;
}