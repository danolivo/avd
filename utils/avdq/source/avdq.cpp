#include <stdio.h>
#include <stdlib.h>
#include <avdheatflux.h>

int main(int argc, char* argv[]) {
	char str[256];
	if (argc != 2) {
		printf("USAGE: %s <source filename>\n", argv[0]);
		exit(-1);
	}
	
	sprintf(str, "r%s", argv[1]);
	FILE *fin = fopen(argv[1], "rt");
	if (fin == 0) {
		printf("[EE]: File %s can't be opened!\n", argv[1]);
		exit(-1);
	}
	
	FILE *fout = fopen(str, "wt");
	if (fout == 0) {
		printf("[EE]: File %s can't be created!\n", str);
		exit(-1);
	}
	
	fprintf(fout, "%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\t%11.11s\t%11.11s\t%6.6s\n", "qw", "Ie", "Iw", "alc", "I0", "Cappa", "P0", "P1", "Vstar");
	double angle, Td;
	double Height, Vel, len, Xeff, PdivP0, p0, Temp;
	int isTurbulent;
	
	printf("%5.5s\t%5.5s\t%8.8s\t%6.6s\t%5.5s\t%6.6s\t%6.6s\t%6.6s\t%4s\n", "ANGLE", "Tdest", "HEIGHT", "VEL", "XgenT", "XEFF", "PdivP0", "Tw", "TURB");
	while (!feof(fin)) {
		int res = fscanf(fin, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%1d\n", &angle, &Td, &Height, &Vel, &len, &Xeff, &PdivP0, &Temp, &isTurbulent);
		printf("%5.2lf\t", angle);
		printf("%6.1lf\t", Td);
		printf("%8.1lf\t", Height);
		printf("%6.1lf\t", Vel);
		printf("%5.2lf\t", len);
		printf("%6.4lf\t",  Xeff);
		printf("%6.4lf\t", PdivP0);
		printf("%6.1lf\t", Temp);
		printf("%1d\n", isTurbulent);
		if (res == EOF) {
			printf("[EE]: Parameters at source file is incorrect. Check File!\n");
			exit(-1);
		}
		AVD_Init(Height, Vel, Xeff, PdivP0, Temp, Td, isTurbulent, len, angle);
		fprintf(fout, "%6.1lf\t%6.1lf\t%6.1lf\t%6.4lf\t%6.1lf\t%6.1lf\t%4.2lf\t%6.4E\n", AVD_Qw(), AVD_Ie(), AVD_Iw(), AVD_acp(), AVD_I0(), AVD_Cappa(), AVD_P0(), AVD_P1());
	}
	system("pause");
	return 0;
}
