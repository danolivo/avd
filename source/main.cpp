#include <common.h>
#include <avdtparser.h>
#include <bluntedcone.h>
//#include <avdheatflux.h>
//#include <atmosphere.h>
#include <thsolver.h>
//#include <unos.h>
//#include <avdbc.h>
#include <boundary.h>
#include <avdsolver.h>

int main(int argc, char *argv[])
{
	char iTRFilename[FILENAME_MAX_LEN];
	char iTPSFilename[FILENAME_MAX_LEN];
	char rFilename[FILENAME_MAX_LEN];

	switch (argc) {
	case 0:
		assert(0);
	case 1: /* Ввод данных будет выполнен с консоли. */
		printf("TRAJECTORY FILENAME:");
		scanf("%s", iTRFilename);
		printf("TPS FILENAME:");
		scanf("%s", iTPSFilename);
		break;
	default: /* Нестандартное количество аргументов. */
		printf("INCORRECT PROGRAM USAGE!\n");
		exit(-1);
	}
	sprintf(rFilename, "%s-%s.res", iTRFilename, iTPSFilename);
	FILE* fout = fopen(rFilename, "wt");
	assert(fout != 0);
	/* Разбор файлов ИД. */
	fprintf(fout, "\n--- SOURCES ---\n");
	trm_t* trm = trm_parse(iTRFilename, fout);
	CBluntedCone *BCone;
	gasdynamics_t gd;
	print_t prn;
	thm_t* thm = thm_parse(iTPSFilename, fout, &BCone, &gd, &prn);
//	CAVDBoundary *avdbc = new CAVDBoundary(thm, trm, &gd, BCone, 0);
//	thm->setLBC(avdbc);
	CSOBoundary *bc = new CSOBoundary(0., 0., 0., 0., 0.);
	thm->setRBC(bc);
	AVDSolver solver(thm, trm, &gd, BCone);
	fprintf(fout, "\n--- RESULTS ---\n");
	fprintf(fout, "\n%7.7s\t%7.7s\t%7.7s\t%7.7s\t%7.7s\t%7.7s\t%7.7s\t%7.7s\t%7.7s\t%7.7s\t%7.7s\t%7.7s\t%7.7s\t%5.5s\t%7.7s\t", "TIME", "QCONV", "DY", "T1", "T2", "P", "ALC", "IE", "IW", "FI", "ALF", "H", "V", "MACH", "XEF");
	printf("\n%7.7s\t%7.7s\t%7.7s\t%7.7s\t%7.7s\t%7.7s\t%7.7s\t", "TIME", "QCONV", "DY", "T1", "T2", "H", "V");
	for (int i=0; i<prn.PRINT_CELLS_NUM; i++) {
		fprintf(fout, "%5.5s%6.1lf\t", "CELL", thm->T[prn.PRINT_CELLS[i]]);
		printf("%5.5s%6.1lf\t", "CELL", thm->T[prn.PRINT_CELLS[i]]);
	}
	fprintf(fout, "\n");
	printf("\n");
	thm->CURRENT_TIME = trm->BEGIN_TIME;
	double TIMESTEP = thm->INIT_TIMESTEP;
	/* --- Основной цикл расчёта - итерации по шагу печати --- */
	for (double time = trm->BEGIN_TIME+prn.print_interval; (thm->CURRENT_TIME < trm->END_TIME); time += min(prn.print_interval, trm->END_TIME-time)) {
		avdsolver_t info = solver.Solve(time);
		double Qw = info.avd.QCONV;
		fprintf(fout, "%7.2lf\t", info.time); printf("%7.2lf\t", info.time);
		fprintf(fout, "%7.1lf\t", Qw/4186.8); printf("%7.1lf\t", Qw/4186.8);
		fprintf(fout, "%7.3lf\t", thm->LDEL*1000.); printf("%7.3lf\t", thm->LDEL*1000.);
		fprintf(fout, "%7.1lf\t", thm->TWL); printf("%7.1lf\t", thm->TWL);
		fprintf(fout, "%7.1lf\t", thm->TWR); printf("%7.1lf\t", thm->TWR);
		fprintf(fout, "%7.4lf\t", info.avd.P1);
		fprintf(fout, "%7.4lf\t", info.avd.ALC);
		fprintf(fout, "%7.1lf\t", info.avd.IE/4186.8);
		fprintf(fout, "%7.1lf\t", info.avd.IW/4186.8);
		fprintf(fout, "%7.2lf\t", info.phi);
		fprintf(fout, "%7.2lf\t", info.al);
		fprintf(fout, "%7.2lf\t", info.H/1000.); printf("%7.2lf\t", info.H/1000.);
		fprintf(fout, "%7.1lf\t", info.V); printf("%7.1lf\t", info.V);
		fprintf(fout, "%5.2lf\t", info.mach);
		fprintf(fout, "%7.5lf\t", info.XEF);

		for (int i=0; i<prn.PRINT_CELLS_NUM; i++) {
			fprintf(fout, "%7.1lf\t", thm->T[prn.PRINT_CELLS[i]-1]);
			printf("%7.1lf\t", thm->T[prn.PRINT_CELLS[i]-1]);
		}
		fprintf(fout, "\n");
		printf("\n");
		fflush(NULL);
	}
	fclose(fout);
//	solver.print();
	fflush(NULL);
	system("pause");
	return 0;
}
