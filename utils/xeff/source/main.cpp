#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <assert.h>
#include <xeff_table.h>
#include <vector>

using namespace std;

static int isTurbulent = 1;
static double M;

/**
 * @brief Выдать X эффективное в точках боковой поверхности ЛА.
 * @details Текстовый файл содержит параметры ЛА в следующем формате:
 * До конца файла, строки вида: <Length> <P/P0'> <rm>
 * Здесь Length - длина вдоль образующей
 * P/P0' - коэффициент давления.
 * rm - радиус миделя.
 * Файл не должен заканчиваться пустой строкой.
 * @param filename - имя файла, из которого считываются параметры ЛА и газодинамика.
 */
void calculateXeff(char* filename, bool HasVelocityData) {
	int PRESSURE_POINTS_NUM;
	vector<pcell_t*> pressure;
	char outFilename[256];
	FILE* fout;
	sprintf(outFilename, "r%s", filename);
	fout = fopen(outFilename, "wt");
	assert(fout != 0);
	
	fstream fp(filename, ios::in);
	/* Считать параметры для xeff */
	cout<<"M: "<<M<<"\n";
	printf("Source data:\n");
	printf("LENGTH\tPdivP0\tMIDSHIP_RADIUS");
	if (HasVelocityData)
		printf("\tVdivVh\n");
	else
		printf("\n");
	int j=0;
	while (!fp.eof())
	{
		pressure.push_back(new pcell_t);
		fp >> pressure[j]->length;
		if (!fp)
			break;
		fp >> pressure[j]->PdivP0;
		fp >> pressure[j]->rm;
		if (HasVelocityData)
			fp >> pressure[j]->VdivVh;
		printf("%lf\t%lf\t%lf", pressure[j]->length, pressure[j]->PdivP0, pressure[j]->rm);
		if (HasVelocityData)
			printf("\t%lf\n", pressure[j]->VdivVh);
		if (pressure[j]->length < 0.) {
			printf("[EE]: Value of Length Along Gen is'nt Correct at line %d!\n", j+1);
			exit(-1);
		}
		if ((pressure[j]->PdivP0 <= 0.) || (pressure[j]->PdivP0 > 1.2)) {
			printf("[EE]: Value of PdivP0 is'nt Correct at line %d!\n", j+1);
			exit(-1);
		}
		if (pressure[j]->rm < 0.) {
			printf("[EE]: Value of Midship Radius is'nt Correct at line %d!\n", j+1);
			exit(-1);
		}
		assert(pressure[j]->VdivVh > 0);
		j++;
	}
	printf("\nSource data loaded\n");
	PRESSURE_POINTS_NUM = j;
	
	CXeff_table xeff(M, pressure, PRESSURE_POINTS_NUM, HasVelocityData);
	fprintf(fout,"%-10.10s\t%-8.8s\t%-8.8s\n", "length", "Xeff^0.2", "Xeff^0.5");
	for (j=0; j<PRESSURE_POINTS_NUM; j++) {
		double xel = pow(xeff.calculate(pressure[j]->length, 0), 0.5);
		double xet = pow(xeff.calculate(pressure[j]->length, 1), 0.2);
		fprintf(fout,"%10.4lf\t%8.5lf\t%8.5lf\n", pressure[j]->length, xet, xel);
	}
	
	fclose(fout);
	fp.close();
	printf("\nXeff Calculation finished.\n");
	return;
}

int main(int argc, char* argv[])
{
	bool HasVelocityData = false;

	M = atoi(argv[2]);
	if (strcmp(argv[3], "-V") == 0)
		HasVelocityData = true;
	
	calculateXeff(argv[1], HasVelocityData);
	
	return 0;
}
