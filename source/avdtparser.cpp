#include <avdtparser.h>
#include <material.h>

/**
 * Прочитать данные с устройства
 * @param f - указатель на устройство ввода
 * @param fmt - формат данных.
 * @param par - указатель на память, в которую надо записать данные
 * @param isNewStr - признак новой строки.
 * @return Результат операции
 */
#include <string>
using namespace std;
static char str[STRING_MAX_LEN];
string format;
CUserMaterial* m[MATERIALS_MAX_NUM];
static int read(FILE* f, const char* fmt, void* par, bool isNewStr = false)
{
	if (isNewStr) {
		format = "";
		char* res;
		if (feof(f))
			return -1;
		res = fgets(str, STRING_MAX_LEN, f);
		if (res == 0)
			return -1;
	}
	string tmp = format + "%";
	tmp += fmt;
//	printf("SSS: %s\n", tmp.c_str());
	int r = sscanf(str, tmp.c_str(), par);
	format +="%*";
	format += fmt;
	return r;
}
trm_t* trm_parse(const char* filename, FILE* device)
{
	trm_t* trm = new trm_t;
	FILE* ftr = fopen(filename, "rt");
	if (ftr == 0) {
		printf("FILE NOT FOUND: %s\n", filename);
		system("pause");
		exit(-1);
	}
	/* --- Прочитать файл ИД с траекторией --- */
	read(ftr, "d", &(trm->POINTS_NUM), true); fprintf(device, "%d\t", trm->POINTS_NUM);
	if (trm->POINTS_NUM <= 0) {
		printf("[EE]: Incorrect number of trajectory points: %d\n", trm->POINTS_NUM);
		exit(-1);
	}
	if (trm->POINTS_NUM >= TRAJECTORY_MAX_LEN) {
		printf("[EE]: Number of trajectory points [%d] exceed the up limit [%d]!\n", trm->POINTS_NUM, TRAJECTORY_MAX_LEN);
		exit(-1);
	}
	read(ftr, "lf", &(trm->BEGIN_TIME)); fprintf(device, "%8.2lf\t", trm->BEGIN_TIME);
	read(ftr, "lf", &(trm->END_TIME)); fprintf(device, "%8.2lf\t", trm->END_TIME);
	if (trm->BEGIN_TIME >= trm->END_TIME) {
		printf("[EE]: BEGIN calculation time [%lf] more than END time [%lf]!\n", trm->BEGIN_TIME, trm->END_TIME);
		exit(-1);
	}
	read(ftr, "lf", &(trm->THETA)); fprintf(device, "%8.2lf\n", trm->THETA);
	if ((trm->THETA < -90.) || (trm->THETA > 90.)) {
		printf("[EE]: THETA value [%lf] is out of range!\n", trm->THETA);
		exit(-1);
	}
	double time[TRAJECTORY_MAX_LEN]; /* Массив опорных моментов времени. */
	/* Прочитать опорные моменты времени. */
	for (int i=0; i<trm->POINTS_NUM;) {
		bool r = true;
		for (int j=0; ((i<trm->POINTS_NUM) && (j < 6)); i++, j++) {
			read(ftr, "lf", &(time[i]), r);
			r = false;
			fprintf(device, "%11.5lf\t", time[i]);
		}
		fprintf(device, "\n");
	}
	/* Прочитать опорные скорости. */
	fprintf(device, "--- VELOCITY ---\n");
	for (int i=0; i<trm->POINTS_NUM; ) {
		double V;
		bool r = true;
		for (int j=0; ((i<trm->POINTS_NUM) && (j < 6)); i++, j++) {
			read(ftr, "lf", &(V), r);
			r = false;
			if (V < 0.) {
				printf("[EE]: VELOCITY value [%lf] is out of range\n", V);
				exit(-1);
			}
			fprintf(device, "%11.5lf\t", V);
			trm->V.add(time[i], V);
		}
		fprintf(device, "\n");
	}
//	trm->V.print();
	/* Прочитать опорные высоты. */
	fprintf(device, "--- HEIGHT ---\n");
	for (int i=0; i<trm->POINTS_NUM; ) {
		double H;
		bool r = true;
		for (int j=0; ((i<trm->POINTS_NUM) && (j < 6)); i++, j++) {
			read(ftr, "lf", &(H), r);
			r = false;
			if (H < 0.) {
				printf("[EE]: HEIGHT value [%lf] is out of range\n", H);
				exit(-1);
			}
			fprintf(device, "%11.5lf\t", H);
			trm->H.add(time[i], H);
		}
		fprintf(device, "\n");
	}
	/* Прочитать опорные углы атаки. */
	fprintf(device, "--- ANGLE OF ATTACK ---\n");
	for (int i=0; i<trm->POINTS_NUM; ) {
		double AL;
		bool r = true;
		for (int j=0; ((i<trm->POINTS_NUM) && (j < 6)); i++, j++) {
			read(ftr, "lf", &(AL), r);
			r = false;
			if ((AL < -75.) || (AL > 75.)) {
				printf("[EE]: ANGLE of ATTACK value [%lf] is out of range\n", AL);
				exit(-1);
			}
			fprintf(device, "%11.5lf\t", AL);
			trm->AL.add(time[i], AL);
		}
		fprintf(device, "\n");
	}
	/* Прочитать опорные углы поворота вокруг оси. */
	fprintf(device, "--- PHI ---\n");
	for (int i=0; i<trm->POINTS_NUM; ) {
		double PHI;
		bool r = true;
		for (int j=0; ((i<trm->POINTS_NUM) && (j < 6)); i++, j++) {
			read(ftr, "lf", &(PHI), r);
			r = false;
			fprintf(device, "%11.5lf\t", PHI);
			trm->PHI.add(time[i], PHI);
		}
		fprintf(device, "\n");
	}
	fclose(ftr);
	return trm;
}
thm_t* thm_parse(const char* filename, FILE* device, CBluntedCone** BCone, gasdynamics_t* gd,  print_t* prn)
{
	assert(filename != 0);
	assert(device != 0);
	assert(BCone != 0);
	assert(gd != 0);
	assert(prn != 0);

	thm_t* thm = new thm_t;
	FILE* ftps = fopen(filename, "rt");
	if (ftps == 0) {
		printf("FILE NOT FOUND: %s\n", filename);
		system("pause");
		exit(-1);
	}
	/* --- Прочитать файл ИД с пакетами материалов --- */
	int LAYERS; /* Количество слоев материалов в расчетной области. */
	double T0; /* Начальная температура. */
	func_points_t T0_FUNC;
	read(ftps, "d", &(LAYERS), true); fprintf(device, "LAYERS=%d\t", LAYERS);
	assert((LAYERS > 0) && (LAYERS < LAYERS_MAX_NUM));
	read(ftps, "lf", &(T0)); fprintf(device, "T0=%8.2lf\t", T0);
	assert(T0 > 0.);
	read(ftps, "d", &(prn->PRINT_CELLS_NUM)); fprintf(device, "PRINT_CELLS_NUM=%d\t", prn->PRINT_CELLS_NUM);
	read(ftps, "d", &(T0_FUNC.count)); fprintf(device, "T0_POINTS=%d\n", T0_FUNC.count);
	/* --- Интерполяционная функция по температуре. --- */
	if (T0_FUNC.count > 0) {
		T0_FUNC.x = new double[T0_FUNC.count];
		T0_FUNC.y = new double[T0_FUNC.count];
		assert(T0_FUNC.x != 0); assert(T0_FUNC.y != 0);
		bool r = true;
		fprintf(device, "X=\t");
		for (int i=0; i<T0_FUNC.count; i++) {
			read(ftps, "lf", &(T0_FUNC.x[i]), r); fprintf(device, "%5.3E ", T0_FUNC.x[i]);
			r = false;
		}
		fprintf(device, "\nT=\t");
		r = true;
		for (int i=0; i<T0_FUNC.count; i++) {
			read(ftps, "lf", &(T0_FUNC.y[i]), r); fprintf(device, "%10.4lf ", T0_FUNC.y[i]);
			r = false;
		}
		fprintf(device, "\n");
	}
	fflush(NULL);
	int isComplexTFH[LAYERS_MAX_NUM]; /* Признак табличного задания ТФХ материала на i-том слое. */
	bool r = true;
	for (int i=0; i<LAYERS; i++) {
		read(ftps, "d", &(isComplexTFH[i]), r);
		r = false;
		fprintf(device, "%d\t", isComplexTFH[i]);
	}
	fprintf(device, "\n");
	fflush(NULL);
	int AblationType[LAYERS_MAX_NUM]; /* Тип уноса для материала на i-том слое. */
	r = true;
	for (int i=0; i<LAYERS; i++) {
		read(ftps, "d", &(AblationType[i]), r);
		r = false;
		fprintf(device, "%d\t", AblationType[i]);
	}
	fprintf(device, "\n");
	fflush(NULL);
	 /* --- Сканирование теплофизики материалов. --- */
	double LAYER_DX[LAYERS_MAX_NUM]; /* Толщина ячейки в слое материала, м. */
	double LAYER_CP[LAYERS_MAX_NUM]; /* Теплоёмкость материала, Дж/кг*К. */
	double LAYER_D[LAYERS_MAX_NUM]; /* Плотность материала, кг/м3 */
	double LAYER_L[LAYERS_MAX_NUM]; /* Теплопроводность материала, Вт/м*К. */
	double LAYER_A[LAYERS_MAX_NUM]; /* Коэффициент А в уравнении уноса. */
	double LAYER_B[LAYERS_MAX_NUM]; /* Коэффициент В в уравнении уноса. */
	double LAYER_TU[LAYERS_MAX_NUM]; /* Температура уноса, К. */
	double LAYER_EPS[LAYERS_MAX_NUM]; /* Степень внешней поверхности материала. */
	double LAYER_AT[LAYERS_MAX_NUM]; /* Ламинарный коэффициент А. */
	fprintf(device, "%10.10s\t%8.8s\t%8.8s\t%8.8s\t%8.8s\t%8.8s\t%8.8s\t%8.8s\t%8.8s\n", "DX", "CP", "DENS", "L", "A", "B", "Tdestr", "EPS", "AT");
	for (int i=0; i<LAYERS; i++) {
		read(ftps, "lf", &(LAYER_DX[i]), true); fprintf(device, "%5.3E\t", LAYER_DX[i]);
		read(ftps, "lf", &(LAYER_CP[i])); fprintf(device, "%8.2lf\t", LAYER_CP[i]);
		read(ftps, "lf", &(LAYER_D[i])); fprintf(device, "%8.2lf\t", LAYER_D[i]);
		read(ftps, "lf", &(LAYER_L[i])); fprintf(device, "%8.2lf\t", LAYER_L[i]);
		read(ftps, "lf", &(LAYER_A[i])); fprintf(device, "%8.2lf\t", LAYER_A[i]);
		read(ftps, "lf", &(LAYER_B[i])); fprintf(device, "%8.2lf\t", LAYER_B[i]);
		read(ftps, "lf", &(LAYER_TU[i])); fprintf(device, "%8.2lf\t", LAYER_TU[i]);
		read(ftps, "lf", &(LAYER_EPS[i])); fprintf(device, "%8.2lf\t", LAYER_EPS[i]);
		read(ftps, "lf", &(LAYER_AT[i])); fprintf(device, "%8.2lf\n", LAYER_AT[i]);
		if (isComplexTFH[i] == 1)
			m[i] = new CUserMaterial("USERMAT", stdout, LAYER_L[i], LAYER_D[i], LAYER_CP[i], LAYER_EPS[i], LAYER_TU[i], LAYER_A[i], LAYER_B[i], AblationType[i]);
		else
			m[i] = new CUserMaterial("CONSTMAT", stdout, LAYER_L[i], LAYER_D[i], LAYER_CP[i], LAYER_EPS[i], LAYER_TU[i], LAYER_A[i], LAYER_B[i], AblationType[i]);
		fflush(NULL);
	}
	int LAYER_CELLS[LAYERS_MAX_NUM];
	/* --- Количество слоев в каждом материале. --- */
	r = true;
	for (int i=0; i<LAYERS; i++) {
		read(ftps, "d", &(LAYER_CELLS[i]), r); fprintf(device, "%d\t", LAYER_CELLS[i]);
		r = false;
	}
	fprintf(device, "\n");
	fflush(NULL);
	r = true;
	/* --- Сканирование номеров ячеек, температуры которых будут выводиться на печать. --- */
	for (int i=0; i<prn->PRINT_CELLS_NUM; i++) {
		read(ftps, "d", &(prn->PRINT_CELLS[i]), r); fprintf(device, "PrintCell%d=%d\n", i, prn->PRINT_CELLS[i]);
		r = false;
	}
	fprintf(device, "\n");
	/* --- Сканирование специальных параметров. --- */
	fprintf(device, "--- SPECIAL PARAMETERS ---\n");
	fprintf(device, "%11.11s\t%11.11s\t%11.11s\t%11.11s\t%11.11s\t%11.11s\t%11.11s\n", "AXIS_LEN", "R0", "TH0", "TIMESTEP", "TIME_PRINT", "TR_HEIGHT", "PHI0");
	/* Расстояние вдоль оси симметрии. */
	read(ftps, "lf", &(gd->X), true); fprintf(device, "%11.5lf\t", gd->X);
	double R0; /* Радиус притупления, м. */
	read(ftps, "lf", &(R0)); fprintf(device, "%11.5lf\t", R0);
	double TH; /* Угол полураствора, град. */
	read(ftps, "lf", &(TH)); fprintf(device, "%11.5lf\t", TH);
	/* Начальный шаг счёта, с. */
	read(ftps, "lf", &(thm->INIT_TIMESTEP)); fprintf(device, "%11.5lf\t", thm->INIT_TIMESTEP);
	/* Шаг печати, с. */
	read(ftps, "lf", &(prn->print_interval)); fprintf(device, "%11.5lf\t", prn->print_interval);
	/* Высота перехода, м. */
	read(ftps, "lf", &(gd->HT)); fprintf(device, "%11.5lf\t", gd->HT);
	/* Начальный угол проворота образующей с исследуемой точкой. */
	read(ftps, "lf", &(gd->PHI0)); fprintf(device, "%11.5lf\n", gd->PHI0);
	*BCone = new CBluntedCone(R0, TH);
	/* --- Сканирование базы интерполяции по газодинамическим параметрам. --- */
	fprintf(device, "--- MACH, ALPHA AND PHI ---\n");
	double MACHS[6];
	r = true;
	for (int i=0; i<6; i++) {
		read(ftps, "lf", &(MACHS[i]), r); fprintf(device, "%11.5lf\t", MACHS[i]);
		r = false;
	}
	fprintf(device, "\n");
	r = true;
	double ALPHAS[7];
	for (int i=0; i<7; i++) {
		read(ftps, "lf", &(ALPHAS[i]), r); fprintf(device, "%11.5lf\t", ALPHAS[i]);
		r = false;
	}
	fprintf(device, "\n");
	r = true;
	double PHIS[3];
	for (int i=0; i<3; i++) {
		read(ftps, "lf", &(PHIS[i]), r); fprintf(device, "%11.5lf\t", PHIS[i]);
		r = false;
	}
	fprintf(device, "\n");
	/* --- Сканирование газодинамических параметров. --- */
	/* -- Коэффициент давления P/P0 */
	fprintf(device, "--- PP0 ---\n");
//	double PP0[3][7][6];
	for (int k=0; k<3; k++) {
		ITable * tbl = new ITable();
		for (int i=0; i<7; i++) {
			IFunc* func = new IFunc();
			double pp0;
			r = true;
			for (int j=0; j<6; j++) {
				read(ftps, "lf", &(pp0), r); fprintf(device, "%11.5lf\t", pp0);
				fflush(NULL);
				func->add(MACHS[j], pp0);
				r = false;
			}
			tbl->add(func, ALPHAS[i]);
			fprintf(device, "\n");
		}
		gd->PP0.add(tbl, PHIS[k]);
	}
	/* -- Эффективная длина для турбулентного режима Xeft^0.2 */
	fprintf(device, "--- XEFT ---\n");
	for (int k=0; k<3; k++) {
		ITable * tbl = new ITable();
		
		for (int i=0; i<7; i++) {
			IFunc* func = new IFunc();
			double xet;
			r = true;
			for (int j=0; j<6; j++) {
				read(ftps, "lf", &(xet), r); fprintf(device, "%11.5lf\t", xet);
				
				r = false;
				func->add(MACHS[j], xet);
			}
			tbl->add(func, ALPHAS[i]);
			fprintf(device, "\n");
		}
		gd->XET.add(tbl, PHIS[k]);
	}
	fflush(NULL);
	/* -- Эффективная длина для ламинарного режима Xefl^0.5 */
	fprintf(device, "--- XEFL ---\n");
	for (int k=0; k<3; k++) {
		ITable * tbl = new ITable();
		for (int i=0; i<7; i++) {
			IFunc* func = new IFunc();
			double xel;
			r = true;
			for (int j=0; j<6; j++) {
				read(ftps, "lf", &(xel), r); fprintf(device, "%11.5lf\t", xel);
				r = false;
				func->add(MACHS[j], xel);
			}
			tbl->add(func, ALPHAS[i]);
			fprintf(device, "\n");
		}
		gd->XEL.add(tbl, PHIS[k]);
	}
	fflush(NULL);
	/* --- Сканирование переменной теплофизики. --- */
	func_points_t tfh;
	tfh.count = 10;
	tfh.x = new double[10];
	tfh.y = new double[10];
	for (int i=0; i<LAYERS; i++) {
		if (isComplexTFH[i] == 1) {
			r = true;
			for (int j=0; j<10; j++) {
				read(ftps, "lf", &(tfh.x[j]), r); fprintf(device, "%11.5lf\t", tfh.x[j]);
				r = false;
			}
			fprintf(device, "\n");
			fflush(NULL);
			r = true;
			for (int j=0; j<10; j++) {
				read(ftps, "lf", &(tfh.y[j]), r); fprintf(device, "%11.5lf\t", tfh.y[j]);
				r = false;
			}
			m[i]->update("CP", &tfh);
			fprintf(device, "\n");
			fflush(NULL);
			r = true;
			for (int j=0; j<10; j++) {
				read(ftps, "lf", &(tfh.y[j]), r); fprintf(device, "%11.5lf\t", tfh.y[j]);
				r = false;
			}
			m[i]->update("L", &tfh);
			fprintf(device, "\n");
			fflush(NULL);
		}
	}
	fclose(ftps);
	fprintf(device, "\n%s\t%s\t%s\n", "LAYER_DX[i]", "m[i]", "T0");
	for (int i=0; i<LAYERS; i++) { /* Создать расчётную модель изделия. */
		thm->add(LAYER_CELLS[i], LAYER_DX[i], m[i], T0);
		fprintf(device, "%lf\t%s\t%lf\n", LAYER_DX[i], m[i]->name(), T0);
	}
	/* Установка начального поля температур, переменного по времени. */
	if (T0_FUNC.count > 0) {
		thm->setTemperature(&T0_FUNC);
		fprintf(device, "SET TEMPERATURE PROFILE:\n");
		double x = 0.;
		fprintf(device, "%14.14s\t%14.14s\n", "COORDINATE", "TEMPERATURE");
		for (int i=thm->fcnum; i<=thm->lcnum; i++) {
			x += thm->width[i]/2.;
			fprintf(device, "%14.8lf\t%14.4lf\n", x, thm->T[i]);
			x += thm->width[i]/2.;
		}
	}
	fflush(NULL);
	return thm;
}
