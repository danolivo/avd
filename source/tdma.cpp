#include <cmath>
#include <common.h>
#include <tdma.h>

/* ----- DATA ----- */

static double* T;
static double* l;
static double* r;
static double* c;
static double* w;
static double* qv;
static double A[CELLS_MAX_NUM];
static double B[CELLS_MAX_NUM];
static double C[CELLS_MAX_NUM];
static double F[CELLS_MAX_NUM];
static double alpha[CELLS_MAX_NUM];
static double betta[CELLS_MAX_NUM];
static int FIRST;
static int LAST;

/* ----- FUNCTIONS ----- */

/**
 * @brief Метод прогонки.
 * @param T1 - указатель на массив температур, в который будет занесен результат расчета.
 * Размерность T1 должна быть не менее last
 */
static void TDMA()
{
	double hi2 = A[LAST]/C[LAST];
	double mu2 = F[LAST]/C[LAST];
	alpha[FIRST] = 0.;
	betta[FIRST] = 0.;
	/* Прямой ход прогонки. */
	for (int i=0; i<LAST; i++) {
		alpha[i+1] = B[i]/(C[i]-alpha[i]*A[i]);
		betta[i+1] = (A[i]*betta[i]+F[i])/(C[i]-alpha[i]*A[i]);
	}
	/* Обратный ход прогонки. */
	T[LAST] = (mu2+hi2*betta[LAST])/(1-alpha[LAST]*hi2);
	for (int i=LAST-1; i>=0; i--)
		T[i] = alpha[i+1]*T[i+1]+betta[i+1];
}

static double Twl()
{
	return T[FIRST];
}
static double Twr()
{
	return T[LAST];
}

void CalculateTDMA(	double* Temp, double* Width, double* VolHeatSrc,
			double* HeatCond, double* Density, double* SpecHeat,
			double Eps[2], double thau, int size)
{
	T = Temp;
	w = Width;
	qv = VolHeatSrc;
	l = HeatCond;
	r = Density;
	c = SpecHeat;
	FIRST = 0;
	LAST = size-1;
	
	assert((Eps[0] >= 0.) && (Eps[0] <= 1.));
	assert((Eps[1] >= 0.) && (Eps[1] <= 1.));
	assert(thau > 0.); assert(size < CELLS_MAX_NUM);
	assert(T != 0); assert(w != 0);
	assert(qv != 0); assert(l != 0);
	assert(r != 0); assert(c != 0);

	for (int i=FIRST; i<=LAST; i++) {
		double CpRho = c[i]*r[i];
		if (i == FIRST) {
			A[FIRST] = 0.;
			qv[FIRST] += 3.*(Eps[0]*5.67/w[FIRST])*pow(T[FIRST]/100., 4.); // С нормализацией
//			qv[FIRST] -= (Eps[0]*5.67/w[FIRST])*pow(T[FIRST]/100., 4.);
		} else {
			double ll = (w[i-1]+w[i])/(w[i-1]/l[i-1]+w[i]/l[i]);
			A[i] = 2.*thau*(ll/CpRho)/(w[i]*(w[i-1]+w[i]));
		}
		if (i == LAST)
			B[i] = 0.;
		else {
			double lr = (w[i+1]+w[i])/(w[i+1]/l[i+1]+w[i]/l[i]);
			B[i] = 2.*thau*(lr/CpRho)/(w[i]*(w[i+1]+w[i]));
		}
		C[i] = 1.+A[i]+B[i];
		if (i == FIRST)
			C[i] += 4.*Eps[0]*5.67E-02*(thau/CpRho)*pow(T[FIRST]/100., 3.)/w[FIRST]; // С нормализацией
		F[i] = T[i] + thau*qv[i]/CpRho;
		if (i == LAST-1) {
			double eps = Eps[1];
			double A4 = 4*eps*5.67E-08*(thau/CpRho)*w[LAST]*pow(Twr(), 3.)/(w[LAST-1]*l[LAST]*(w[LAST-1]/l[LAST-1]+w[LAST]/l[LAST]));
			double A5 = -4*eps*5.67E-08*(thau/CpRho)*pow(Twr(), 3.)/(l[LAST-1]*(w[LAST-1]/l[LAST-1]+w[LAST]/l[LAST]));
			double A6 = -4*eps*5.67E-08*(thau/CpRho)*(pow(Twr(), 4.)/(4.*w[LAST-1])-pow(Twr(), 4.)/w[LAST-1]);
			A[i] += 0.;
			C[i] += A4;
			B[i] += A5;
			F[i] += A6;
		}
	}
	/* Итеративное выполнение прогонки до тех пор, пока не устоится температура внешней стенки. */
	TDMA();
}
