#include <assert.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <math_addons.h>

double MIN(double a, double b) {
	if (a<=b)
		return a;
	else
		return b;
}

double MAX(double a, double b) {
	if (a<=b)
		return b;
	else
		return a;
}

// Перевод угла из градусов в радианы
// angle - угол в градусах
// return: угол радианах
double rad(double angle) {
	return PI*angle/180;
}

// Перевод угла из радиан в градусы
// angle - угол в радианах
// return: угол в градусах
double grad(double angle) {
	return 180*angle/PI;
}

/*
 * Расчёт определенного интеграла по линейной функции
 * Значения функции задаются в пределах интегрирования
 */
double ma_2pIntegral(const double* x, const double* y) {
	double N = 3;
	double dx, _integral;
	double sum=0;
	int i;
	
	assert(x != 0);
	assert(y != 0);
	if (x[0] == x[1]) // Может возникнуть при автоматизированной обработке больших массивов данных из-за ошибок округления
		return 0;
//	assert((x[0] != x[1]) || !printf("x[0]: %f\t x[1]: %f\n", x[0], x[1]));

	sum = -1; // фиктивное значение для первого шага	
	do {
		_integral = sum;
		N *=2;
		dx = fabs((x[1] - x[0])/N);
		sum = 0;
		for (i = 0; i < N; i++)
			sum += fabs(in_linear(x, y, (x[0]+dx*i)+(dx/2))*dx);

		assert(N < 10000);
	} while (fabs(_integral - sum) > 0.001);

	return sum;
}

/*
 * Вычислить значение определённого интеграла в интервале [0...x]
 * func - набор известных значений функции
 * Для расчёта промежуточных значений используем линейную интерполяцию
 */
double ma_integral(const func_points_t* points, const double x) {
	int i=0, j;
	double sum=0;
	
	/* Защитные инструкции */
	assert(points != NULL);
	assert(points->x != NULL);
	assert(points->y != NULL);
	assert(points->count > 0);
	/* ------------------- */

	/* Посчитать первый отрезок интеграла (от нуля до первой точки интерполяции)*/
	sum = points->x[0]*points->y[0];
	/*
	 * Поиск позиции точки x в наборе известных точек (слева, справа, между)
	 */
	for (i=0; (i < points->count) && (x > points->x[i]); i++);
	
	if (i == 0) { // Точка слева от известного набора точек.
		sum -= (points->x[0]-x)*points->y[0];
		return sum;
	}

	/* Посчитать внутреннюю область интеграла. */
	for (j = 1; j < i; j++)
		sum += (points->x[j]-points->x[j-1])*in_LinearFunc(points, 0.5*(points->x[j]+points->x[j-1]));
	/* Посчитать последний отрезок. */
	sum += (x - points->x[j-1])*in_LinearFunc(points, 0.5*(x+points->x[j-1]));
	
	return sum;
}

/* 
 * Квадратичная интерполяция
 * x0[], y[0] - массив известных значений x и соответствующих им значений y
 * x - точка, для которой необходимо вычислить y путем квадратичной интерполяции
 * return значение y
 */
double in_Quadratic(const double x0[3], const double y0[3], const double x) {

	return (x-x0[1])*(x-x0[2])/(x0[0]-x0[1])*(x0[0]-x0[2])*y0[0]+(x-x0[0])*(x-x0[2])/(x0[1]-x0[0])*(x0[0]-x0[2])*y0[1]+(x-x0[0])*(x-x0[1])/(x0[1]-x0[0])*(x0[2]-x0[1])*y0[2];
/**	
	a = ((y0[2]-y0[0])*(x0[1]-x0[0])-(y0[1]-y0[0])*(x0[2]-x0[0]))/((x0[2]*x0[2]-x0[0]*x0[0])*(x0[1]-x0[0])-(x0[1]*x0[1]-x0[0]*x0[0])*(x0[2]-x0[0])); //d*x0[2]-y0[0]+d*x0[0])/(x0[2]*x0[2]-x0[0]*x0[0]+(x0[1]-x0[0])*(x0[0]-x0[2]));
	b = (y0[1] - y0[0] - a*(x0[1]*x0[1]-x0[0]*x0[0]))/(x0[1]-x0[0]);
	c = y0[0] - a*x0[0]*x0[0] - b*x0[0];
	printf("a=%f, b=%f, c=%f\n",a, b, c);
	return a*x*x + b*x + c; */
}

double in_Quadratic2(double x0, double y0, double x1, double y1, double x2, double y2, double x)
{
	assert(x1 > x0); assert(x2 > x1);
	double tmp = (y1-y0)/(x1-x0);
	double a = (y2-tmp*x2-y0+tmp*x0) / (x2*x2-x2*(x1+x0)-x0*x0+(x1+x0)*x0);
	double b = tmp - a*(x1+x0);
	double c = y0 - a*x0*x0 - b*x0;
/*	if ((y0 < y1) || (y1 < y2)) {
		printf("x0=%E y0=%lf x1=%lf y1=%lf x2=%lf y2=%lf x=%E f=%E\n", x0, y0, x1, y1, x2, y2, x, a*x*x+b*x+c);
		exit(-1);
	} else */
		return a*x*x+b*x+c;
}
/*
 * Линейная интерполяция
 * x0[], y[0] - массив известных значений x и соответствующих им значений y
 * x - точка, для которой необходимо вычислить y путем линейной интерполяции
 * return значение y
*/
double in_linear(const double x0[2], const double y0[2], const double x) {

	if ((x0[1]-x0[0]) == 0) {
	   assert((y0[1]-y0[0]) == 0);	
	   return y0[0];
 	   }
//	printf("x=%f, y0=%f, y1=%f x0=%f, x1=%f\n", x, y0[0],y0[1], x0[0], x0[1]);
	return y0[0]+((y0[1]-y0[0])*(x-x0[0]))/(x0[1]-x0[0]);
}
double in_linear(double x0, double y0, double x1, double y1, double x) {

	if ((x1-x0) == 0) {
	   assert((y1-y0) == 0);	
	   return y0;
 	   }
	return y0+((y1-y0)*(x-x0))/(x1-x0);
}

/**
 * @brief Линейная интерполяция по заданному набору точек
 * Для корректной работы функции, точки предварительно должны быть отсортированы
 * по x
 */
double in_LinearFunc(const func_points_t* points, const double x, int method) {
	int i=0;
	int i0;

	/* Защитные инструкции */
	assert(points != NULL);
	assert(points->x != NULL);
	assert(points->y != NULL);
	assert(points->count > 0);
	/* ------------------- */
	
	if (points->count == 1)
		return points->y[0];
	/*
	 * Найти интервал, в который попадает точка x
	 */
	while ((i<points->count) && (x > points->x[i])) i++;
	
	/* Если точка левее последней известной точки - используем первую точку */
	if (i == 0)
		return points->y[0];
		
	/* Если точка правее последней известной точки - используем последнюю точку. */
	else if (i == points->count)
		return points->y[points->count-1];
		
	/* Если точка находится в интервале между двух известных точек */
	else
		i0 = i-1;
		
	double i_M1, i_M2, i_M;
	switch (method) {
	case 0:
		return in_linear(&(points->x[i0]), &(points->y[i0]), x);
	case 1:
		i_M1 = 1./pow(points->x[i0], 2);
		i_M2 = 1./pow(points->x[i0+1], 2);
		i_M = 1./pow(x, 2);
		return points->y[i0] + (points->y[i0+1]-points->y[i0])*(i_M-i_M1)/(i_M2-i_M1);
	default:
		assert(0);
		break;
	}
	return 0.;
}

int func_points_cpy(const func_points_t* src, func_points_t* dst) {
	if (dst == 0)
		return -1;
	if (src == 0)
		return 1;
	dst->count = src->count;
	dst->x = (double *)calloc(dst->count, sizeof(double));
	dst->y = (double *)calloc(dst->count, sizeof(double));
	for (int i=0; i<dst->count; i++) {
		dst->x[i] = src->x[i];
		dst->y[i] = src->y[i];
	}
	return 0;
}
IFunc::IFunc()
{
}
void IFunc::add(double valX, double valF)
{
	func_t* t = new func_t;
	t->valX = valX;
	t->valF = valF;
	func1D.push_back(t);
}
double IFunc::val(double x)
{
	assert(func1D.size() > 0);
	if ((func1D[0]->valX-x)*(func1D[func1D.size()-1]->valX-x) >= 0.) // Значение выходит за границу диапазона интерполяции
		if (fabs(func1D[0]->valX-x) < fabs(func1D[func1D.size()-1]->valX-x))
			return func1D[0]->valF;
		else
			return func1D[func1D.size()-1]->valF;
	if (func1D.size() == 1)
		return func1D[0]->valX;
	for (int i=1; i<(int)func1D.size(); i++) {
//		printf("\nVX=%lf\t", func1D[i]->valX);
		if ((func1D[i]->valX-x)*(func1D[i-1]->valX-x) <= 0.) {
//			printf("x=%lf valX=%lf valX2=%lf i=%d\n", x, func1D[i-1]->valX, func1D[i]->valX, i);
			double P1, P2, valX1, valX2;
			P1 = func1D[i-1]->valF;
			valX1 = func1D[i-1]->valX;
			P2 = func1D[i]->valF;
			valX2 = func1D[i]->valX;
			return P1 + (P2-P1)*(x-valX1)/(valX2-valX1);
		}
	}
	assert(0);
	return 0;
}
void IFunc::print()
{
	printf("X=\t");
	for (int i=0; i<(int)func1D.size(); i++)
		printf("%11.5lf\t", func1D[i]->valX);
	printf("\nF=\t");
	for (int i=0; i<(int)func1D.size(); i++)
		printf("%11.5lf\t", func1D[i]->valF);
	printf("\n");
}
IFunc::~IFunc()
{
	func1D.clear();
}
ITable::ITable()
{
}
void ITable::add(IFunc* func, double valY)
{
	assert(func != 0);
	table_t *t = new table_t;
	t->func = func;
	t->valY = valY;
	func2D.push_back(t);
}
double ITable::val(double x, double y)
{
	assert(func2D.size() > 0);
//	printf("x=%lf y=%lf %lf %lf\n", x,y, func2D[0]->valY, func2D[func2D.size()-1]->valY);
	if ((func2D[0]->valY-y)*(func2D[func2D.size()-1]->valY-y) >= 0.) {
		if (fabs(func2D[0]->valY-y) < fabs(func2D[func2D.size()-1]->valY-y))
			return func2D[0]->func->val(x);
		else
			return func2D[func2D.size()-1]->func->val(x);
	}
//	assert(func2D[func2D.size()-1]->valY >= y);
	
	if (func2D.size() == 1)
		return func2D[0]->func->val(x);
//	printf("START\n");
	for (int i=1; i<(int)func2D.size(); i++)
		if ((func2D[i]->valY - y)*(func2D[i-1]->valY - y) <= 0.) {
			double P1, P2, valY1, valY2;
			P1 = func2D[i-1]->func->val(x);
			valY1 = func2D[i-1]->valY;
			P2 = func2D[i]->func->val(x);
			valY2 = func2D[i]->valY;
//			printf("valY1=%lf valY2=%lf P1=%lf P2=%lf i=%d\n", valY1, valY2, P1, P2, i);
			return P1 + (P2-P1)*(y-valY1)/(valY2-valY1);
		}
	assert(0);
	return 0;
}
ITable::~ITable()
{
	func2D.clear();
}
ICube::ICube()
{
}
void ICube::add(ITable* tbl, double valZ)
{
	assert(tbl != 0);
	cube_t *t = new cube_t;
	t->func = tbl;
	t->valZ = valZ;
	func3D.push_back(t);
}
double ICube::val(double x, double y, double z)
{
	assert(func3D.size() > 0);
	if (func3D[0]->valZ > z) {
		for (int i=0; i<(int)func3D.size(); i++)
			printf("valZ %d=%lf", i, func3D[i]->valZ);
		printf("\n");
	}
	assert((func3D[0]->valZ <= z) || (!printf("z=%lf (func3D[0]->valZ=%lf func3D.size() = %ld\n", z, func3D[0]->valZ, func3D.size())));
	assert(func3D[func3D.size()-1]->valZ >= z);
	if (func3D.size() == 1)
		return func3D[0]->func->val(x, y);
	for (int i=1; i<(int)func3D.size(); i++)
		if (func3D[i]->valZ >= z) {
			double P1, P2, valZ1, valZ2;
			P1 = func3D[i-1]->func->val(x, y);
			valZ1 = func3D[i-1]->valZ;
			P2 = func3D[i]->func->val(x, y);
			valZ2 = func3D[i]->valZ;
			return P1 + (P2-P1)*(z-valZ1)/(valZ2-valZ1);
		}
	assert(0);
	return 0;
}
ICube::~ICube()
{
	func3D.clear();
}
