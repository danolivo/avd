#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <math_addons.h>

double min(double a, double b) {
	if (a<=b)
		return a;
	else
		return b;
}

double max(double a, double b) {
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
	
	assert(x != NULL);
	assert(y != NULL);
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
	}
	return 0.;
}
/*
int main() {
	func_points_t points;
	
	points.x = (double*)calloc(3, sizeof(double));
	points.y = (double*)calloc(3, sizeof(double));
	points.count = 3;
	points.x[0] = -5;
	points.x[1] = 1;
	points.y[0] = -5;
	points.y[1] = 1;

	points.x[2] = 2;
	points.y[2] = 3;
//	printf("F(x) = %f\n", in_LinearFunc(&points, 1.5));
	printf("integral = %f\n", ma_integral(&points, 2));
	
	return 0;
}
*/
