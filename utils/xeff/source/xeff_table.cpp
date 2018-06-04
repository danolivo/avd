#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include <xeff_table.h>
#include <cmath>
#include <assert.h>

using namespace std;

CXeff_table::CXeff_table(double M, vector<pcell_t*> pressure, int PRESSURE_POINTS_NUM, bool HasVelocityData) {
	const double hi = 1.4;
	
	this->HasVelocityData = HasVelocityData;
	// 1. Выделить память под массив точек интегрирования для расчёта xeff
	turb_points.x = new double[PRESSURE_POINTS_NUM];
	turb_points.y = new double[PRESSURE_POINTS_NUM];
	turb_points.count = PRESSURE_POINTS_NUM;
	lam_points.x = new double[PRESSURE_POINTS_NUM];
	lam_points.y = new double[PRESSURE_POINTS_NUM];
	lam_points.count = PRESSURE_POINTS_NUM;
	assert(PRESSURE_POINTS_NUM > 0);
	assert((turb_points.x != 0) && (turb_points.y != 0));
	assert((lam_points.x != 0) && (lam_points.y != 0));
	
	// 2. Подготавливаем подинтегральную функцию
	double PdivP0_touch;
	double rm;
	
	// Заголовок печатаемой таблицы
	printf("\npoint\tVdivVh\tPdivP0\tft\tfl\n");
	for (int i=0; i<PRESSURE_POINTS_NUM; i++)
	{
		double VdivVh;
		
		turb_points.x[i] = pressure[i]->length;
		lam_points.x[i] = pressure[i]->length;
		assert(pressure[i]->length >= 0.);
		assert((pressure[i]->PdivP0 > 0.) && (pressure[i]->PdivP0 <= 1.2));
		/*
		 * 3. Расчёт отношения скорости набегающего на тело потока к скорости на
		 * бесконечности
		*/	
		PdivP0_touch = pressure[i]->PdivP0;
		rm = pressure[i]->rm;
		if (HasVelocityData)
			VdivVh = pressure[i]->VdivVh;
		else
			VdivVh = sqrt(((2./(hi-1))*(1/(M*M))+1)*(1-pow(PdivP0_touch, (hi-1)/hi)));

		if (VdivVh < 0.)
		{
			printf("[EE]: Velocity is negative! Check input pressure table!\n");
			exit(-1);
		}
		assert(VdivVh >=0);
		/*
		 * 4. Расчёт подинтегральной функции.
		 */
		turb_points.y[i] = PdivP0_touch * VdivVh * pow(rm, 5./4.);
		lam_points.y[i] = PdivP0_touch * VdivVh * pow(rm, 2.);
		printf("%d\t%lf\t%lf\t%lf\t%lf\n", i, VdivVh, pressure[i]->PdivP0, turb_points.y[i], lam_points.y[i]);
	}
	printf("\nXeff function loaded\n");

	return;
};

double CXeff_table::calculate(length_t length, bool IsTurbulent) {

	assert((length >=0) );
//	printf("length=%E turb_points.count=%d x=%E\n", length, turb_points.count, turb_points.x[turb_points.count-1]);
	if (turb_points.x[turb_points.count-1] < length)
		return 0;
	if (lam_points.x[lam_points.count-1] < length)
		return 0;
		
	if (length == 0)
		return 0;
//	printf("in_LinearFunc(&lam_points, length)=%E\n", in_LinearFunc(&lam_points, length));
	if (IsTurbulent)	
		return ma_integral(&turb_points, length)/in_LinearFunc(&turb_points, length);
	else
		return ma_integral(&lam_points, length)/in_LinearFunc(&lam_points, length);
}
