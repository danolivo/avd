/**
 * @file xeff_table.h
 * @brief Расчёт х эффективного для метода эффективной длины В.С. Авдуевского
 * @details Данный класс предназначен для расчёта x эффективного для случая, когда
 * имеется массив давлений в точках боковой поверхности ЛА.
 * Применение: для ЛА сложной формы, для которых нельзя использовать эмпирические
 * методики расчёта давления и трудно программировать класс геометрии.
 * @version 0.1
 * @copyright MIT License
 * @author А.В. Лепихов
 * @date 07.12.2011
 */

#ifndef _XEFF_TABLE_H
#define _XEFF_TABLE_H

#include <vector>
#include <math_addons.h>
#include <magnitude.h>

using namespace std;

/**
 * @brief Структура данных для параметров интегрирования X эффективного.
 */
typedef struct {
	/** Расстояние вдоль образующей от критической точки */
	length_t length;
	/** Относительное давление в точке */
	double PdivP0;
	/** Радиус миделя в точке */
	length_t rm;
	
	double VdivVh;
} pcell_t;

/**
 * @brief Предоставляет инструмент расчёта X эффективного по параметрам ЛА, записанным в таблице.
 */
class CXeff_table {
	/** Точки интегрирования */
	func_points_t turb_points;
	func_points_t lam_points;
	bool HasVelocityData;
public:
	/**
	 * @brief Конструктор класса.
	 * @param M - число Маха.
	 * @param pressure - таблица параметров расчёта X эффективного в точке ЛА.
	 * @param PRESSURE_POINTS_NUM - количество точек интерполяции.
	 */
	CXeff_table(double M, vector<pcell_t*> pressure, int PRESSURE_POINTS_NUM, bool HasVelocityData);
	
	/**
	 * @brief Расчёт Х эффективного для заданной точки боковой поверхности.
	 * @details Экстраполяция Х эффективного не допускается. В случае попытки экстраполяции возвращается 0.
	 * @param length - расстояние от критической точки вдоль боковой поверхности.
	 * @param IsTurbulent - Режим течения. true - турбулентный режим, иначе ламинарный.
	 * @return Х эффективное.
	 */
	double calculate(length_t length, bool IsTurbulent = true);
};

#endif /* _XEFF_TABLE_H */
