/**
 * @file math_addons.h
 * @brief Математические операции и константы.
 * @version 0.1
 * @copyright MIT License
 * @author A.V. Lepikhov
 * @date 2012
 */

#ifndef _MATH_ADDONS_H 
#define _MATH_ADDONS_H

#include <math.h>

/** Значение числа Пи, полученное из стандартного калькулятора MS Windows XP*/
#define PI	(3.1415926535897932384626433832795)

/** Постоянная Стефана-Больцмана */
#define SBconst 5.67032 // [Вт*м^-2*К^-4]

/**
 * @brief Тип данных общего назначения для хранения функций интерполяции.
 * @details За корректностью создания и заполнения структуры следит пользователь.
 */
typedef struct {
	/** Массив значений аргумента. */
	double *x;
	/** Массив значений функции. */
	double *y;
	/** количество точек интерполяции. */
	int count;
} func_points_t;

/* Для сокращения записи */
/** Тангенс угла в градусах */
#define TAN(theta)	tan(rad(theta))
/** Косинус угла в градусах */
#define COS(theta)	cos(rad(theta))
/** Синус угла в градусах */
#define SIN(theta)	sin(rad(theta))
/* --------------------- */

/** 
 * @brief Функция вычисления минимума
 * @param a - первый аргумент.
 * @param b - второй аргумент.
 */
extern double min(double a, double b);
/** 
 * @brief Функция вычисления максимума
 * @param a - первый аргумент.
 * @param b - второй аргумент.
 */
extern double max(double a, double b);

/**
 * @brief Угол angle в радианах.
 * @param angle - угол в градусах.
 */
extern double rad(double angle);

/**
 * @brief Угол angle в градусах.
 * @param angle - угол в радианах.
 */
extern double grad(double angle);

/** 
 * @brief Квадратичная интерполяция
 * @details x0[], y[0] - массив известных значений x и соответствующих им значений y
 * x - точка, для которой необходимо вычислить y путем квадратичной интерполяции
 * @return значение y
 */
extern double in_Quadratic(const double x0[3], const double y0[3], const double x);

/**
 * @brief Линейная интерполяция
 * @details x0[], y[0] - массив известных значений x и соответствующих им значений y
 * x - точка, для которой необходимо вычислить y путем линейной интерполяции
 * @return значение y
*/
extern double in_linear(const double x0[2], const double y0[2], const double x);

/**
 * @brief Линейная интерполяция по заданному набору точек.
 * @details Для корректной работы функции, точки предварительно должны быть отсортированы.
 * Экстраполяция влево и вправо от диапазона известных значений не допускается (возвращается крайнее левое и крайнее
 * правое значение соответственно).
 * @param points - указатель на таблицу интерполяции. Здесь x - значение аргумента; y - значение функции; count - количество элементов в таблице.
 * @param x - значение аргумента, для которого необходимо получить значение функции.
 * @param method - метод интерполяции:
 * 0 - линейная;
 * 1 - по квадрату числа маха (М^-2)
 */
extern double in_LinearFunc(const func_points_t* points, double x, int method = 0);

/**
 * @brief Вычислить значение определённого интеграла в интервале [0...x].
 * @param points - набор известных значений функции.
 * @param x - верхний предел интегрирования.
 * @details Для расчёта промежуточных значений используется линейная интерполяция.
 */
extern double ma_integral(const func_points_t* points, const double x);

#endif /* MATH_ADDONS_H */

