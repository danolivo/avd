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

#include <cmath>

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
extern double MIN(double a, double b);
/** 
 * @brief Функция вычисления максимума
 * @param a - первый аргумент.
 * @param b - второй аргумент.
 */
extern double MAX(double a, double b);

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

extern double in_Quadratic2(double x0, double y0, double x1, double y1, double x2, double y2, double x);
/**
 * @brief Линейная интерполяция
 * @details x0[], y[0] - массив известных значений x и соответствующих им значений y
 * x - точка, для которой необходимо вычислить y путем линейной интерполяции
 * @return значение y
*/
extern double in_linear(const double x0[2], const double y0[2], const double x);
extern double in_linear(double x0, double y0, double x1, double y1, double x);
/**
 * @brief Copy interpolation table from src to dst
 * @param src - source interpolation table
 * @param dst - destination interpolation table
 * @return -1 - if dst is NULL; 1 - if src is NULL; 0 - if copying was successed.
 */
extern int func_points_cpy(const func_points_t* src, func_points_t* dst);
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

#include <vector>
using namespace std;

typedef struct {
	double valF;
	double valX;
} func_t;
/** Класс предоставляет инструменты управления функцией интерполяции */
class IFunc {
	/** Cтруктура для хранения значений функции */
	vector<func_t *> func1D;
public:
	/** Конструктор класса. */
	IFunc();
	/**
	 * @brief Добавить точку интерполяции для функции
	 * @param valX - значение аргумента.
	 * @param valF - значение функции
	 */
	void add(double valX, double valF);
	/** Выдать значение функции интерполяции
	 * @param x - значение аргумента функции интерполяции.
	 */
	double val(double x);
	/** Деструктор класса. */
	void print();
	~IFunc();
};

typedef struct {
	/** Указатель на функцию интерполяции одного переменного. */
	IFunc* func;
	/** Значение аргумента. */
	double valY;
} table_t;
/** Класс предоставляет инструменты интерполяции физической величины по двум параметрам. */
class ITable {
	/** Хранилище аргументов и значений функции. */
	vector<table_t*> func2D;
public:
	/** Конструктор класса. */
	ITable();
	/**
	 * @brief Добавить набор значений функции для второго аргумента
	 * @param func - указатель на функцию одного аргумента.
	 * @param valY - значение второго аргумента
	 */
	void add(IFunc* func, double valY);
	/** Выдать значение функции интерполяции
	 * @param x - значение первого аргумента функции интерполяции.
	 * @param y - значение второго аргумента функции интерполяции.
	 */
	double val(double x, double y);

	/** Деструктор класса. */
	~ITable();
};
typedef struct {
	/** Указатель на функцию интерполяции двух переменных. */
	ITable* func;
	/** Значение аргумента. */
	double valZ;
} cube_t;
class ICube {
	/** Хранилище аргументов и значений функции. */
	vector<cube_t*> func3D;
public:
	/** Конструктор класса. */
	ICube();
	/**
	 * @brief Добавить набор значений функции для второго аргумента
	 * @param valZ - значение третьего аргумента.
	 * @param tbl - набор значений двумерной табличной функции
	 */
	void add(ITable* tbl, double valZ);
	/** Выдать значение функции интерполяции
	 * @param x - значение первого аргумента функции интерполяции.
	 * @param y - значение второго аргумента функции интерполяции.
	 * @param y - значение третьего аргумента функции интерполяции.
	 */
	double val(double x, double y, double z);

	/** Деструктор класса. */
	~ICube();
};
#endif /* MATH_ADDONS_H */

