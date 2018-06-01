/**
 * @brief Решатель тепловой задачи с уносом, аналогичный по логике работы программе Авдуевский.
 * @author А.В. Лепихов
 * @date 2013 г.
 * @copyright MIT License
 */
#ifndef _AVDSOLVER_H_
#define _AVDSOLVER_H_

#include <common.h>
#include <avdtparser.h> // Парсер файла ИД
#include <thsolver.h> // тепловой решатель
#include <model.h> // тепловая модель
#include <bluntedcone.h>


typedef struct {
	double I0, IW, IE, ISTAR, ALC, ALC1, QCONV, XAP1, P0, P1, V0, F1, F2, F3, KDIS, KENTH;
} avd_t;

/** Моментальный снимок состояния расчёта. */
typedef struct {
	/** Момент времени, для которого сделан снимок. */
	double time;
	/** Высота, м */
	double H;
	/** Скорость, м/с*/
	double V;
	/** Угол атаки, град. */
	double al;
	/** Угол проворота, град. */
	double phi;
	/** Число Маха */
	double mach;
	/** Эффективная длина ^0.2 */
	double XEF;
	/** Коэффициент давления P/P0. */
	double PP0;
	avd_t avd;
	/** Коэффициент массообмена. */
//	double acp;
	/** Определяющая энтальпия, Дж/кг. */
//	double Ie;
	/** Энтальпия на стенке, Дж/кг. */
//	double Iw;
	/** Полная энтальпия, Дж/кг. */
//	double I0;
	/** Давление на стенке, Па. */
//	double Ps;
	/** Безразмерный параметр уноса. */
	double G;
	/** Результат расчёта теплопередачи. */
	solve_result_t srt;
} avdsolver_t;
/** Решатель одномерной тепловой задачи с уносом, аналогичный по логике работы программе Авдуевский. */
class AVDSolver {
	/** Указатель на газодинамические параметры в окрестности точки. */
	gasdynamics_t* gd;
	/** Указатель на структуру с траекторией. */
	trm_t* trm;
	/** Указатель на тепловую модель. */
	thm_t* thm;
	/** Указатель на класс теплового решателя. */
	CTHSolver* thsolver;
	/** Шаг счёта. */
	double TIMESTEP;

	CBluntedCone* BCone;
public:
	/** Конструктор класса. */
	AVDSolver(thm_t* thm, trm_t* trm, gasdynamics_t* gd, CBluntedCone* BCone);
	/**
	 * @brief Выполнить расчет.
	 * @param time - момент времени, до которого необходимо выполнить расчет.
	 */
	avdsolver_t Solve(double time);
	/** Деструктор класса. */
	~AVDSolver();
	/** Вывод состояния класса. */
	void print();
};

#endif /* _AVDSOLVER_H_ */
