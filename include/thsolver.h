#ifndef _THSOLVER_H_
#define _THSOLVER_H_

#include <common.h>
#include <model.h>

/** Структура с индикаторами результата расчета. */
typedef struct {
	/** Интегральное количество теплоты (конвективная составляющая), которое пришло на левую стенку за посчитанный отрезок расчетного времени, [Дж*м^-2]. */
	double QLconv;
	/** Интегральное количество теплоты (радиационная составляющая), которое пришло на левую стенку за посчитанный отрезок расчетного времени, [Дж*м^-2]. */
	double QLrad;
	/** Интегральное количество теплоты (конвективная составляющая), которое пришло на правую стенку за посчитанный отрезок расчетного времени, [Дж*м^-2]. */
	double QRconv;
	/** Интегральное количество теплоты (радиационная составляющая), которое пришло на правую стенку за посчитанный отрезок расчетного времени, [Дж*м^-2]. */
	double QRrad;
	/** Изменение количества теплоты внутри расчетной области за посчитанный отрезок расчетного времени, [Дж*м^-2]. */
	double dHeatQty;
	/** Последний использованный шаг счёта. */
	double LAST_TIMESTEP;
	/** Maximum temperature change at the timestep. */
	double CURRENT_DT_MAX;
} solve_result_t;

/** Реализует расчет теплопередачи в твердом теле в одномерной постановке. */
class CTHSolver {
	/** Указатель на тепловую модель. */
	thm_t *thm;
	/** Массив температур, К. */
	double T[CELLS_MAX_NUM+2];
	/** Массив длин ячеек, м.*/
	double w[CELLS_MAX_NUM+2];
	/** Массив теплопроводностей в ячейках. */
	double l[CELLS_MAX_NUM+2];
	/** массив теплоёмкостей. */
	double c[CELLS_MAX_NUM+2];
	/** Массив плотностей. */
	double r[CELLS_MAX_NUM+2];
	/** Массив источников тепла по ячейкам. */
	double qv[CELLS_MAX_NUM+2];
	/** Количество ячеек в расчётной модели. */
	int SIZE;
	/** Степень черноты на левой и правой границе модели. */
	double eps[2];
	/** Минимальный шаг счёта. */
	double TIMESTEP_MIN;
	/** Максимальный шаг счёта. */
	double TIMESTEP_MAX;
	/** Maximum temperature change per timestep or 0 */
	double DT_MAX;
	/** Result of solving process. */
	solve_result_t sres;
public:
	/** 
	 * @brief Class constructor. 
	 * @param thm - pointer to the thermal model.
	 */
	CTHSolver(thm_t* thm);
	/** Деструктор класса. */
	~CTHSolver();
	/**
	 * @brief Solve thermal task from CURRENT_TIME to time
	 * @param time - next stop point of timeline
	 * @return Heat quantity growth at the model during calculation time, [J]
	 */
	solve_result_t Solve(double time);
	/** Подготавливает тепловую модель к расчету. Подразумевается, что здесь считается только прогрев. */
	void Prepare();
	/**
	 * @brief Operations before solving iteration
	 * @return Current value of heat quantity [J]
	 */
	double Pre(double time);
	/**
	 * @brief Operations after solving iteration
	 * @return Current value of heat quantity [J]
	 */
	double Post(double time);
	/** Установка настроек солвера.
	 * @param TIMESTEP_MIN - minimum timestep
	 * @param TIMESTEP_MAX - maximum timestep
	 * @param DT_MAX - maximum temperature change per timestep or 0 another
	 */
	void setPrefs(double TIMESTEP_MIN, double TIMESTEP_MAX, double DT_MAX = 1.);
	/**
	 * @brief Do calculation of time step.
	 * @param TIMESTEP - desired time step.
	 * @return Result of iteration - SUCCESS or NO.
	 */
	int DoIteration(double TIMESTEP);
	/**
	 * Init left and right boundary values.
	 * @param lbc - pointer to the left boundary.
	 * @param lbc - pointer to the right boundary.
	 * @param time - current time for the boundary values.
	 */
	void setBoundaries(CBoundary *lbc, CBoundary *rbc, double time);
	/**
	 * Init left boundary.
	 * @param cbc - pointer to the boundary.
	 * @param time - current time for the boundary value.
	 */
	void setLeftBoundary(CBoundary *cbc, double time);
	/**
	 * Init right boundary.
	 * @param cbc - pointer to the boundary.
	 * @param time - current time for the boundary value.
	 */
	void setRightBoundary(CBoundary *cbc, double time);
	/**
	 * return pointer to the thermal model.
	 */
	thm_t* getTHM();
	/** Print preferences */
	void printPreferences();
	/** Return heat quantity. */
	double currentHeatQty();
	/** Print thermal solver state. */
	virtual void print(FILE* out = stdout);
protected:
	/** Return left boundary temperature, K. */
	double Twl();
	/** Return right boundary temperature, K. */
	double Twr();
};

#endif /* _THSOLVER_H_ */
