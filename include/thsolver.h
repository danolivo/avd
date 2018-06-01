#ifndef _THSOLVER_H_
#define _THSOLVER_H_

#include <common.h>
#include <model.h>

/** ��������� � ������������ ���������� �������. */
typedef struct {
	/** ������������ ���������� ������� (������������ ������������), ������� ������ �� ����� ������ �� ����������� ������� ���������� �������, [��*�^-2]. */
	double QLconv;
	/** ������������ ���������� ������� (������������ ������������), ������� ������ �� ����� ������ �� ����������� ������� ���������� �������, [��*�^-2]. */
	double QLrad;
	/** ������������ ���������� ������� (������������ ������������), ������� ������ �� ������ ������ �� ����������� ������� ���������� �������, [��*�^-2]. */
	double QRconv;
	/** ������������ ���������� ������� (������������ ������������), ������� ������ �� ������ ������ �� ����������� ������� ���������� �������, [��*�^-2]. */
	double QRrad;
	/** ��������� ���������� ������� ������ ��������� ������� �� ����������� ������� ���������� �������, [��*�^-2]. */
	double dHeatQty;
	/** ��������� �������������� ��� �����. */
	double LAST_TIMESTEP;
	/** Maximum temperature change at the timestep. */
	double CURRENT_DT_MAX;
} solve_result_t;

/** ��������� ������ ������������� � ������� ���� � ���������� ����������. */
class CTHSolver {
	/** ��������� �� �������� ������. */
	thm_t *thm;
	/** ������ ����������, �. */
	double T[CELLS_MAX_NUM+2];
	/** ������ ���� �����, �.*/
	double w[CELLS_MAX_NUM+2];
	/** ������ ����������������� � �������. */
	double l[CELLS_MAX_NUM+2];
	/** ������ ������������. */
	double c[CELLS_MAX_NUM+2];
	/** ������ ����������. */
	double r[CELLS_MAX_NUM+2];
	/** ������ ���������� ����� �� �������. */
	double qv[CELLS_MAX_NUM+2];
	/** ���������� ����� � ��������� ������. */
	int SIZE;
	/** ������� ������� �� ����� � ������ ������� ������. */
	double eps[2];
	/** ����������� ��� �����. */
	double TIMESTEP_MIN;
	/** ������������ ��� �����. */
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
	/** ���������� ������. */
	~CTHSolver();
	/**
	 * @brief Solve thermal task from CURRENT_TIME to time
	 * @param time - next stop point of timeline
	 * @return Heat quantity growth at the model during calculation time, [J]
	 */
	solve_result_t Solve(double time);
	/** �������������� �������� ������ � �������. ���������������, ��� ����� ��������� ������ �������. */
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
	/** ��������� �������� �������.
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
