/**
 * @file boundary.h
 * @brief Boundary interface
 * @version 0.1
 * @copyright MIT License
 * @author A.V. lepikhov
 * @date 2012
 */

#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include <common.h>
#include <math_addons.h>
#include <model.h>
/**
 * @brief Структура данных, возвращаемых классом ГУ.
 */
typedef struct {
	/** Heat transfer coefficient */
	double acp;
	/** Характерная энтальпия. */
	double Ie;
	/** Wall enthalpy. */
	double Iw;
	/** Total enthalpy, [J/kg]. */
	double I0;
	/** Boundary emissivity */
	double eps;
	/** Wall temperature */
	double Tw;
	/** Pressure at the wall, [Pa]. */
	double P;
} BC_t;
/**
 * @brief Abstract boundary
 * @details Is an ancestor for any boundaries
 */
class CBoundary {
	int Type;
protected:
	/** Base parameters of second-order boundary condition */
	BC_t BC;
public:
	/** Output device for logging */
	FILE* flog;
	/** Class Constructor.
	 * @param Type -type of boundary.
	 * @param flog - pointer to the output device
	 */
	CBoundary(int Type, FILE* flog = stdout);
	/** Class constructor
	 * @param boundary - pointer to the class, which will be assigned to the new class.
	 */
	CBoundary(CBoundary *boundary);
	/** Returns class type. */
	int type();
	/** Calculate boundary parameters.
	 * @param Time - current time, s.
	 */
	virtual BC_t GetBC(double Time) = 0;
	/**
	 * @brief Set time-dependent boundary parameter changes by multipliers.
	 * @param pname - boundary parameter name
	 * @param table - pointer to interpolation table for multiplier
	 * @return Opcode. NULL_VALUE - interpolation table can't be used. -2 - pname is not valid
	 */
	virtual int SetTimeChangeCoefs(const char* pname, const func_points_t* table) = 0;
	/**
	 * @brief Print class preferences.
	 */
	virtual void print() = 0;
};
/** Second-order boundary. */
class CSOBoundary : public CBoundary {
public:
	/** 
	 * @brief Class constructor
	 * @param eps - base value of emissivity
	 * @param Ie - base value of local boundary layer edge enthalpy
	 * @param acp - base value of GHTC
	 * @param flog - output device for logging
	 */
	CSOBoundary(double acp, double I0, double Ie, double Iw, double Ps, double eps = 0., FILE* flog = stdout);
	/**
	 * @brief Выдать параметры граничного условия.
	 * param Time - текущий момент расчётного времени.
	 */
	virtual BC_t GetBC(double Time);
	/**
	 * @brief Set time-dependent boundary parameter changes by multipliers.
	 * @param pname - boundary parameter name
	 * @param table - pointer to interpolation table for multiplier
	 * @return Opcode. NULL_VALUE - interpolation table can't be used. -2 - pname is not valid
	*/
	virtual int SetTimeChangeCoefs(const char* pname, const func_points_t* table);
	/** Распечатать внутренние данные класса. */
	virtual void print();
};
#endif /* _BOUNDARY_H_ */
