/**
 * @file material.h
 *
 * @date 27.08.2013
 * @author: A.V. Lepikhov
 */

#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include <common.h>
/** Size of material short name. */
#define MATERIAL_NAME_MAX_LEN		(20)
/** Interface to the material properties. */
class CMaterial {
	/** Material short name. */
	char Name[MATERIAL_NAME_MAX_LEN];
protected:
	/** Log device pointer. */
	FILE* flog;
public:
	/** Class constructor.
	 * @param title - material short name.
	 * @param flog - pointer to the output device.
	 */
	CMaterial(const char* title, FILE* flog);
	/** Return short material name. */
	char* name();
	/** Heat Conductivity.
	 * @param T - current temperature.
	 */
	virtual double l(double T) = 0;
	/** Heat Specific heat.
	 * @param T - current temperature.
	 */
	virtual double c(double T) = 0;
	/** density.
	 * @param T - current temperature.
	 */
	virtual double r(double T) = 0;
	/** Emissivity.
	 * @param T - current temperature.
	 */
	virtual double eps(double T) = 0;
	/** Ablation A-coefficient.
	 * @param T - current temperature.
	 */
	virtual double a(double T) = 0;
	/** Ablation B-coefficient.
	 * @param T - current temperature.
	 */
	virtual double b(double T) = 0;
	virtual double Td(double T) = 0;
	virtual int at(void) = 0;
	/** Ablation Heat Integral from T=0 temperature.
	 * @param T - current temperature.
	 */
	virtual double heatQuantity(double T) = 0;
	/** Print material info.
	 * @param pointer to the output device.
	 */
	virtual void print(FILE *out) = 0;
};
/** Material with constant properties. */
class ConstMaterial: public CMaterial {
	double L;
	double RHO;
	double CP;
	double EPS;
	double Tdestr;
	double A;
	double B;
	double type;
public:
	/** Class constructor
	 * @param title - material short name
	 * @param flog - pointer to the output device
	 * @param L - heat conductivity
	 * @param RHO - density
	 * @param CP - specific heat
	 * @param EPS - emissivity
	 * @param Tdestr - destruction temperature
	 * @param A - ablation coefficient
	 * @param B - ablation coefficient
	 * @param Type - Ablation type.
	 */
	ConstMaterial(const char* title, FILE* flog, double L, double RHO, double CP, double EPS,
			double Tdestr, double A, double B, int Type);
	/** Heat Conductivity.
	 * @param T - current temperature.
	 */
	virtual double l(double T);
	/** Heat Specific heat.
	 * @param T - current temperature.
	 */
	virtual double c(double T);
	/** density.
	 * @param T - current temperature.
	 */
	virtual double r(double T);
	/** Emissivity.
	 * @param T - current temperature.
	 */
	virtual double eps(double T);
	/** Ablation A-coefficient.
	 * @param T - current temperature.
	 */
	virtual double a(double T);
	/** Ablation B-coefficient.
	 * @param T - current temperature.
	 */
	virtual double b(double T);
	virtual double Td(double T);
	virtual int at(void);
	/**
	 * @brief Heat Quantity at the reference temperature
	 * @param T - reference temperature
	 * @return Heat Quantity, [J*m^-3]
	 */
	virtual double heatQuantity(double T);
	/**
	 * @brief Check correctess of material thermophysics
	 * SUCCESS - if all defined correctly; NULL_VALUE, overwise
	 */
	virtual int check();
	/** Print material info.
	 * @param pointer to the output device.
	 */
	virtual void print(FILE *out = stdout);
};
class CUserMaterial : public ConstMaterial {
	/** Heat conductivity interpolation function. */
	func_points_t fL;
	/** Density interpolation function. */
	func_points_t fRHO;
	/** Specific heat interpolation function. */
	func_points_t fCP;
	/** emissivity interpolation function. */
	func_points_t fEPS;
	/** Destruction temperature interpolation function. */
	func_points_t fTdestr;
	/** Ablation A-coefficient interpolation function. */
	func_points_t fA;
	/** Ablation B-coefficient interpolation function. */
	func_points_t fB;
public:
	/** Class constructor
	 * @param title - material short name
	 * @param flog - pointer to the output device
	 * @param L - heat conductivity
	 * @param RHO - density
	 * @param CP - specific heat
	 * @param EPS - emissivity
	 * @param Tdestr - destruction temperature
	 * @param A - ablation coefficient
	 * @param B - ablation coefficient
	 * @param Type - boundary type.
	 */
	CUserMaterial(const char* title, FILE* flog, double L, double RHO, double CP, double EPS,
			double Tdestr, double A, double B, int Type);
	/**
	 * Set interpolation function to the property
	 * @param prm - material property
	 * @param table - interpolation fuction for the property
	 */
	int update(const char* prm, const func_points_t *table);
	/** Heat Conductivity.
	 * @param T - current temperature.
	 */
	virtual double l(double T);
	/** Heat Specific heat.
	 * @param T - current temperature.
	 */
	virtual double c(double T);
	/** density.
	 * @param T - current temperature.
	 */
	virtual double r(double T);
	/** Emissivity.
	 * @param T - current temperature.
	 */
	virtual double eps(double T);
	/** Ablation A-coefficient.
	 * @param T - current temperature.
	 */
	virtual double a(double T);
	/** Ablation B-coefficient.
	 * @param T - current temperature.
	 */
	virtual double b(double T);
	virtual double Td(double T);
	virtual int at(void);
	/**
	 * @brief Heat Quantity at the reference temperature
	 * @param T - reference temperature
	 * @return Heat Quantity, [J*m^-3]
	 */
	virtual double heatQuantity(double T);
	/**
	 * @brief Check correctess of material thermophysics
	 * SUCCESS - if all defined correctly; NULL_VALUE, overwise
	 */
	virtual int check();
	/** Print class state.
	 * @param out - pointer to the output device.
	 */
	virtual void print(FILE *out = stdout);
};
#endif /* _MATERIAL_H_ */
