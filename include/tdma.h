#ifndef _TDMA_H_
#define _TDMA_H_
/**
 * –еализует метод прогонки.
 * @param Temp - массив температур. 
 * @param Width - массив длин €чеек.
 * @param VolHeatSrc -массив источников тепловыделени€.
 * @param HeatCond - массив теплопроводностей.
 * @param Density - массив плотностей.
 * @param SpecHeat - массив теплоемкостей. 
 * @param eps - степень черноты правого и левого √”.
 * @param thau - расчЄтное врем€.
 * @param size - количество €чеек в модели.
 */
extern void CalculateTDMA(	double* Temp, double* Width, double* VolHeatSrc,
					double* HeatCond, double* Density, double* SpecHeat,
					double Eps[2], double thau, int size);

#endif /* _TDMA_H_ */
