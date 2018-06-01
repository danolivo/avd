#ifndef _TDMA_H_
#define _TDMA_H_
/**
 * ��������� ����� ��������.
 * @param Temp - ������ ����������. 
 * @param Width - ������ ���� �����.
 * @param VolHeatSrc -������ ���������� ��������������.
 * @param HeatCond - ������ �����������������.
 * @param Density - ������ ����������.
 * @param SpecHeat - ������ �������������. 
 * @param eps - ������� ������� ������� � ������ ��.
 * @param thau - ��������� �����.
 * @param size - ���������� ����� � ������.
 */
extern void CalculateTDMA(	double* Temp, double* Width, double* VolHeatSrc,
					double* HeatCond, double* Density, double* SpecHeat,
					double Eps[2], double thau, int size);

#endif /* _TDMA_H_ */
