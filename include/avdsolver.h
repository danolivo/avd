/**
 * @brief �������� �������� ������ � ������, ����������� �� ������ ������ ��������� ����������.
 * @author �.�. �������
 * @date 2013 �.
 * @copyright MIT License
 */
#ifndef _AVDSOLVER_H_
#define _AVDSOLVER_H_

#include <common.h>
#include <avdtparser.h> // ������ ����� ��
#include <thsolver.h> // �������� ��������
#include <model.h> // �������� ������
#include <bluntedcone.h>


typedef struct {
	double I0, IW, IE, ISTAR, ALC, ALC1, QCONV, XAP1, P0, P1, V0, F1, F2, F3, KDIS, KENTH;
} avd_t;

/** ������������ ������ ��������� �������. */
typedef struct {
	/** ������ �������, ��� �������� ������ ������. */
	double time;
	/** ������, � */
	double H;
	/** ��������, �/�*/
	double V;
	/** ���� �����, ����. */
	double al;
	/** ���� ���������, ����. */
	double phi;
	/** ����� ���� */
	double mach;
	/** ����������� ����� ^0.2 */
	double XEF;
	/** ����������� �������� P/P0. */
	double PP0;
	avd_t avd;
	/** ����������� �����������. */
//	double acp;
	/** ������������ ���������, ��/��. */
//	double Ie;
	/** ��������� �� ������, ��/��. */
//	double Iw;
	/** ������ ���������, ��/��. */
//	double I0;
	/** �������� �� ������, ��. */
//	double Ps;
	/** ������������ �������� �����. */
	double G;
	/** ��������� ������� �������������. */
	solve_result_t srt;
} avdsolver_t;
/** �������� ���������� �������� ������ � ������, ����������� �� ������ ������ ��������� ����������. */
class AVDSolver {
	/** ��������� �� ���������������� ��������� � ����������� �����. */
	gasdynamics_t* gd;
	/** ��������� �� ��������� � �����������. */
	trm_t* trm;
	/** ��������� �� �������� ������. */
	thm_t* thm;
	/** ��������� �� ����� ��������� ��������. */
	CTHSolver* thsolver;
	/** ��� �����. */
	double TIMESTEP;

	CBluntedCone* BCone;
public:
	/** ����������� ������. */
	AVDSolver(thm_t* thm, trm_t* trm, gasdynamics_t* gd, CBluntedCone* BCone);
	/**
	 * @brief ��������� ������.
	 * @param time - ������ �������, �� �������� ���������� ��������� ������.
	 */
	avdsolver_t Solve(double time);
	/** ���������� ������. */
	~AVDSolver();
	/** ����� ��������� ������. */
	void print();
};

#endif /* _AVDSOLVER_H_ */
