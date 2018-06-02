#include <cstdio>
#include <cstring>
#include <assert.h>
#include <boundary.h>

CBoundary::CBoundary(int Type, FILE* flog)
{

	this->Type = Type;
	this->flog = flog;
}
CBoundary::CBoundary(CBoundary *boundary)
{
	this->Type = boundary->type();
	this->flog = boundary->flog;
	BC_t BCtemp = boundary->GetBC(0.);
	BC.eps = BCtemp.eps;
	BC.Tw = BCtemp.Tw;
	BC.Ie = BCtemp.Ie;
	BC.Iw = BCtemp.Iw;
	BC.acp = BCtemp.acp;
}
int CBoundary::type()
{
	return this->Type;
}

CFOBoundary::CFOBoundary(double Tw, FILE* flog) : CBoundary(1, flog)
{
	assert(Tw >= 0);
	BC.eps = 0;
	BC.I0 = -1;
	BC.Ie = -1;
	BC.Iw = -1;
	BC.P = -1;
	BC.acp = 1e12;
	BC.Tw = Tw;
	this->flog = flog;
}
BC_t CFOBoundary::GetBC(double Time)
{
	return BC;
}

void CFOBoundary::print()
{
	printf("Tw=%lf\n", BC.Tw);
	return;
}
int CFOBoundary::SetTimeChangeCoefs(const char* pname, const func_points_t* table)
{
	assert(0);
	return NULL_VALUE;
}

CSOBoundary::CSOBoundary(double acp, double I0, double Ie, double Iw,  double Ps, double eps, FILE* flog) : CBoundary(2, flog)
{
	assert(Ie >= 0.);
//	assert(acp >= 0.);
	assert(eps >= 0.);
	BC.eps = eps;
	BC.I0 = I0;
	BC.Ie = Ie;
	BC.Iw = Iw;
	BC.P = Ps;
	BC.acp = acp;
	BC.Tw = -1.;
	this->flog = flog;
}
BC_t CSOBoundary::GetBC(double Time) {

	return BC;
}

void CSOBoundary::print() {
	printf("\nEPS, IE, ACP\n");
	printf("eps=%5.3lf\tIe=%lf\tacp=%lf\n", BC.eps, BC.Ie, BC.acp);
	return;
}
int CSOBoundary::SetTimeChangeCoefs(const char* pname, const func_points_t* table)
{
	assert(0);
	return NULL_VALUE;
}
