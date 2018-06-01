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
