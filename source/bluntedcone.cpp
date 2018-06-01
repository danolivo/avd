#include <assert.h>
#include <cstdlib>
#include <cmath>
#include <bluntedcone.h>
#include <cstdio>

#define PI				(3.1415926)

CBluntedCone::CBluntedCone(double R0, double th) {
	this->R0 = R0;
	this->th = PI*th/180.; // В радианах
	this->w = PI/2. - this->th; // Центральный угол.
}

double CBluntedCone::x0() const {
	return R0*(1. - cos(w));
}

double CBluntedCone::y0() const {
	return R0*sin(w);
}

double CBluntedCone::mr(double x) const {
	if (x<=x0()) {
		double a;
		a = acos((R0-x)/R0);
		return R0*sin(a);
	} else {
		return y0() + (x-x0())*tan(th);
	}
}

double CBluntedCone::theta(double x) const {
	if (x >= x0())
		return th;
	else
		return PI/2. - acos((R0-x)/R0);
}


double CBluntedCone::x(double y) const {
	
	if (y <= y0())
		return R0-sqrt(pow(R0, 2.)-pow(y, 2.));
	else
		return x0() + (y-y0())/tan(th);
}

double CBluntedCone::xgen(double x) const {
	double angle;
	if (x <= x0()) {
		angle = acos((R0-x)/R0);
		return R0*angle;
	} else {
		angle = acos((R0-x0())/R0);
		return R0*angle + (x-x0())/cos(th);
	}
}

double CBluntedCone::Sbottom(double x) const {

	return PI*pow(mr(x), 2.);
}

double CBluntedCone::R()
{
	return R0;
}
