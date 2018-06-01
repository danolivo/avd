#include <material.h>
#include <cstring>

CMaterial::CMaterial(const char* title, FILE* flog)
{
	strcpy(this->Name, title);
	this->flog = flog;
}

char* CMaterial::name()
{
	return this->Name;
}

ConstMaterial::ConstMaterial(const char* title, FILE* flog, double L, double RHO, double CP, double EPS,
		double Tdestr, double A, double B, int Type) : CMaterial(title, flog)
{
	this->L = L;
	this->RHO = RHO;
	this->CP = CP;
	this->EPS = EPS;
	this->A = A;
	this->B = B;
}

double ConstMaterial::l(double T)
{
	return L;
}
double ConstMaterial::c(double T)
{
	return CP;
}
double ConstMaterial::r(double T)
{
	return RHO;
}
double ConstMaterial::eps(double T)
{
	return EPS;
}
double ConstMaterial::a(double T)
{
	return A;
}
double ConstMaterial::b(double T)
{
	return B;
}
double ConstMaterial::heatQuantity(double T)
{
	return RHO*CP*T;
}
int ConstMaterial::check()
{
	if (L <= 0.)
		return NULL_VALUE;
	if (RHO <= 0.)
		return NULL_VALUE;
	if (CP <= 0.)
		return NULL_VALUE;
	if (EPS < 0.)
		return NULL_VALUE;
	if (Tdestr <= 0.)
		return NULL_VALUE;
	if (A <= 0.)
		return NULL_VALUE;
	if (B <= 0.)
		return NULL_VALUE;
	if (type < 0)
		return NULL_VALUE;
	return SUCCESS;
}
void ConstMaterial::print(FILE *out)
{

}
CUserMaterial::CUserMaterial(  const char* title, FILE* flog, double L, double RHO, double CP, double EPS,
				double Tdestr, double A, double B, int Type)
:ConstMaterial(title, flog, L, RHO, CP, EPS, Tdestr, A, B, Type)
{
	fL.count = 0;
	fCP.count = 0;
	fRHO.count = 0;
	fEPS.count = 0;
	fTdestr.count = 0;
	fA.count = 0;
	fB.count = 0;
}
int CUserMaterial::update(const char* prm, const func_points_t *table)
{
	if (strcmp(prm, "L") == 0) {
		func_points_cpy(table, &fL);
	} else
	if (strcmp(prm, "CP") == 0) {
		func_points_cpy(table, &fCP);
	} else
	if (strcmp(prm, "RHO") == 0) {
		func_points_cpy(table, &fRHO);
	} else
	if (strcmp(prm, "EPS") == 0) {
		func_points_cpy(table, &fEPS);
	} else
		return NULL_VALUE;
	return SUCCESS;
}
int CUserMaterial::check()
{
	if ((ConstMaterial::l(0.) <= 0.) && (fL.count == 0))
		return NULL_VALUE;
	if ((ConstMaterial::c(0.) <= 0.) && (fCP.count == 0))
			return NULL_VALUE;
	if ((ConstMaterial::r(0.) <= 0.) && (fRHO.count == 0))
			return NULL_VALUE;
	if ((ConstMaterial::eps(0.) <= 0.) && (fEPS.count == 0))
		return NULL_VALUE;
	return SUCCESS;
}
double CUserMaterial::l(double T)
{
	if (fL.count != 0)
		return in_LinearFunc(&fL, T);
	else
		return ConstMaterial::l(T);
}
double CUserMaterial::c(double T)
{
	if (fCP.count != 0)
		return in_LinearFunc(&fCP, T);
	else
		return ConstMaterial::c(T);

}
double CUserMaterial::r(double T)
{
	if (fRHO.count != 0)
		return in_LinearFunc(&fRHO, T);
	else
		return ConstMaterial::r(T);

}
double CUserMaterial::eps(double T)
{
		if (fEPS.count != 0)
		return in_LinearFunc(&fEPS, T);
	else
		return ConstMaterial::eps(T);
}
double CUserMaterial::a(double T)
{
	if (fA.count != 0)
		return in_LinearFunc(&fA, T);
	else
		return ConstMaterial::a(T);
}
double CUserMaterial::b(double T)
{
	if (fB.count != 0)
		return in_LinearFunc(&fB, T);
	else
	return ConstMaterial::b(T);
}
double CUserMaterial::heatQuantity(double T)
{
	if (fCP.count != 0)
		return r(T)*ma_integral(&fCP, T);
	else
		return ConstMaterial::heatQuantity(T);
}
void CUserMaterial::print(FILE *out)
{

}
