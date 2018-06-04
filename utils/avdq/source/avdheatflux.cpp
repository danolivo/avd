#include <assert.h>
#include <cmath>
#include <cstdio>

#include <magnitude.h>
#include <atmosphere.h>
#include <math_addons.h>
#include <avdheatflux.h>

//static double Istar(double hii);
static double omega(double hii);
static double Ie(double w);
static double calculate_hi(void);

static temperature_t Tw; /* Температура стенки */
static double cappa; /* Расчётный показатель адиабаты */
static bool regime; /* Режим течения: ламинарный или турбулентный */
static double P1divP0, total_pressure;
static velocity_t V;
static length_t H;
static double xeff;
static double XgenT;
static double theta;
static temperature_t Tdestr;
static double Rinf, Minf, Tinf;

static bool initSign=false; // Признак того, что функция инициализации подсистемы была вызвана как минимум 1 раз.

static double RHO = 0.0604;
static double MACH = 6.;
static double TEMP = 100.;

/**
 * @brief Давление торможения, Па.
 * @return Выдает давление торможения, рассчитанное по упрощенной методике, полученной из отчёта ИТМО СО РАН, либо константное
 * значение, заданное пользователем.
 */
double AVD_P0() {
	if (total_pressure < 0.) {
		double cp00 = 1.166+0.108*pow(2.5-1.5/pow(Minf, 2.), 2.);
		double pT0 = cp00/2. + 101325./(Rinf*pow(V, 2.));
		return pT0*Rinf*V*V;
	} else
		return total_pressure;
}

void AVD_Init(	length_t Height, velocity_t Vel, double Xeff, double PdivP0,
				double Twall, temperature_t Td, bool isTurbulent,
				length_t len, angle_t angle, pressure_t p0) {
	initSign = true;
	assert((PdivP0 <= 1.2) && (PdivP0 > 0.));
	assert(Twall > 0.);
	assert(Xeff > 0.);
	assert(Vel > 0.);
	assert(Vel*sa_SoundSpeed(Height) > 1.);
	assert((Height > -2000.) && (Height < 1.2E+06));
	if (Height < 0.) {
		Rinf = RHO;
		Tinf = TEMP;
		Minf = Vel/(20.04*sqrt(TEMP));
		printf("Spec. params: Rinf=%lf\tTinf=%lf\tMinf=%lf\n", RHO, TEMP, Minf);
	} else {
		Rinf = sa_Density(Height);
		Tinf = sa_Temperature(Height);
		Minf = Vel/sa_SoundSpeed(Height);
	}
	assert(len > 0.);
	initSign = true;
	theta = angle;
	XgenT = len;	
	Tdestr = Td;
	H = Height;
	V = Vel;
	xeff = Xeff;
	P1divP0 = PdivP0;
	total_pressure = p0;
	Tw = Twall;
	regime = isTurbulent;

	cappa = calculate_hi();
	AVD_Qw();
	
	return;
}

double AVD_P1() {
	return P1divP0*AVD_P0();
}

double AVD_Ie() {
	return Ie(omega(cappa));
}

static double x1_cap() {
	
	return (XgenT-1.57*(1.-theta/(90.)))*COS(theta); // здесь точно 90 град
}

static double x_cap_min() {

	return 3.*pow(24./(theta+4.), 3.-9./Minf);
}

/*
 * Функция w
*/
static double omega(double hii) {

	return 1./pow( P1divP0, (hii-1.)/hii ) - 1.;
}

static double bba_hi(double Istar) {

	assert(Istar >= 0.);

	if (Istar > 500.)
		return 1.23-95.6/Istar + 6.0E+04/pow(Istar, 2.);
	else
		return 1.4-0.12*(Istar/500.);
}

/* Итеративное вычисление показателя адиабаты */
static double calculate_hi(void) {
	double _hi, hii; // старое значение hi
	double I;
	
	hii = 1.23; // начальное значение в соответствии с методикой Авдуевского
/*	do {
		_hi = hii;
		I = Istar(_hi);
		hii = bba_hi(I);
		
	} while ((1.-hii/_hi) >= 0.01);
*/
	return hii;
}

static double Kdis() {
/*
	if (V>1000.)
		return 3.2*(pow(1.+(hi-1.)*Minf*Minf/2., 0.6)*pow(Istar(cappa), 0.05))/(1.375*pow(Minf, 1.2)*pow((hi+1.)/2., 0.6));
	else
		if (0.78*pow(Istar(cappa), 0.05) < 1.)
			return 1.;
		else
			return 0.78*pow(Istar(cappa), 0.05); */
	return 1;
}

static double Ksher() {
	if (Tw >= Tdestr)
		return 1.38+1.2*pow(P1divP0, 2.5);
	else
		return 1;
}

static double Kenth() {
/*	double K_Max_Enth;
	double K_Turb_Enth;
	
	K_Max_Enth = 1.11+0.055*(theta/10.);
		
	if (x1_cap() > x_cap_min())
		K_Turb_Enth = ((x1_cap()-x_cap_min())/(2.*x_cap_min())) * (K_Max_Enth-1.) + 1.;
	else
		K_Turb_Enth = 1.;
	
	if (K_Turb_Enth > K_Max_Enth) {
		return K_Max_Enth;
	} else
		return K_Turb_Enth; */
	return 1.;
}

static double Kdis_lam() {
/*	double K;
	if (Istar(cappa) > 1800.)
		K = 6.4/pow(Istar(cappa), 0.25);
	else
		// K = 3.42/pow(Istar(cappa), 0.49); // Было до 21.03.13
		K = 3.42/pow(Istar(cappa), 0.19);
	
	if (K > 1.)
		return 1.;
	else
		return K; */
	return 1;
}

/* Энтальпия торможения, [ккал/кг] */
double AVD_I0() {
	assert(initSign);
	
	return 0.24*Tinf + (V*V)/8370.;
}

/*
 * Значение энтальпии при температуре стенки, [ккал/кг]
 */
double AVD_Iw() {
	assert(initSign);
	assert(Tw > 0.);
	
	if (Tw <= 1000.)
		return 0.245*Tw;
	else if ((Tw > 1000.) && (Tw <= 2700. /*2500*/))
		return 245.+310.*(Tw-1000.)/1000.;
	else if (Tw < 4500.)
		return (Tw*Tw)/(8800.*pow(1.+0.05*log10(AVD_P1()), 2.));
	else 
		assert(0 || !printf("[EE]: Tw=%lf\n", Tw));
}

double AVD_Cappa() {
	return cappa;
}
/*
double AVD_Istar() {
	return Istar(cappa);
}
*/
/*double AVD_VdivVh() {
	double cap = calculate_hi();
	double VdivVh = sqrt(((2./(1.4-1.))*(1./(Minf*Minf))+1.)*(1.-pow(P1divP0, (1.4-1.)/1.4)));
	return VdivVh;
} */

/*
 * Характерная энтальпия, [ккал/кг]
 */
 /*
static double Istar(double hii) {
	double m, b, w;
	double V1divVNsharp;

	assert(initSign);
	assert(hii > 0.);
	assert((theta >= 0.) && (theta <= 30.));
	
	V1divVNsharp = 1.15*COS(theta)-0.055*pow(theta/10., 0.97)/Minf;
	if (V1divVNsharp < 0.)
		V1divVNsharp = 1.0;
	
	m = (((hi-1.)/2.) * Minf*Minf * V1divVNsharp*V1divVNsharp)/((1.-((hi-1.)/2.)*Minf*Minf*(1.-V1divVNsharp*V1divVNsharp))*(1./(1./pow(P1divP0, (hii-1.)/hii) - 1.)));
	assert(m >= 0.);
	b = ((x1_cap() - x_cap_min())/(2.*x_cap_min())) * (m-1.)+1.;
	if (x1_cap() < x_cap_min())
		if (fabs(b) > fabs(m))
			b = m;
		else
			b = 1.;

	w = omega(hii);
	assert(w >= 0.);
	
	if (b*w > ((1.-AVD_Iw()/AVD_I0())/(1.+AVD_Iw()/AVD_I0())))
		return AVD_I0() * 0.25 * ((b*w+1.)/(b*w)) * pow(1.-AVD_Iw()/Ie(w), 2.) + AVD_Iw();
	else
		return AVD_I0()/(1.+b*w);
} */

/*
 * Эффективная энтальпия потока в данном сечении тела вне погранслоя
 */
static double Ie(double w) {
	
	assert(initSign);
	
	return AVD_I0() * (1.+0.89*w)/(1.+w);
}

/*
 * Коэфициент теплоотдачи, [кг/(м^2*с)]
 */
double AVD_acp() {
	double K1, F1, F2, F3, c;
	
	assert(initSign);
	
	if (XgenT == 0.) {
		printf("You have'nt calculate A/Cp at critical point!\n");
		return 0.;
	}
	/* Параметр c используется в оригинальной методике Авдуевского
	 * Здесь, при нулевом угле атаки, он не нужен
	 */
	if (regime == TURBULENT) {
		c = 0.;
		K1 = (-1.+sqrt( 1.+(4.*(cappa-1.)/cappa)*c*c/(SIN(theta)*omega(cappa)) ))/3.;
		F1 = 1.+0.137*SIN((90.)*(1.-(Tw-273.)/1000.));
//		F2 = pow((hi+1.)/(hi-1.), 0.4) * pow(1.-pow(P1divP0, (hi-1.)/hi), 0.4) * pow(P1divP0, 0.8) * pow(AVD_I0()/Istar(cappa), 0.6);
		F3 = pow(pow((hi+1.)/2., (hi+1.)/(hi-1.))*pow(2./(hi-1.), 1./(hi-1.))*pow(Minf, 2.*hi/(hi-1.))/pow(2.*hi*Minf*Minf/(hi-1.)-1., 1./(hi-1.)), 0.8) * pow((1+(hi-1)*Minf*Minf/2)/(Minf*Minf*(hi+1)/2), 0.4) * (1/pow(hi*Minf*Minf, 0.2)) * (1/pow(1+Minf*Minf*(hi-1)/2, 0.6));
//		printf("ADC: %lf %lf %lf %lf %lf %lf %lf %lf\n", xeff, K1, F1, F2, F3, cappa, theta, SIN(theta));
		return 0.243e-2 * (pow(V, 1.2)*pow(Rinf/gc, 0.8)/xeff) * pow(1.+0.835*K1, 0.2)*F1*F2*F3*Kdis()*Kenth()*Ksher(); // [кг/(м^2*сек)
	} else {
/*		if (Istar(cappa) < 1800.)
			F1 = 1.45/pow(AVD_Iw(), 0.06);
		else*/
			F1 = 1.;
			
		F2 = pow((hi+1.)/(hi-1.), 0.25) * pow(1.-pow(P1divP0, (hi-1.)/hi), 0.25) * pow(P1divP0, 0.5);
		F3 = pow( pow((hi+1.)/2., (hi+1.)/(hi-1.)) * pow(2./(hi-1.), 1./(hi-1.)) * pow(Minf, 2.*hi/(hi-1.))/pow(2.*hi*Minf*Minf/(hi-1.)-1., 1./(hi-1.)), 0.5) * (1./Minf) * pow((1.+(hi-1.)*Minf*Minf/2.)/(Minf*Minf*(hi+1.)/2.), 0.25);
		
		return 0.0000108*(pow(V, 1.5)*pow(Rinf/gc, 0.5)/xeff)*F1*F2*F3*Kdis_lam()*Kenth();
		/* Для меня остался неясным вопрос: каким должно быть Kenth при 
		 ламинарном режиме. В определении расчёта ламинарного режима в док., в
		 заголовке прописано определение Kenth, однако в тексте определение
		 отсутствует. Использую то же, что и для турбулентного режима. */
	}

}

/*
 * Расчёт конвективного теплового потока к поверхности обтекаемого тела, [ккал/(м^2*c)]
 */
double AVD_Qw() {

	assert(initSign);
	/* Итеративно рассчитать показатель адиабаты */
	cappa = calculate_hi();

	double adc = AVD_acp();

	assert((adc > 0.)||!printf("AVD_Ie=%lf AVD_Iw=%lf I0=%lf V=%lf H=%lf\n", AVD_Ie(), AVD_Iw(), AVD_I0(), V, H));
//	printf("AVD_Ie=%lf\tIw=%lf\n",AVD_Ie(omega(cappa)), AVD_Iw());
	return adc*(Ie(omega(cappa)) - AVD_Iw());
}
