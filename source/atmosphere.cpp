/*
 * Параметры стандартной атмосферы
 * Реализация ГОСТ 4401-81 с изменениями от 2004 г.
 * Разработчик А.В. Лепихов
*/

/*
 * Описание.
 * Настоящий стандарт устанавливает числовые значения основных параметров
 * атмосферы для высот от минус 2000 до 120 000 м. Средние значения атмосферы
 * устанавливаются для широты 45 град. 32' 33'', соответствующие среднему уровню
 * солнечной активности.
 * Стандарт предназначен для использования при расчётах и проектировании
 * летательных аппаратов, при обработке геофизических и метеорологических
 * наблюдений и для приведения результатов испытаний ЛА и их элементов к
 * одинаковым условиям и не предназначен для баллистических расчётов
 * искусственных спутников земли.
*/

/*
 * В качестве входных данных для вычисления того или иного параметра стандартной
 * атмосферы задается геометрическая высота H
 */

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <atmosphere.h>

/* Параметры атмосферы на среднем уровне моря */
const double ac		= 340.294; // скорость звука [м/с]
const double gc		= 9.80665; // ускорение свободного падения [м/с^2]
const double Hpc	= 8434.5; // масштаб высоты по давлению [м]
const double lc		= 66.328e-9; // средняя длина свободного пробега частиц воздуха [м]
const double Mc		= 28.96442; // молярная масса воздуха [кг/кмоль]
const double nc		= 25.471e24; // концентрация частиц [m^-3]
const double Pc		= 101325; // давление [Па]
const double Tc		= 288.15; // температура Кельвина [К]
const double Vc		= 458.94; // средняя скорость частиц воздуха [м/с]
const double Gammac	= 12.013; // удельный вес [H/м^3]
const double nuc	= 1.4607e-5; // кинематическая вязкость [м^2/c]
const double muc	= 1.7e-5; // динамическая вязкость [Па*с]
const double lambdac = 0.025343; // теплопроводность [Вт/(м*К)]
const double omegac = 6.9193e9; // частота соударений частиц воздуха [с^-1]
const double rhoc	= 1.225; // плотность [кг/м^3]
/* ----------------------------------------- */

/* Базовые константы */
const double Na		= 602.257e24; // число Авогадро, кмоль^-1
const double Rstar	= 8314.32; // универсальная газовая постоянная [Дж*К^-1*кмоль^-1]
const double R		= 287.05287; // удельная газовая постоянная [Дж*кг^-1*K^-1]
const double r		= 6356767.; // условный радиус Земли [м]

/*
 * Эмпирические коэффициенты Сатерленда в уравнении для определения динамической
 * вязкости
*/
const double S		= 110.4; // [K]
const double BethaS = 1.458e-6; // [кг*с^-1*м^-1*K^-0.5]
/* -------------- */

const double hi		= 1.4; // показатель адиабаты для воздуха
const double a		= 0.365e-9; // эффективный диаметр молекул воздуха при столкновении [м]

/*
 * Молярная масса воздуха М [кг/кмоль]
 * Air Molar Mass
 * H - высота над уровнем моря, [м]
 * return: значение молярной массы воздуха для высот от минус 2000 м до 1200 км.
 * Описание:
 * Атмосфера Земли представляет собой смесь газов, водяного пара и некоторого
 * количества аэрозолей. В определённых условиях в составе воздуха меняется
 * концентрация водяного пара, углекислого газа, озона и некоторых других
 * составляющих, содержание которых в атмосфере незначительно. Больше других
 * подвержено изменению содержание водяного пара, концентрация которого у
 * поверхности Земли при высокой температуре может достигать 4%, а с увеличением
 * высоты и понижением температуры быстро падает. Состав сухого воздуха в
 * области гомосферы до высоты 90-95 км остается практически постоянным.
*/
double sa_MolarMass(double H) {
	const double B[4][6] = {
		{46.9083, 40.4668, 6.3770, 75.6896, 112.4838, 9.8970}, // B0
		{-29.71210e-5, -15.52722e-5, 6.25497e-5, -17.61243e-5, -30.68086e-5, -1.19732e-5}, // B1
		{12.08693e-10, 3.55735e-10, -1.10144e-10, 1.33603e-10, 2.90329e-10, 7.78247e-12}, // B2
		{-1.85675e-15, -3.02340e-16, 3.36907e-17, -2.87884e-17, -9.20516e-17, -1.77541e-18} // B3
	};
	
	assert((H >= -2000.) && (H <= 1200000.));
	
	if (H < 94000.)
		return Mc;
	if ((H > 94000.) && (H <= 97000.))
		return 28.82+0.158 * pow(1-7.5e-08*pow(H-94000,2), 0.5) - 2.479e-4 * pow(97000-H, 0.5);
		
	// Линейный закон
	if ((H > 97000.) && (H <= 97500.))
		return -0.00012*(H-97000.) + sa_MolarMass(97000.);
	if ((H > 97500.) && (H <= 120000.))
		return -0.0001511*(H-97500.) + sa_MolarMass(97500.);

	// Полиномы
	if ((H > 120000.) && (H <= 250000.))
		return B[0][0] + B[1][0]*H + B[2][0]*pow(H,2) + B[3][0]*pow(H,3);
		
	if ((H > 250000.) && (H <= 400000.))
		return B[0][1] + B[1][1]*H + B[2][1]*pow(H,2) + B[3][1]*pow(H,3);
		
	if ((H > 400000.) && (H <= 650000.))
		return B[0][2] + B[1][2]*H + B[2][2]*pow(H,2) + B[3][2]*pow(H,3);
		
	if ((H > 650000.) && (H <= 900000.))
		return B[0][3] + B[1][3]*H + B[2][3]*pow(H,2) + B[3][3]*pow(H,3);
		
	if ((H > 900000.) && (H <= 1050000.))
		return B[0][4] + B[1][4]*H + B[2][4]*pow(H,2) + B[3][4]*pow(H,3);
		
	if ((H > 1050000.) && (H <= 1200000.))
		return B[0][5] + B[1][5]*H + B[2][5]*pow(H,2) + B[3][5]*pow(H,3);

}

/*
 * Возвращает параметры, принятые в данном ГОСТе для расчета свойств
 * стандартной атмосферы на высотах от минус 2000 до 120 000 м.
 * H - геопотенциальная высота над уровнем моря [м]
 * out H0 - опорная величина геопотенциальной высоты для данного H
 * out Tm0 - опорная величина молярной температуры для данного H
 * out BethaM - опорная величина градиента молярной температуры для данного H
 *
 * Параметры "Молярная температура" и "градиент молярной температуры"
 * совпадают с обычной температурой и градиентом до высоты 94 000 м.
*/
static double standardAtmosphereRanges(double H, double* H0, double* Tm0, double* BethaM) {

	assert((H >= -2000.) && (H <= 1200000.));
	assert(H0 != 0);
	assert(Tm0 != 0);
	assert(BethaM != 0);
	
	if (H > 117777.) {
		*BethaM = 0.011;
		*H0 = 117777.;
		*Tm0 = 380.6;
	}
	
	if ((H >= -2000.) && (H <= 0.)) {
		*BethaM = -0.0065;
		*H0 = -2000.;
		*Tm0 = 301.15;
	}
	if ((H > 0.) && (H <= 11000.)) {
		*BethaM = -0.0065;
		*H0 = 0.;
		*Tm0 = 288.15;
	}
	if ((H > 11000.) && (H <= 20000.)) {
		*BethaM = 0.0;
		*H0 = 11000.;
		*Tm0 = 216.65;
	}
	if ((H > 20000.) && (H <= 32000.)) {
		*BethaM = 0.001;
		*H0 = 20000.;
		*Tm0 = 216.65;
	}
	if ((H > 32000.) && (H <= 47000.)) {
		*BethaM = 0.0028;
		*H0 = 32000.;
		*Tm0 = 228.65;
	}
	if ((H > 47000.) && (H <= 51000.)) {
		*BethaM = 0.0;
		*H0 = 47000.;
		*Tm0 = 270.65;
	}
	if ((H > 51000.) && (H <= 71000.)) {
		*BethaM = -0.0028;
		*H0 = 51000.;
		*Tm0 = 270.65;
	}
	if ((H > 71000.) && (H <= 85000.)) {
		*BethaM = -0.002;
		*H0 = 71000.;
		*Tm0 = 214.65;
	}
	if ((H > 85000.) && (H <= 94000.)) {
		*BethaM = 0.0;
		*H0 = 85000.;
		*Tm0 = 186.65;
	}
	if ((H > 94000.) && (H <= 102450.)) {
		*BethaM = 0.003;
		*H0 = 94000.;
		*Tm0 = 186.65;
	}
	if ((H > 102450.) && (H <= 117777.)) {
		*BethaM = 0.011;
		*H0 = 102450.;
		*Tm0 = 212.00;
	}
	return 0;
}

/*
 * Возвращает температуру T [град К]
 * Temperature
 * H - высота над уровнем моря, [м]
 * return: значение температуры для высот от минус 2000 м до 1200 км.
*/
double sa_Temperature(double H) {
	return sa_MolarTemperature(H)*sa_MolarMass(H)/Mc;
}

/*
 * Возвращает молярную температуру T [град К]
 * Molar Temperature
 * H - высота над уровнем моря, [м]
 * return: значение молярной температуры для высот от минус 2000 м до 1200 км.
*/
double sa_MolarTemperature(double H) {
	double Hgeo, T0, betha, h, h0;
	
	assert((H >= -2000.) && (H <= 1200000.));
	
	// Перевод геометрической высоты в геопотенциальную
	Hgeo = (6356767.*H)/(6356767.+H);
	
	// По геопотенциальной высоте до H=120 км
	if ((H >= -2000.) && (H <= 120000.)) {
		standardAtmosphereRanges(Hgeo, &h0, &T0, &betha);
		h = Hgeo;
//		printf("h=%f\t h0=%f\t T0=%f\t betha=%f\n", h, h0, T0, betha);
	}
	// Кусочно-линейное изменение температуры по геометрической высоте
	if ((H > 120000.) && (H <= 140000.)) {
		T0 = 334.42;
		betha = 0.011259;
		h0 = 120000.;
		h = H;
	}
	if ((H > 140000.) && (H <= 160000.)) {
		T0 = 559.6;
		betha = 0.0068;
		h0 = 140000.;
		h = H;
	}
	if ((H > 160000.) && (H <= 200000.)) {
		T0 = 695.6;
		betha = 0.00397;
		h0 = 160000.;
		h = H;
	}
	if ((H > 200000.) && (H <= 250000.)) {
		T0 = 834.4;
		betha = 0.00175;
		h0 = 200000.;
		h = H;
	}
	if ((H > 250000.) && (H <= 325000.)) {
		T0 = 941.9;
		betha = 0.00057;
		h0 = 250000.;
		h = H;
	}
	if ((H > 325000.) && (H <= 400000.)) {
		T0 = 984.65;
		betha = 0.00015;
		h0 = 325000.;
		h = H;
	}
	if ((H > 400000.) && (H <= 600000.)) {
		T0 = 995.9;
		betha = 0.00002;
		h0 = 400000.;
		h = H;
	}
	if ((H > 600000.) && (H <= 800000.)) {
		T0 = 999.9;
		betha = 0.0000005;
		h0 = 600000.;
		h = H;
	}
	if ((H > 800000.) && (H <= 1200000.)) {
		T0 = 1000.;
		betha = 0.0;
		h0 = 800000.;
		h = H;
	}

	return T0 + betha*(h-h0);
}

/*
 * Рекурсивное вычисление давления
 * Hgeo - геопотенциальная высота [м]
 * return - давление на геопотенциальной высоте Hgeo
*/
static double static_pressure_calc(double Hgeo) {
	double H0, Tm0, BethaM;
	
	// Высоты от минус 2000 до 120 000 м
	
	// База рекурсии
	if (Hgeo == -2000.)
		return 127774.;
	// Опорные параметры стандартной атмосферы по табл. 5
	standardAtmosphereRanges(Hgeo, &H0, &Tm0, &BethaM);
	
	if (BethaM !=0.)
		return pow(10., log10(static_pressure_calc(H0)) - (gc/(BethaM*R)) * log10((Tm0+BethaM*(Hgeo-H0))/Tm0));
	else
		return pow(10., log10(static_pressure_calc(H0)) - (0.434294*gc/(sa_Temperature((6356767.*Hgeo)/(6356767.-Hgeo))*R)) * (Hgeo-H0));
}

/*
 * Статическое давление [Па]
 * H - геометрическая высота [м]
 * return: статическое давление на высоте от минус 2000 до 1 200 000 м
*/
double sa_StaticPressure(double H) {

	assert((H>=-2000.) && (H<=1200000.));
	
	if (H <= 120000.)
		return static_pressure_calc((6356767.*H)/(6356767.+H));
	else {
		// Высоты от 120 000 м до 1 200 000 м
		return sa_AirParticlesConcentration(H)*Rstar*sa_Temperature(H)/Na;
	}
}

/*
 * Концентрация частиц воздуха, n
 * Air Particles Concentration
 * H - геометрическая высота [м]
*/
double sa_AirParticlesConcentration(double H) {
	const double A[5][9] = {
		{0.210005867e4,		0.10163937e4,		0.7631575e3,		0.1882203e3,		0.2804823e3,		0.5599362e3,		0.8358756e3,		0.8364965e2,			0.383220e2	},
		{-0.561844757e-1,	-0.2119530830e-1,	-0.1150600844e-1,	-0.2265999519e-2,	-0.2432231125e-2,	-0.3714141392e-2,	-0.4265393073e-2,	-0.3162492458e-3,	-0.50980e-4	},
		{0.5663986231e-6,	0.1671627815e-6,	0.6612598428e-7,	0.1041726141e-7,	0.8055024663e-8,	0.9358870345e-8,	0.8252842085e-8,	0.4602064246e-9,	0.18100e-10	},
		{-0.2547466858e-11,	-0.5894237068e-12,	-0.1708736137e-12,	-0.2155574922e-13,	-0.1202418519e-13,	-0.1058591881e-13,	-0.7150127437e-14,	-0.3021858469e-15,	0.0			},
		{0.4309844119e-17,	0.7826684089e-18,	0.1669823114e-18,	0.1687430962e-19,	0.6805101379e-20,	0.4525531532e-20,	0.2335744331e-20,	0.7512304301e-22,	0.0			}
	};
	const double m[9] = {17, 16, 15, 15, 14, 13, 12, 12, 11};
	int range;
	
	assert((H >= -2000.) && (H <= 1200000.));
	
	if (H < 120000.)
	
		return 7.243611e22*sa_StaticPressure(H)/sa_Temperature(H);
		
	else {
		if ((H >= 120000.) && (H <=150000.))
			range = 0;
		if ((H > 150000.) && (H <=200000.))
			range = 1;
		if ((H > 200000.) && (H <=250000.))
			range = 2;
		if ((H > 250000.) && (H <=350000.))
			range = 3;
		if ((H > 350000.) && (H <=450000.))
			range = 4;
		if ((H > 450000.) && (H <=600000.))
			range = 5;
		if ((H > 600000.) && (H <=800000.))
			range = 6;
		if ((H > 800000.) && (H <=1000000.))
			range = 7;
		if ((H > 1000000.) && (H <=1200000.))
			range = 8;
		return (A[0][range] + A[1][range]*H + A[2][range]*pow(H,2) + A[3][range]*pow(H,3) + A[4][range]*pow(H,4))*pow(10.,m[range]);
	}
}

/*
 * Плотность воздуха [кг/м^3]
 * H - геометрическая высота над уровнем моря [м]
 * return: плотность воздуха на высотах от минус 2000 до 1 200 000 м
*/
double sa_Density(double H) {
	
	return (sa_StaticPressure(H) * sa_MolarMass(H))/(sa_Temperature(H)*Rstar);
}

/*
 * Масштаб высоты (шкала высоты) по давлению [м]
 * H - геометрическая высота над уровнем моря [м]
 * return: масштаб высоты на высотах от минус 2000 до 1 200 000 м
 * Описание: высота однородной атмосферы
*/
double sa_Hp(double H) {
	return Rstar*sa_Temperature(H)/(sa_MolarMass(H)*sa_g(H));
}

/*
 * Возвращает скорость звука [м/с]
 * Sound Speed
 * H - высота над уровнем моря, [м]
 * return: значение скорости звука для высот от минус 2000 м до 1200 км.
*/
double sa_SoundSpeed(double H) {

	assert((H>=-2000.) && (H<=1200000.));
	
	return 20.046796*sqrt(sa_Temperature(H));
}

/*
 * Возвращает динамическую вязкость [Па*с]
 * Dynamic Viscosity
 * H - высота над уровнем моря, [м]
 * return: значение динамической вязкости по методу Сатерленда для высот
 * от минус 2000 м до 90 км.
*/
double sa_DynamicViscosity(double H) {

//	assert((H>=-2000) && (H<=90000));
	
	return BethaS*pow(sa_Temperature(H), 1.5)/(sa_Temperature(H)+S);
}

/*
 * Геопотенциальная высота [м]
 * H - геометрическая высота над уровнем моря [м]
 * return значение геопотенциальной высоты в диапазоне от минус 2000 до 1 200 000 м
*/
double sa_GeopotentialHeight(double H) {

	assert((H>=-2000.) && (H<=1200000.));
	
	return (6356767.*H)/(6356767.+H);
}

/*
 * Ускорение свободного падения [м/с^2]
 * H - геометрическая высота над уровнем моря [м]
 * return значение ускорения свободного падения в диапазоне от минус 2000 до 1 200 000 м

*/
double sa_g(double H) {

	assert((H>=-2000) && (H<=1200000.));
	
	return gc*pow(r/(r+H),2.);
}

/*
 * Удельный вес [Н/м^-3]
 * Specific Gravity
 * H - геометрическая высота над уровнем моря [м]
 * return значение удельного веса в диапазоне от минус 2000 до 1 200 000 м

*/
double sa_SpecificGravity(double H) {

	assert((H>=-2000.) && (H<=1200000.));
	
	return sa_Density(H)*sa_g(H);
}

/*
 * Кинематическая вязкость [м^2*с]

 * Kinematic Viscosity
 * H - высота над уровнем моря, [м]
 * return: значение кинематической вязкости по методу Сатерленда для высот
 * от минус 2000 м до 1200 км.
*/
double sa_KinematicViscosity(double H) {

//	assert((H>=-2000) && (H<=90000));
	
	return sa_DynamicViscosity(H)/sa_Density(H);
}

/*
 * Средняя скорость частиц воздуха [м/с]
 * H - высота над уровнем моря, [м]
 * return: значение средней скорости для высот от минус 2000 м до 1200 км.
*/
double sa_AverageAirParticlesVelocity(double H) {

	assert((H>=-2000.) && (H<=1200000.));
	
	return sa_AirPartsShockFrequency(H)*sa_AverageAirParticlesPassing(H);
}

/*
 * Коэффициент теплопроводности [Вт/(м*К)]

 * Thermal Conductivity
 * H - высота над уровнем моря, [м]
 * return: значение коэффициента теплопроводности для высот
 * от минус 2000 м до 1200 км.

*/
double sa_ThermalConductivity(double H) {
	assert((H>=-2000.) && (H<=1200000.));
	
	return (2.648151e-3 * pow(sa_Temperature(H), 1.5))/(sa_Temperature(H)+245.4*pow(0.1, 12/sa_Temperature(H)));
}

/*
 * Средняя длина свободного пробега частиц воздуха [м]
 * H - высота над уровнем моря, [м]
 * return: значение средней длины свободного пробега для высот от минус 2000 м до 1200 км.
*/
double sa_AverageAirParticlesPassing(double H) {

	assert((H>=-2000.) && (H<=1200000.));
	
	return 2.332376e-5*sa_Temperature(H)/sa_StaticPressure(H);
}

/*
 * Средняя частота соударений частиц воздуха [1/с]
 * H - высота над уровнем моря, [м]
 * return: значение средней частоты соударений для высот от минус 2000 м до 1200 км.
*/
double sa_AirPartsShockFrequency(double H) {

	assert((H>=-2000.) && (H<=1200000.));
	
	return 6.238629e6*sa_StaticPressure(H)/(sqrt(sa_Temperature(H)*sa_MolarMass(H)));
}

/*
 * Тест
*/
/*
int main(void) {
	// Высота 5000
	// Молярная масса
	printf("H=5 000\t M=%2.3f\n", sa_MolarMass(5000));
	// Высота 95 000, M=28.96
	printf("H=95 000\t M=%2.3f\n", sa_MolarMass(95000));
	// Высота 97 000, M=28.910
	printf("H=97 000\t M=%2.3f\n", sa_MolarMass(97000));
	// Высота 98 000, M=28.774
	printf("H=98 000\t M=%2.3f\n", sa_MolarMass(98000));

	// Высота 200 000, 20.98
	printf("H=200 000\t M=%2.3f\n", sa_MolarMass(200000));
	// Высота 300 000, M=17.74
	printf("H=300 000\t M=%2.3f\n", sa_MolarMass(300000));
	// Высота 500 000, M=14.33
	printf("H=500 000\t M=%2.3f\n", sa_MolarMass(500000));
	// Высота 700 000, M=7.99
	printf("H=700 000\t M=%2.3f\n", sa_MolarMass(700000));
	// Высота 950 000, M=4.11
	printf("H=950 000\t M=%2.3f\n", sa_MolarMass(950000));
	// Высота 1 195 500, M=3.67
	printf("H=1 195 000\t M=%2.3f\n", sa_MolarMass(1195000));
	// Температура
	// Высота -1000, T=294.651
	printf("H=-1 000\t T=%2.3f\n", sa_Temperature(-1000));
	// T=255.650
	printf("H=5 004\t T=%2.3f\n", sa_Temperature(5004));
	// T=216.650
	printf("H=15 035\t T=%2.3f\n", sa_Temperature(15035));
	// T=221.650
	printf("H=25 099\t T=%2.3f\n", sa_Temperature(25099));
	// T=250.350
	printf("H=40 000\t T=%2.3f\n", sa_Temperature(40000));
	// T=270.650
	printf("H=50 396\t T=%2.3f\n", sa_Temperature(50396));
	// T=233.292
	printf("H=65 000\t T=%2.3f\n", sa_Temperature(65000));
	// T=196.650
	printf("H=81 020\t T=%2.3f\n", sa_Temperature(81020));
	// a=282.538
	printf("H=81 020\t a=%2.3f\n", sa_SoundSpeed(80000));
	// mu=1.3208
	printf("H=81 020\t mu=%e\n", sa_DynamicViscosity(80000));
	// lambda=1.7987
	printf("H=81 020\t lambda=%e\n", sa_ThermalConductivity(80000));
	// static pressure=113931
	printf("H=-1 000\t Static Pressure=%e\n", sa_StaticPressure(-1000));
	// static pressure=54019.9
	printf("H=5 004\t Static Pressure=%e\n", sa_StaticPressure(5004));

	// static pressure=0.886272
	printf("H=81 020\t Static Pressure=%e\n", sa_StaticPressure(81020));
	printf("H=1000000\t Static Pressure=%e\n", sa_StaticPressure(1000000)); // =7.51043e-9
	// Концентрация частиц
	printf("H=50396\t Static Pressure=%e\n", sa_StaticPressure(50396)); // =7.59+1
	printf("H=50000\t Static Pressure=%e\n", sa_StaticPressure(50000));
	printf("H= 50396\t T=%2.3f\n", sa_Temperature(50396)); // =270.65
	printf("H= 50000\t T=%2.3f\n", sa_Temperature(50000));
	printf("H=50000\t n=%e\n", sa_AirParticlesConcentration(50000)); // =2.1352+22
	printf("H=130000\t n=%e\n", sa_AirParticlesConcentration(130000)); // =2.367+17
	printf("H=160000\t n=%e\n", sa_AirParticlesConcentration(160000)); // =3.162+16 
	printf("H=220000\t n=%e\n", sa_AirParticlesConcentration(220000)); // =4.037+15
	printf("H=300000\t n=%e\n", sa_AirParticlesConcentration(300000)); // =6.507+14
	printf("H=400000\t n=%e\n", sa_AirParticlesConcentration(400000)); // =1.057+14
	printf("H=500000\t n=%e\n", sa_AirParticlesConcentration(500000)); // =2.189+13
	printf("H=700000\t n=%e\n", sa_AirParticlesConcentration(700000)); // =2.312+12
	printf("H=900000\t n=%e\n", sa_AirParticlesConcentration(900000)); // =7.873+11
	printf("H=1100000\t n=%e\n", sa_AirParticlesConcentration(1100000)); // =4.145+11
	// Плотность
	printf("H=1100000\t rho=%e\n", sa_Density(1100000)); // =2.60170e-15
	printf("H=110000\t rho=%e\n", sa_Density(110000)); // =9.34035e-8
	printf("H=110000\t T=%e\n", sa_Temperature(110000)); // =255.487
	printf("H=110000\t P=%e\n", sa_StaticPressure(110000)); // =7.359e-3
	printf("H=75000\t nu=%e\n", sa_KinematicViscosity(75000)); // =3.4465e-1
	printf("H=71802\t Hgeo=%e\n", sa_GeopotentialHeight(71802)); // =71000
	printf("H=-1999\t Hgeo=%e\n", sa_GeopotentialHeight(-1999)); // =2000
	printf("H=71802\t g=%e\n", sa_g(71802)); // =9.5888
	printf("H=71000\t sg=%e\n", sa_SpecificGravity(71000)); // =6.9023e-4
	// Масштаб высоты
	printf("H=71000\t Hp=%e\n", sa_Hp(71000)); // =6489.9
	// Скорость частиц воздуха
	printf("H=71000\t Vavg=%e\n", sa_AverageAirParticlesVelocity(71000)); // =398.13
	printf("H=71000\t T=%e\n", sa_Temperature(71000)); // =216.846
	printf("H=71000\t M=%e\n", sa_MolarMass(71000)); // =Mc
	printf("H=80000\t Vavg=%e\n", sa_AverageAirParticlesVelocity(80000)); // =379.14
	printf("H=80000\t x=%e\n", 398.13*sa_MolarMass(80000)/sa_Temperature(80000));
	printf("H=71000\t x=%e\n", 379.14*sa_MolarMass(71000)/sa_Temperature(71000));
	printf("H=1000\t x=%e\n", 453.74*sa_MolarMass(1000)/sa_Temperature(1000));
	printf("H=1000\t Pass=%e\n", sa_AverageAirParticlesPassing(1000)); // =7.3090e-8
	printf("H=1000\t Shock=%e\n", sa_AirPartsShockFrequency(1000)); // =6.2079e9
	printf("H=80000\t Shock=%e\n", sa_AirPartsShockFrequency(80000)); // =8.6564e4
	printf("H=80000\t Vavg=%e\n", sa_AirPartsShockFrequency(80000)*sa_AverageAirParticlesPassing(80000)); // =

	return 0;
}
*/
