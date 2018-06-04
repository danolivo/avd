/**
 * @file atmosphere.h
 * @brief Параметры стандартной атмосферы
 * @details Реализация ГОСТ 4401-81 с изменениями от 2004 г.
 * @author Разработчик А.В. Лепихов
 * @date 2011 г.
 */

#ifndef _ATMOSPHERE_H
#define _ATMOSPHERE_H

/*
 * H - геометрическая высота над уровнем моря, [м]
*/

/* Параметры атмосферы на среднем уровне моря */
/** скорость звука [м/с] */
extern const double ac;
/** ускорение свободного падения [м/с^2] */
extern const double gc;
/** масштаб высоты по давлению [м] */
extern const double Hpc;
/** средняя длина свободного пробега частиц воздуха [м] */
extern const double lc;
/** молярная масса воздуха [кг/кмоль] */
extern const double Mc;
/** концентрация частиц [m^-3] */
extern const double nc;
/** давление [Па] */
extern const double Pc;
/** температура Кельвина [К] */
extern const double Tc;
/** средняя скорость частиц воздуха [м/с] */
extern const double Vc;
/** удельный вес [H/м^3] */
extern const double Gammac;
/** кинематическая вязкость [м^2/c] */
extern const double nuc;
/** динамическая вязкость [Па*с] */
extern const double muc;
/** теплопроводность [Вт/(м*К)] */
extern const double lambdac;
/** частота соударений частиц воздуха [с^-1] */
extern const double omegac;
/** плотность [кг/м^3] */
extern const double rhoc;
/* ----------------------------------------- */

/* Базовые константы */
/** число Авогадро, кмоль^-1 */
extern const double Na;
/** универсальная газовая постоянная [Дж*К^-1*кмоль^-1] */
extern const double Rstar;
/** удельная газовая постоянная [Дж*кг^-1*K^-1] */
extern const double R;
/** условный радиус Земли [м] */
extern const double r;
/** показатель адиабаты для воздуха */
extern const double hi;
/** эффективный диаметр молекул воздуха при столкновении [м] */
extern const double a;
/**
 * Эмпирический коэффициент Сатерленда в уравнении для определения динамической
 * вязкости
*/
extern const double S; // [K]
/**
 * Эмпирический коэффициент Сатерленда в уравнении для определения динамической
 * вязкости
*/
extern const double BethaS; // [кг*с^-1*м^-1*K^-0.5]
/* ----------------- */

/* Функции */

/**
 * @brief Молярная масса воздуха М [кг/кмоль]
 * @details Высоты: от минус 2000 м до 1 200 000 м.
 * @param H - высота над уровнем моря.
 */
double sa_MolarMass(double H);

/**
 * @brief Температура T [град К]
 * @details Высоты: от минус 2000 м до 1 200 000 м.
 * @param H - высота над уровнем моря.
 */
double sa_Temperature(double H);

/**
 * @brief Молярная температура T [град К]
 * @details Высоты: от минус 2000 м до 1 200 000 м.
 * @param H - высота над уровнем моря.
 */
double sa_MolarTemperature(double H);

/**
 * @brief Концентрация частиц воздуха, n
 * @details Высоты: от минус 2000 м до 1 200 000 м.
 * @param H - высота над уровнем моря.
 */
double sa_AirParticlesConcentration(double H);

/**
 * @brief Статическое давление [Па]
 * @details Высоты: от минус 2000 м до 1 200 000 м.
 * @param H - высота над уровнем моря.
 */
double sa_StaticPressure(double H);

/**
 * @brief Плотность воздуха [кг/м^3]
 * @details Высоты: от минус 2000 м до 1 200 000 м.
 * @param H - высота над уровнем моря.
 */
double sa_Density(double H);

/**
 * @brief Скорость звука [м/с]
 * @details Высоты: от минус 2000 м до 1 200 000 м.
 * @param H - высота над уровнем моря.
 */
double sa_SoundSpeed(double H);

/**
 * @brief Динамическая вязкость [Па*с]
 * @details Высоты: от минус 2000 м до 90 000 м.
 * @param H - высота над уровнем моря.
 */
double sa_DynamicViscosity(double H);

/**
 * @brief Кинематическая вязкость [м^2*с]
 * @details Высоты: от минус 2000 м до 90 000 м.
 * @param H - высота над уровнем моря.
 */
double sa_KinematicViscosity(double H);

/**
 * @brief Коэффициент теплопроводности [Вт/(м*К)]
 * @details Высоты: от минус 2000 м до 1 200 000 м.
 * @param H - высота над уровнем моря.
 */
double sa_ThermalConductivity(double H);

/**
 * @brief Геопотенциальная высота [м]
 * @details Высоты: от минус 2000 м до 1 200 000 м.
 * @param H - высота над уровнем моря.
 */
double sa_GeopotentialHeight(double H);

/**
 * @brief Ускорение свободного падения [м/с^2]
 * @details Высоты: от минус 2000 м до 1 200 000 м.
 * @param H - высота над уровнем моря.
 */
double sa_g(double H);

/**
 * @brief Удельный вес [Н/м^-3]
 * @details Высоты: от минус 2000 м до 1 200 000 м.
 * @param H - высота над уровнем моря.
 */
double sa_SpecificGravity(double H);

/**
 * @brief Масштаб высоты (шкала высоты) по давлению [м]
 * @details Высоты: от минус 2000 м до 1 200 000 м.
 * @param H - высота над уровнем моря.
 */
double sa_Hp(double H);

/**
 * @brief Средняя скорость частиц воздуха [м/с]
 * @details Высоты: от минус 2000 м до 1 200 000 м.
 * @param H - высота над уровнем моря.
 */
double sa_AverageAirParticlesVelocity(double H);

/**
 * @brief Средняя длина свободного пробега частиц воздуха [м]
 * @details Высоты: от минус 2000 м до 1 200 000 м.
 * @param H - высота над уровнем моря.
 */
double sa_AverageAirParticlesPassing(double H);

/**
 * @brief Средняя частота соударений частиц воздуха [1/с]
 * @details Высоты: от минус 2000 м до 1 200 000 м.
 * @param H - высота над уровнем моря.
 */
double sa_AirPartsShockFrequency(double H);

#endif /* _ATMOSPHERE_H */
