
/**
 * @file magnitude.h
 * @brief Список используемых физических величин с указанием размерности.
 * 
 * @sa common.h
 * @version 0.1
 * @copyright MIT License
 * @author А.В. Лепихов
 * @date 2012
 */
 
#ifndef _MAGNITUDE_H
#define _MAGNITUDE_H

/**
 * @brief Количество Паскаль в одной технической атмосфере.
 */
#define AT	(98066.5)

/**
 * @brief Температура, [К]
 */
typedef double temperature_t;

/**
 * @brief Удельный тепловой поток, [Вт*м^-2]
 */
typedef double heatflux_t;

/**
 * @brief Удельный объёмный тепловой поток, [Вт*м^-3]
 */
typedef double vheatflux_t;

/**
 * @brief Энтальпия (теплосодержание), [Дж/кг]
 */
typedef double enthalpy_t; 

/**
 * @brief Коэффициент теплоотдачи, [Вт*м^-2*К^-1]
 */
typedef double alpha_t;

/**
 * @brief Угол, [град.]
 */
typedef double angle_t;

/**
 * @brief Угол, [рад.]
 */
typedef double radians_t;

/**
 * @brief Длина, [м]
 */
typedef double length_t;

/**
 * @brief Коэффициент теплопроводности материала, [Вт*м^-1*К^-1]
 */
typedef double heatconductivity_t;

/**
 * @brief Коэффициент теплоёмкости материала, [Дж*кг^-1*К^-1]
 */
typedef double specificheat_t;

/**
 * @brief Плотность материала, [кг*м^-3]
 */
typedef double density_t;

/**
 * @brief Масса материала, [кг]
 */
typedef double mass_t;

/**
 * @brief Скорость, [м*с^-1]
 */
typedef double velocity_t;

/**
 * @brief Давление, [Па]
 */
typedef double pressure_t;

#endif /* _MAGNITUDE_H */
