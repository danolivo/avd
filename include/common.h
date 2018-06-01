/**
 * @brief Общие типы и структуры данных
 * @author А.В. Лепихов
 * @date 2013 г.
 * @copyright MIT License
 */

#ifndef _COMMON_H_
#define _COMMON_H_

#include <assert.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <math_addons.h>

/** Пустое значение. */
#define NULL_VALUE		(-1)
/** Успешное выполнение работы. */
#define SUCCESS			(0)
/** Максимальное количество ячеек в расчетной модели. */
#define CELLS_MAX_NUM		(1000)
/** Максимальное количество слоев материалов в расчетной модели. */
#define LAYERS_MAX_NUM		(20)
/** Максимальная длина строки, содержащей имя файла. */
#define FILENAME_MAX_LEN	(256)
/** Максимальная длина строки во входном файле. */
#define STRING_MAX_LEN		(1024)
/** Максимальное количество опорных точек траектории. */
#define TRAJECTORY_MAX_LEN	(1000)
/** Максимальное количество материалов. */
#define MATERIALS_MAX_NUM	(20)
/** Стандартный минимальный шаг счёта, если в файле ИД не указано иное. */
#define STD_TIMESTEP_MIN	(1.0E-4)
/** Стандартный максимальный шаг счёта, если в файле ИД не указано иное. */
#define STD_TIMESTEP_MAX	(0.01)

#endif /* _COMMON_H_ */

