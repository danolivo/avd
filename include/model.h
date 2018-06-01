/**
 * @brief Структура расчётной модели
 * @author А.В. Лепихов
 * @date 2013 г.
 * @copyright MIT License
 */

#ifndef _MODEL_H_
#define _MODEL_H_

#include <common.h>
#include <material.h>
/** Прототип граничного условия. */
class CBoundary;

/**
 * @brief Тепловая модель пакета материалов
 */
class thm_t {
public:
	/** Pointer to the left boundary */
	CBoundary *LBC;
	/** Pointer to the right boundary */
	CBoundary *RBC;
	/** Pointer to material at each cell. */
	CMaterial *m[CELLS_MAX_NUM];
	/** Cell temperature */
	double T[CELLS_MAX_NUM];
	/** Cell width */
	double width[CELLS_MAX_NUM];
	/** First cell number at the Model (connected to the left boundary) */
	int fcnum;
	/** Last cell number at the Model (connected to the left boundary) */
	int lcnum;
	/** left boundary movement distance during solution */
	double LDEL;
	/** Using for limiting minimum cell size during ablation process. */
	double PrimaryLeftCellSize;
	/** Using for limiting minimum cell size during ablation process. */
	double PrimaryRightCellSize;
	/** Начальный шаг счёта. */
	double INIT_TIMESTEP;
	/** Текущее расчетное время. */
	double CURRENT_TIME;
	/** Температура на левой границе, К. */
	double TWL;
	/** Температура на правой границе, К. */
	double TWR;
	/**
	 * @brief Class constructor.
	 */
	thm_t();
	/**
	 * @brief Add new layer to the right side of TPS
	 * @param Cells - Number of cells in the TPS layer
	 * @param Width - Width of the TPS layer
	 * @param Material - Pointer to material at each cell of TPS layer
	 * @param Temp - Temperature at each cell of the TPS layer
	 * @return Operation result code
	 */
	int add(int Cells, double Width, CMaterial* Material, double Temp);
	/**
	 * @brief Crop a small piece of model.
	 * @param side - side of model. 0 - left side.
	 * @param dx - value of deleting width, [m].
	 */
	int crop(int side, double dx);
	/** Returns length of the model. */
	double length();
	/** Set left boundary.
	 * @param BC - pointer tj the boundary.
	 */
	void setLBC(CBoundary* BC);
	/** Set right boundary.
	 * @param BC - pointer tj the boundary.
	 */
	void setRBC(CBoundary* BC);
	/** Set non-uniform tempareture by interpolation function
	 * @param Tf - pointer to the temperature-by-coordinate interpolation function
	 */
	void setTemperature(func_points_t * Tf);
	void print();
};

#endif /* _MODEL_H_ */
