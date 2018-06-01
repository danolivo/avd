#ifndef BLUNTEDCONE_H
#define BLUNTEDCONE_H

/*
 * Геометрические параметры затупленного конуса
 */
class CBluntedCone {
	/**  радиус сферического затупления. */
	double R0;
	/** угол полураствора конуса. */
	double th;
	/** центральный угол. */
	double w;

public:
	/** Координата x стыка притупления и конуса. */
	double x0() const;
	/** Координата y стыка притупления и конуса. */
	double y0() const;
	/**
	 * Конструктор класса.
	 * @param R0 -радиус притупления.
	 * @param Th - угол полураствора конуса, град
	 */
	CBluntedCone(double R0, double Th);
	/**
	 * Рассчитать радиус миделева сечения в заданной по оси точке
	 * @param x - координата точки вдоль оси симметрии.
	 */
	double mr(double x) const;
	/**
	 * Рассчитать угол между образующей тела и набегающим потоком в заданной точке, в радианах.
	 * @param x - координата точки вдоль оси симметрии.
	 */
	double theta(double x) const;
	/**
	 * Координату x по известной координате y точки.
	 * @param y - координата точки вдоль оси симметрии.
	 */
	double x(double y) const;
	/**
	 * Рассчитать расстояние от критической точки вдоль образующей до точки с координатой x вдоль оси симметрии.
	 * @param x - координата точки вдоль оси симметрии.
	 */
	double xgen(double x) const;
	/**
	 * Рассчитать площадь миделева сечения.
	 * @param x - координата точки вдоль оси симметрии.
	 */
	double Sbottom(double x) const;
	double R();
};

#endif /* BLUNTEDCONE_H */