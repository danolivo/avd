#include <model.h>
#include <cstring>

thm_t::thm_t()
{
	lcnum = -1;
	fcnum = -1;
	CURRENT_TIME = 0.;
	INIT_TIMESTEP = 1.0;
	LDEL = 0.;
	TWL = NULL_VALUE;
	TWR = NULL_VALUE;
	LBC = 0;
	RBC = 0;
}
int thm_t::add(int Cells, double Width, CMaterial* Material, double Temp)
{
	if (lcnum+Cells >= CELLS_MAX_NUM)
		return NULL_VALUE;
	
	if (fcnum == NULL_VALUE) {
		fcnum = 0;
		lcnum = -1;
	}
	/* Init layer properties */
	for (int i = 1; i<=Cells; i++) {
		lcnum++;
		width[lcnum] = Width/Cells;
		m[lcnum] = Material;
		T[lcnum] = Temp;
	}
	TWL = T[fcnum];
	TWR = T[lcnum];
	PrimaryLeftCellSize = width[fcnum];
	PrimaryRightCellSize = width[lcnum];
	return SUCCESS;
}
void thm_t::setLBC(CBoundary* BC)
{
	this->LBC = BC;
}
void thm_t::setRBC(CBoundary* BC)
{
	this->RBC = BC;
}
double thm_t::length() {
	int i;
	double len;

	assert(fcnum >= 0);
	assert((lcnum >= 0) && (lcnum < CELLS_MAX_NUM));
	assert((lcnum - fcnum) >= 0);

	for (i=fcnum, len=0.; i<=lcnum; i++) {
		assert(width[i] > 0.);
		len += width[i];
	}

	return len;
}
int thm_t::crop(int side, double dx) {
	double del;
	int i;
	double x0 = 0.;
	double y0 = TWL;
	double x1 = width[fcnum]/2.;
	double y1 = T[fcnum];
	double x2 = width[fcnum]+width[fcnum+1]/2.;
	double y2 = T[fcnum+1];
	assert((dx >= 0.) || !printf("current value of dx is: %f\n", dx));
	assert(side == 0);
	if (dx == 0.)
		return 0;
	assert(side == 0);
	// Отрезать слева
	for (i=fcnum, del=dx; del >= width[i]; i++) {
		del -= width[i];
		assert(i < lcnum);
	}
	/* Интерполировать (экстраполировать) температуру по квадратичному закону. */
	double T2;

		
	if (i > fcnum) {
		TWL = in_linear(x1, y1, x2, y2, dx);
		T2 = in_linear(x1, y1, x2, y2, x2+del/2.);
	} else {
		double a = y0;
		double b = log(y1/a)/x1;
		TWL = a*exp(b*dx);	
		assert(TWL > 0.);
		T2 = in_linear(x1, y1, x2, y2, x1+del/2.);
	}
	width[i] -= del;
	// Коррекция погрешностей
	if (width[i]/PrimaryLeftCellSize < 0.01) {
		i++;
		width[i] += width[i-1];
		assert(T2 > 0.);
		T2 = in_linear(x1, y1, x2, y2, x2-(width[i-1])/2.);
		assert(T2 > 0.);
	}
	// Удалить воздушный зазор.
	while (strcmp(m[i]->name(), "Air") == 0) {
		LDEL += width[i];
		i++;
	}
	if (i > fcnum)
		PrimaryLeftCellSize = width[i];
	assert(TWL > 0.);
	assert(T2 > 0.);
	fcnum = i;
	LDEL += dx;
	T[i] = T2; 
	return 0;
}

void thm_t::setTemperature(func_points_t * Tf)
{
	double x = 0.;
	for (int i=fcnum; i<=lcnum; i++) {
		x += width[fcnum]/2.;
		T[i] = in_LinearFunc(Tf, x);
		x += width[i]/2.;
	}
}
void thm_t::print()
{
	printf("T= ");
	for (int i=fcnum; i<=lcnum; i++) {
		printf("%lf ", T[i]);
	}
	printf("\n");
	printf("w= ");
	for (int i=fcnum; i<=lcnum; i++) {
		printf("%5.3E ", width[i]);
	}
	printf("\n");
}

