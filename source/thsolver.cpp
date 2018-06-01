#include <tdma.h>
#include <thsolver.h>
#include <boundary.h>
#include <cstring>

#define NEED_SMALLER_STEP	(1)
#define TIMESTEP_MIN_CALCS	(2)

CTHSolver::CTHSolver(thm_t* thm)
{
	assert(thm != 0);
	this->thm = thm;
	this->TIMESTEP_MIN = STD_TIMESTEP_MIN;
	this->TIMESTEP_MAX = STD_TIMESTEP_MAX;
	this->T[0] = thm->T[thm->fcnum];
	for (int i=1; i< CELLS_MAX_NUM; i++) {
		T[i] = -1.;
		l[i] = -1.;
		c[i] = -1.;
		r[i] = -1.;
		qv[i] = 0.;
	}
}

CTHSolver::~CTHSolver()
{
}
	
solve_result_t CTHSolver::Solve(double time)
{
	sres.QLconv = 0.;
	sres.QLrad = 0.;
	sres.QRconv = 0.;
	sres.QRrad = 0.;
	sres.dHeatQty = 0.;
	double dHeatQty = 0.;
	double residual = 0.;
	int counter = 0;
	double dt = time - thm->CURRENT_TIME;
	Prepare();
	double TIMESTEP = time - thm->CURRENT_TIME;
	bool isPrintWTmin = false;
	assert(time - thm->CURRENT_TIME > 0.);
	while (thm->CURRENT_TIME < time) {
		double HeatQty_Pre = Pre(thm->CURRENT_TIME);
		double twlk = Twl();
		double twrk = Twr();
		int res = DoIteration(TIMESTEP);
		if (res == NEED_SMALLER_STEP) {
			Prepare();
			if (TIMESTEP > TIMESTEP_MIN)
				TIMESTEP = TIMESTEP/2.;
			continue;
		} else if ((res == SUCCESS) || (res == TIMESTEP_MIN_CALCS)) {
			thm->CURRENT_TIME += TIMESTEP;
			counter++;
			dHeatQty = Post(thm->CURRENT_TIME)-HeatQty_Pre;
			/* Calculate heat quantity at boundaries */
			double ql;
			if (thm->LBC->type() != 1) { // Second-order left BC
				BC_t bc = thm->LBC->GetBC(thm->CURRENT_TIME-TIMESTEP);
				double qlconv, qlrad;
				qlconv = bc.acp*(bc.Ie-bc.Iw);
				sres.QLconv+= qlconv*TIMESTEP;
				qlrad = (-bc.eps*5.67E-08*(pow(twlk, 4.)+4*pow(twlk, 3.)*(Twl()-twlk)));
				sres.QLrad += qlrad*TIMESTEP;
				ql = qlconv+qlrad;
			} else
				ql = 2*(twlk-T[1])*l[1]/w[1];
			double qr;
			if (thm->RBC->type() != 1) { // Second-order left BC
				BC_t bc = thm->RBC->GetBC(thm->CURRENT_TIME-TIMESTEP);
				double qrconv, qrrad;
				qrconv = bc.acp*(bc.Ie-bc.Iw);
				sres.QRconv+= qrconv*TIMESTEP;
				qrrad = (-bc.eps*5.67E-08*(pow(twrk, 4.)+4*TIMESTEP*pow(twrk, 3.)*(Twr()-twrk)/TIMESTEP));
				sres.QRrad += qrrad*TIMESTEP;
				qr = qrconv+qrrad;
			} else
				qr = 2*(twrk-T[SIZE-2])*l[SIZE-2]/w[SIZE-2];
			sres.dHeatQty += (ql+qr)*TIMESTEP-dHeatQty;
		} else
			assert(0 || !printf("[EE]: Unknown result value: %d\n", res));
		TIMESTEP = min(TIMESTEP*1.1, TIMESTEP_MAX);
		TIMESTEP = min(TIMESTEP, time-thm->CURRENT_TIME);
	}
	sres.dHeatQty/=(sres.QLconv+sres.QLrad);
	sres.LAST_TIMESTEP = max(dt/counter, STD_TIMESTEP_MIN);
	TIMESTEP_MIN = sres.LAST_TIMESTEP;
	return sres;
}

void CTHSolver::Prepare()
{
	assert(thm->lcnum < CELLS_MAX_NUM+2);
	assert(thm->fcnum >= 0);
	SIZE = thm->lcnum-thm->fcnum+3;
	w[0]= thm->PrimaryLeftCellSize/100.;
	w[SIZE-1]= thm->PrimaryRightCellSize/100.;
	T[0]= thm->TWL;
	T[SIZE-1]= thm->TWR;
	/* Make a internal copy of thermal model */
	memcpy(&(T[1]), &(thm->T[thm->fcnum]), (SIZE-2)*sizeof(double));
	memcpy(&(w[1]), &(thm->width[thm->fcnum]), (SIZE-2)*sizeof(double));
}
double CTHSolver::Pre(double time)
{
//	printf("PRE:\n");
	/* Update thermal properties */
	for (int j=1; j<SIZE-1; j++) {
		int i = j+thm->fcnum-1;
		l[j] = thm->m[i]->l(T[j]);
		c[j] = thm->m[i]->c(T[j]);
		r[j] = thm->m[i]->r(T[j]);
	}
	/* Update boundaries */
	setBoundaries(thm->LBC, thm->RBC, time);
	/** Calculate heat quantity */
	return currentHeatQty();
}

double CTHSolver::Post(double time)
{
	memcpy(&(thm->T[thm->fcnum]), &(T[1]), (SIZE-2)*sizeof(double));
	thm->TWL = Twl();
	if (thm->TWL <= 0.) {
		printf("T0=%lf T1=%lf T2=%lf w0=%E w1=%lf w2=%lf qv0=%lf\n", T[0], T[1], T[2], w[0], w[1], w[2], qv[0]);
		print();
	}
	assert(thm->TWL > 0.);
	thm->TWR = Twr();
	return currentHeatQty();
}
int CTHSolver::DoIteration(double TIMESTEP)
{
	double tempT0 = T[0];
	CalculateTDMA(T, w, qv, l, r, c, eps, TIMESTEP, SIZE);
	if (TIMESTEP <= TIMESTEP_MIN)
		return TIMESTEP_MIN_CALCS;
	if (T[0] <= 0.) {
			print();
			thm->print();
			printf("[EE]: cell %d has T=%lf! Previous T=%lf TIMESTEP=%lf\n", 0, T[0], thm->T[1], TIMESTEP);
			fflush(NULL);
	}
	if ((fabs(T[0]-tempT0) > DT_MAX) && (DT_MAX > 0.))
		return NEED_SMALLER_STEP;
	sres.CURRENT_DT_MAX = fabs(T[0] - tempT0);
	
	for (int i=1; i<SIZE-2; i++) {
		int j=thm->fcnum+i-1;
		if (T[i] <= 0.) {
			printf("i=%d T[i]=%lf\n", i, T[i]);
			print();
			thm->print();
		}
		assert(T[i] >= 0.);
		double DT = fabs(T[i]-thm->T[j]);
		if (DT > sres.CURRENT_DT_MAX)
			sres.CURRENT_DT_MAX = DT;
		if ((DT > DT_MAX) && (DT_MAX > 0.))
			return NEED_SMALLER_STEP;
		if (T[i] < 0.) {
			print();
			thm->print();
			printf("[EE]: cell %d has T=%lf! Previous T=%lf TIMESTEP=%lf\n", i, T[i], thm->T[j], TIMESTEP);
			fflush(NULL);
		}
		assert((T[i] >= 0.));	
	}
	return SUCCESS;
}
void CTHSolver::setBoundaries(CBoundary *lbc, CBoundary *rbc, double time)
{
	setLeftBoundary(lbc, time);
	setRightBoundary(rbc, time);
}
void CTHSolver::setLeftBoundary(CBoundary *cbc, double time)
{
	BC_t bc = cbc->GetBC(time);
	if (cbc->type() == 1) { /* First kind type boundary */
		l[0] = 1.0E+15*l[1];
		c[0] = 1.0E+03*c[1];
		r[0] = 1.0E+03*r[1];
		assert(bc.Tw >= 0.);
		T[0] = bc.Tw;
		qv[0] = 0.;
		eps[0] = 0.;
	} else {
		assert(bc.Ie >= 0.); assert(bc.eps >= 0.);
		double q = bc.acp*(bc.Ie-bc.Iw);
		T[0] = thm->TWL;
		c[0] = c[1];
		r[0] = r[1];
		l[0] = l[1];
		qv[0] = q/w[0];
		eps[0] = bc.eps;
//		assert(bc.eps == 0);
	}
}
void CTHSolver::setRightBoundary(CBoundary *cbc, double time)
{
	BC_t bc = cbc->GetBC(time);
	if (cbc->type() == 1) { /* First kind type boundary */
		l[SIZE-1] = 1.0E+15*l[SIZE-2];
		c[SIZE-1] = 1.0E+03*c[SIZE-2];
		r[SIZE-1] = 1.0E+03*r[SIZE-2];
		assert(bc.Tw >= 0.);
		T[SIZE-1] = bc.Tw;
		qv[SIZE-1] = 0.;
		eps[1] = 0.;
	} else {
		l[SIZE-1] = l[SIZE-2];
		c[SIZE-1] =c[SIZE-2];
		r[SIZE-1] = r[SIZE-2];
		assert(bc.acp >= 0.); assert(bc.Ie >= 0.); assert(bc.eps >= 0.);
		double q = bc.acp*(bc.Ie-bc.Iw);
		qv[SIZE-1] = q/w[SIZE-1];
		T[SIZE-1] = thm->TWR;
		eps[1] = bc.eps;
	}
}

thm_t* CTHSolver::getTHM()
{
	return thm;
}

void CTHSolver::setPrefs(double TIMESTEP_MIN, double TIMESTEP_MAX, double DT_MAX)
{
	this->TIMESTEP_MIN = TIMESTEP_MIN;
	this->TIMESTEP_MAX = TIMESTEP_MAX;
	this->DT_MAX = DT_MAX;
}
void CTHSolver::printPreferences()
{
	printf("SOLVER PREFERENCES:\n");
	printf("CURRENT_TIME=%lf\n", thm->CURRENT_TIME);
	printf("TIMESTEP_MIN=%E\n", TIMESTEP_MIN);
	printf("TIMESTEP_MAX=%E\n", TIMESTEP_MAX);
}
double CTHSolver::Twl()
{
	return T[0];
}
double CTHSolver::Twr()
{
	return T[SIZE-1];
}
double CTHSolver::currentHeatQty()
{
	double sum =  w[0]*thm->m[thm->fcnum]->heatQuantity(T[0]);
	for (int i=1; i<SIZE-1; i++)
		sum += w[i]*thm->m[thm->fcnum+i-1]->heatQuantity(T[i]);
	sum +=  w[SIZE-1]*thm->m[thm->lcnum]->heatQuantity(T[SIZE-1]);
	return sum;
}
void CTHSolver::print(FILE* out)
{
	fprintf(out, "T= ");
	for (int i=0; i<SIZE-1; i++)
		fprintf(out, "%10.4lf ", T[i]);
	fprintf(out, "\n");
	fprintf(out, "w= ");
	for (int i=0; i<SIZE-1; i++)
		fprintf(out, "%5.3E ", w[i]);
	fprintf(out, "\n");
	fprintf(out, "l= ");
	for (int i=0; i<SIZE-1; i++)
		fprintf(out, "%5.3E ", l[i]);
	fprintf(out, "\n");
	fprintf(out, "c= ");
	for (int i=0; i<SIZE-1; i++)
		fprintf(out, "%5.3E ", c[i]);
	fprintf(out, "\n");
	fprintf(out, "r= ");
	for (int i=0; i<SIZE-1; i++)
		fprintf(out, "%5.3E ", r[i]);
	fprintf(out, "\n");
}
