#include <avdsolver.h>
#include <boundary.h>

AVDSolver::AVDSolver(thm_t* thm, trm_t* trm, gasdynamics_t* gd, CBluntedCone* BCone)
{
	assert(thm != 0);
	assert(trm != 0);
	assert(gd != 0);

	this->thm = thm;
	this->trm = trm;
	this->gd = gd;
	this->BCone = BCone;
	this->thsolver = new CTHSolver(thm);
	this->thsolver->setPrefs(STD_TIMESTEP_MIN, STD_TIMESTEP_MAX, 0.);
	this->TIMESTEP = STD_TIMESTEP_MIN;
}


void BCA(double HB, double* TB, double* PHB, double* ROHB, double* F) {
	double G0 = 9.80665;
	double R1 = 8314.32;
	double R = 6356766.0;
	double HM1 = 28.96442;
	double A[13][4] = {	{-2000.,	301.15,	-0.0065,	127774.},
				{-0.,		288.15,	-0.0065,	101325.},
				{11000.,	216.65,	0.0,		22632.},
				{20000.,	216.65,	0.001,		5474.87},
				{32000.,	228.65,	0.0028,		868.014},
				{47000.,	270.65,	0.0,		110.906},
				{51000.,	270.65,	-0.0028,	66.9384},
				{71000.,	214.65,	-0.002,		3.95639},
				{80000.,	196.65,	-0.002,		0.886272},
				{85000.,	186.65,	0.,		0.363403},
				{94000.,	186.65,	0.003,		0.0699754},
				{102450.,	212.,	0.011,		0.0164122},
				{117777.,	380.6,	0.011,		0.00266618}	};

	double H1=R*HB/(R+HB);
	double HM;
	if (H1 <= 85000.0)
		HM=28.96442;
	else if (H1 <= 94000.0)
		HM=28.96442-0.00942*(H1-85000.0)/9000.0;
	else if (H1 <= 102450.0)
		HM=28.955-1.109*(H1-94000.0)/8450.0;
	else
		HM=27.846-2.396*(H1-102450.0)/15327.0;
	double RT=R1/HM;
	int i;
	for (i = 1; (A[i][0] <= H1) && (i<12); i++);
	double T=A[i-1][1]+(A[i][1]-A[i-1][1])*(H1-A[i-1][0])/(A[i][0]-A[i-1][0]);
      *TB=T*HM/HM1;
      *F=20.046796*pow(*TB,0.5);
      if (A[i-1][2] != 0.0) {
		if (A[i-1][2] > 0)
			*PHB=A[i-1][3]/pow(1+A[i-1][2]*(H1-A[i-1][0])/A[i-1][1], G0/(A[i-1][2]*RT));
		else
			*PHB=A[i-1][3]*pow(1+A[i-1][2]*(H1-A[i-1][0])/A[i-1][1], -G0/(A[i-1][2]*RT));
      } else
         *PHB=A[i-1][3]*exp(-G0*(H1-A[i-1][0])/(RT*A[i-1][1]));
      *PHB=*PHB/G0;
      *ROHB=*PHB/(RT*(*TB));
      return;
}

/* 
 * P = P/P0
 * O - ���� ������������, ����
*/
static avd_t old_avd(double TB, double ROH, double VT, double M, double H, double P, double TW, double X, double O, double Rad, double XEF, int LT, int ICC, double GK)
{
//	double TB, PH, ROH, D;
//	BCA(H, &TB,&PH, &ROH, &D);
//	double M = VT/D;
	double ALC;
	double OI, WI, EI, XAP1;
	OI=0.24*TB+VT*VT/8370.0;
	double AM1 = M*M;
	double AM = M;
	double XAP = 1.4;
	double C3=pow(AM1, XAP/(XAP-1.))/pow(((2.*XAP/(XAP-1.))*AM1-1.), (1./(XAP-1.)));
	double CX1=pow(((XAP+1)/2.), ((XAP+1.)/(XAP-1.)));
	double CX2=pow((2./(XAP-1.)), (1./(XAP-1.)));
	double CX3=pow(6, 0.25);
	double CX4=pow(1.2, 6)*pow(5, 2.5);
	double CX5=1.0/729.0/0.85;
	double P11=CX1*CX2*C3;
	double PP=ROH*28.3*TB*9.80665*0.0001;
	double P0=P11*PP;
	double P1=P*P0;
	double ALC0, SIB;
	double DK, EK, F1, F2, F3, V0;
	if (M <= 2.0) {
		ALC0=0.0007*pow(ROH*VT, 0.8)/pow(X, 0.2)*pow(((TB+TW)/(2.0*TB)), 2.137)*pow(0.24*TB, 0.137);
		WI=0.24*TW;
		EI=OI;
	} else {
		C3=1./pow((1.+(XAP-1.)/2.*AM1), 0.6);
		double C1=pow(((1.+(XAP-1.)/2.*AM1)/((XAP+1.)/2.*AM1)), 0.4);
		double C2=1./pow((XAP*AM1), 0.2);
		F3=pow(P11, 0.8)*C1*C2*C3;
		F1=1.+0.137*sin(1.57*(1.-(TW-273.)/1000.));
		if (TW <= 1000.0)
			WI=0.245*TW;
		else if (TW <= 2500.0)
			WI=245.+0.310*(TW-1000.);
		else
			WI=TW*TW/8800./pow((1.+0.05*log10(P0)), 2);
		double XM=3.*pow((24./O+4.), (3.-9./AM));
		double XAP0=1.23;
		XAP1=XAP0;
		V0=1.15*(cos(O*0.0174))-0.055*pow((O/10.), 0.97)/AM;
		double X1=(X/Rad-1.57*(1.0-O/90.))*cos(O*0.0174);
		double OM;
		if (V0 < 1.0)
			V0=1.0;
		do {
			double B0;
			OM=1./(pow(P, ((XAP1-1.)/XAP1)))-1.;
			C1=1.-(XAP-1.)/2.*AM1*(1.-pow(V0, 2));
			C2=(XAP-1.)/2.*AM1*pow(V0, 2);
			double OMM=C2/C1/OM;
			if (X1 < XM)
				B0=1.0;
			else
				B0=(X1-XM)/2./XM*(OMM-1.)+1.;
	
			if (abs(B0) > abs(OMM))
				B0=OMM;
			double B00=(1.-WI/OI)/(1.+WI/OI);
			if (B0*OM > B00)
				SIB=OI/4.*((B0*OM+1.)/B0/OM)*pow((1.-WI/OI), 2)+WI;
			else
				SIB=OI/(1.+OM*B0);
			if (SIB > 500.0)
				XAP0=1.23-96./SIB+6E+4/pow(SIB, 2);
			else
				XAP0=1.4-0.12*SIB/500.;
			XAP1=XAP0;
		} while (abs(XAP0-XAP1) >= 0.01);

		F2=pow(((XAP+1.)/(XAP-1.)), 0.4)*pow((1.-pow(P, ((XAP-1.)/XAP))), 0.4)*pow(P, 0.8)*pow((OI/SIB), 0.6);
		DK=0.78*pow(SIB, 0.05);
		if (VT <= 1000.0)
			DK=1.0;
		double EKM=1.11+0.0055*O;
		EK=(X1-XM)/2./XM*(EKM-1.)+1.;
		if (EK > EKM)
			EK=EKM;
		if (X1 <= XM)
			EK=1.0;
		EI=OI*(1.+0.89*OM)/(1.+OM);
		ALC0=0.00243*pow(VT, 1.2)*pow(ROH, 0.8)*F1*F2*F3*DK*EK/XEF;
		if (LT != 0) {
			double F1L=1.0;
			if (SIB < 1800.0)
				F1L=1.45/pow(WI, 0.06);
			double F2L=CX3*pow((1.-pow(P, 0.288)), 0.25)*pow(P, 0.5);
			double F3L=pow((CX4*pow(AM, 7.)/pow((7.*AM1-1.), 2.5)), 0.5)*pow(((1.+0.2*AM1)/1.2/AM1), 0.25)/AM;
			if (SIB >= 1800.0)
				DK=6.4/pow(SIB, 0.25);
			else
				DK=3.42/pow(SIB, 0.19);
			ALC0=1.08E-5*pow(VT, 1.5)*sqrt(ROH)*F1L*F2L*F3L/XEF*EK*DK;
		}
	}
	ALC=ALC0;
	double QW, ALC1 = ALC;
	if ((ICC == 0) || (TW < 1000.0)) {
		QW=ALC*(EI-WI);
	} else {
		ALC1=ALC*(pow((0.5+sqrt(0.25+GK*CX5)), 0.333)-pow((sqrt(0.25+GK*CX5)-0.5), 0.333));
		QW=1.3*ALC1*(EI-WI);
	}
	avd_t avd;
	avd.I0 = OI*4186.8;
	avd.IW = WI*4186.8;
	avd.IE = EI*4186.8;
	avd.ISTAR = SIB*4186.8;
	avd.ALC = ALC;
	avd.ALC1 = ALC1;
	avd.QCONV = QW*4186.8;
	avd.XAP1 = XAP1;
	avd.P0 = P0;
	avd.P1 = P1;
	avd.V0 = V0;
	avd.F1 = F1;
	avd.F2 = F2;
	avd.F3 = F3;
	avd.KDIS = DK;
	avd.KENTH = EK;
	return avd;
}

avdsolver_t AVDSolver::Solve(double time)
{
	assert(time > thm->CURRENT_TIME);
	double G1 = 0.;
	avdsolver_t info;
	info.G = 0.;
	func_points_t points;
	
	points.count = 20;
	points.x = new double [20];
	points.y = new double [20];
	FILE *fQ = fopen("q.txt", "rt");
	for (int i=0; i< points.count; i++)
		fscanf(fQ, "%lf\t%lf\n", &(points.x[i]), &(points.y[i]));
	fclose(fQ);
	for (; thm->CURRENT_TIME < time; )
	{
		double AT = thm->m[thm->fcnum]->at(); /* Get ablation type of surface */
		CBoundary *bc;
		
		TIMESTEP = min(TIMESTEP, time - thm->CURRENT_TIME);
		/* ��������� ������������ �� ��������� ��������. */
		info.H = trm->H.val(thm->CURRENT_TIME);
		double TB, PH, ROH, D;
		BCA(info.H, &TB, &PH, &ROH, &D);
		info.V =  trm->V.val(thm->CURRENT_TIME);
		info.al =  trm->AL.val(thm->CURRENT_TIME);
		info.phi =  fabs(fmod(trm->PHI.val(thm->CURRENT_TIME), 360.));	
		if (info.phi > 180.)
			info.phi = 360. - info.phi;
		info.mach = info.V/D;
		bool isTurbulent;
		if (info.H < gd->HT) {
			info.XEF = gd->XET.val(info.mach, info.al, info.phi);
			isTurbulent = true;
		} else {
			info.XEF = gd->XEL.val(info.mach, info.al, info.phi);
			isTurbulent = false;
		}
		info.PP0 = gd->PP0.val(info.mach, info.al, info.phi);
		
		info.avd = old_avd(TB, ROH, info.V, info.mach, info.H, info.PP0, thm->TWL, gd->X, BCone->theta(gd->X), BCone->R(), info.XEF, !isTurbulent, AT, info.G);
		info.avd.QCONV = in_LinearFunc(&points, thm->CURRENT_TIME, 0);
		info.avd.ALC = info.avd.QCONV/(info.avd.IE-info.avd.IW);
		info.avd.ALC1 = info.avd.ALC;
		info.srt.QLrad = 0.; info.srt.QRrad = 0.; info.srt.QLconv = 0.; info.srt.QRconv = 0.;
		double Tw = thm->TWL;
		double A = thm->m[thm->fcnum]->a(Tw);
		if (A == 0) {
			printf("Ablation A-coefficient can't be zero. Check material properties at source file!\n");
			exit(-1);
		}
		double B = thm->m[thm->fcnum]->b(Tw);
		CBoundary* tmpbc;
		tmpbc = thm->LBC;
		TIMESTEP = min(TIMESTEP, STD_TIMESTEP_MAX);
		TIMESTEP = max(TIMESTEP, STD_TIMESTEP_MIN);

		if ((AT == 0) && (Tw >= thm->m[thm->fcnum]->Td(Tw)-20.0) && (info.avd.IE > info.avd.IW))
			bc = new CFOBoundary(thm->m[thm->fcnum]->Td(0));
		else
			bc = new CSOBoundary(info.avd.QCONV/(info.avd.IE-info.avd.IW), info.avd.I0, info.avd.IE, info.avd.IW, info.avd.P1*101325., thm->m[thm->fcnum]->eps(thm->TWL));
		thm->setLBC(bc);
		assert(TIMESTEP > 0.);
		solve_result_t srt = thsolver->Solve(thm->CURRENT_TIME+TIMESTEP);
		Tw = thm->TWL;
		info.srt = srt;
		info.srt.QLrad += srt.QLrad; info.srt.QRrad += srt.QRrad; info.srt.QLconv += srt.QLconv; info.srt.QRconv += srt.QRconv;
		/* ������ ����� */
		if (AT == 1)
		{
			if (Tw > 1000.)
			{
				info.G = A+B*info.avd.I0/4186.8;
				double r = thm->m[thm->fcnum]->r(Tw);
				
				if (Tw <= 3000.)
					G1 = 0.043*(A/0.19)*sqrt(16.-pow(6.-(Tw)/500., 2.))*(1.+1.4e7/(pow(info.avd.P1, 0.67)*exp(6.14e4/(Tw))));
				else if (Tw > 3000.) {
					G1 = 0.172*(A/0.19)*(1.+1.4e7/(pow(info.avd.P1, 0.67)*exp(6.14e4/(Tw))));
				}
			
				assert(G1 >= 0.);
				if (G1 < info.G)
					info.G = G1;
				else
					thm->TWL = 6.14E+04/log(2.4e6/(fabs(info.G-0.172*(A/0.19))*pow(info.avd.P1, 0.67)));
				if (Tw > 4500.)
					printf("Tw=%lf T=%lf G=%lf G1=%lf Ps=%lf w=%lf TIMESTEP=%E fcnum=%d TWL=%lf\n", Tw, thm->T[thm->fcnum], info.G, G1, info.avd.P1/101325., thm->width[thm->fcnum], TIMESTEP, thm->fcnum, thm->TWL);
				double Vdx = info.G*(info.avd.ALC1)/r;
				thm->crop(0, Vdx*TIMESTEP);
			}
		} else if ((AT == 0) && (Tw >= thm->m[thm->fcnum]->Td(Tw)-20.0) && (info.avd.IE > info.avd.IW))
		{
			double Vdx;
			double r = thm->m[thm->fcnum]->r(Tw);
			
			info.G = B*4186.8+A*(info.avd.IE-info.avd.IW);
			Vdx = (info.avd.QCONV)/r/info.G;
//			printf("A=%lf B=%lg IE=%lf Vdx=%lf TIMESTEP=%E\n", A, B, info.avd.IE, Vdx, TIMESTEP);
			thm->crop(0, Vdx*TIMESTEP);
		}
		assert(TIMESTEP > 1.0E-20);
		thm->setLBC(tmpbc);
		if (srt.CURRENT_DT_MAX > 10.)
			TIMESTEP /= 2.;
		else
			TIMESTEP *= 1.2;
	}
	
	info.time = thm->CURRENT_TIME;
	return info;
}
AVDSolver::~AVDSolver()
{
}
void AVDSolver::print()
{
	thsolver->print();
}
