#include "stdafx.h"


extern double agree_param;
extern double Eps;


SignChecker::SignChecker()
{
	this->DTWwindow = 1;

	this->SimpleChecks = 4;
	this->DTWChecks = 3;

	this->SimpleTestResult = 0.0;
	this->DTWTestResult = 0.0;

	this->mySimpleCheckFunctions = std::unique_ptr<SimpleCheckFunctions[]>(new SimpleCheckFunctions[SimpleChecks]);
	
	mySimpleCheckFunctions[0] = &RatioXY;
	mySimpleCheckFunctions[1] = &AvgXY;
	mySimpleCheckFunctions[2] = &SinXY;
	mySimpleCheckFunctions[3] = &RatioTouches;

	this->SimpleM = std::unique_ptr<Matrix<double>[]>(new Matrix<double>[this->SimpleChecks]);
	this->SimpleMch = std::unique_ptr<Matrix<double>[]>(new Matrix<double>[this->SimpleChecks]);

	this->myDTWCheckFunctions = std::unique_ptr<DTWCheckFunctions[]>(new DTWCheckFunctions[DTWChecks]);
	
	myDTWCheckFunctions[0] = &CityBlockMetric;
	//myDTWCheckFunctions[1] = &EuclideMetric;
	myDTWCheckFunctions[2] = &VelocityMetric;
	myDTWCheckFunctions[1] = &SinMetric;

	this->DTWch = std::unique_ptr<Matrix<double>[]>(new Matrix<double>[this->DTWChecks]);
	this->DTW = std::unique_ptr<Matrix<double>[]>(new Matrix<double>[this->DTWChecks]);
	
	this->CGch = std::unique_ptr<Matrix<double>[]>(new Matrix<double>[this->DTWChecks]);
	this->CG = std::unique_ptr<Matrix<double>[]>(new Matrix<double>[this->DTWChecks]);
	
	this->CVch = std::unique_ptr<Matrix<double>[]>(new Matrix<double>[this->DTWChecks]);
	
	this->Ncg = std::unique_ptr<Matrix<int>[]>(new Matrix<int>[this->DTWChecks]);
}

void SignChecker::GenerateMatixeByIPen(int icheck, CTab *ipens, int npens)
{
	double argchi, argchj, argval;

	for (int i = 0; i < npens; i++)
	{
		argchi = (*mySimpleCheckFunctions[icheck])(ipens, i);
		for (int j = 0; j < i; j++)
		{
			argchj = (*mySimpleCheckFunctions[icheck])(ipens, j);
			argval = abs(argchi - argchj);
			if (npens < 3) argval /= std::max(argchi, argchj);
			this->SimpleMch[icheck].SetElem(i,j,argval);
		}
	}
}

bool SignChecker::SimpleCheckForRandomForge(CTab *ipens, int npens)
{
	for (int i = 0; i < this->SimpleChecks; i++)
	{
		this->SimpleMch[i].Init(npens, npens);
		this->GenerateMatixeByIPen(i, ipens, npens);
		this->SimpleM[i].Init(this->SimpleMch[i], 1, 1);
		this->SimpleTestResult += 1 - CalcDifference(this->SimpleMch[i], this->SimpleM[i], npens);
	}
	this->SimpleTestResult /= double(this->SimpleChecks);
	return true;
}

bool SignChecker::DTWCheckForSimpleForge(DPoints *dpens, SPoints *spens, int npens)
{

	for (int i = 0; i < this->DTWChecks; i++)
	{
		this->DTWch[i].Init(npens, npens);
		this->CGch[i].Init(npens, npens);
		this->CVch[i].Init(npens, npens);
		this->Ncg[i].Init(npens, npens);
		this->GenerateMatrixByDPen(i, dpens, spens, npens);
		this->DTW[i].Init(this->DTWch[i], 1, 1);
		this->CG[i].Init(this->CGch[i], 1, 1);

		this->DTWTestResult += 1 - this->DTWCalcDifference(i, npens);
	}
	this->DTWTestResult /= double(this->DTWChecks);
	return true;
}

double SignChecker::CalcDifference( Matrix<double>& ratioM_withCheckedSign, Matrix<double>& ratioM_withoutCheckedSign, int npens)
{
	int imax = 0;
	int jmax = 0;
	double elmaxch, elavgch, elavg, ratval;
	
	elmaxch = ratioM_withCheckedSign.GetIndMaxElem(imax, jmax);
	
	elavgch = ratioM_withCheckedSign.AvgByI(npens - 1);
	elavg = ratioM_withoutCheckedSign.AvgAll();

	if (imax == npens - 1)
		elavgch = elmaxch;

	if (npens < 3)
	{
		if (elavgch > 1)
			ratval = 1.0 - (1.0 / elavgch);
		else
			ratval = elavgch;
	}
	else
	{
		ratval = abs(elavgch - elavg) / std::max(elavgch, std::max(elavg, Eps));
	}

	return ratval;
}

double SignChecker::DTWCalcDifference(int icheck, int npens)
{
	double ratval;

	int jCVcheckmax = this->CVch[icheck].GetJMaxElInI(npens - 1);

	double DTWavg = this->DTW[icheck].AvgByI(jCVcheckmax);
	double sigma = sqrt(this->DTW[icheck].SumSqByI(jCVcheckmax, DTWavg));
	double sigmacheck = abs(this->DTWch[icheck](npens - 1, jCVcheckmax) - DTWavg);
	
	if ((sigma == 0) && (sigmacheck == 0))
		ratval = 0;
	else
		ratval = abs(sigma - sigmacheck) / std::max(sigmacheck, sigma);

	double CGavg = this->CG[icheck].AvgByI(jCVcheckmax);
	sigma = sqrt(this->CG[icheck].SumSqByI(jCVcheckmax, CGavg));
	sigmacheck = abs(this->CGch[icheck](npens - 1, jCVcheckmax) - CGavg);

	if ((sigma == 0) && (sigmacheck == 0))
		ratval += 0;
	else
		ratval += abs(sigma - sigmacheck) / std::max(sigmacheck, sigma);

	return ratval/2.0;
}

void SignChecker::GenerateMatrixByDPen(int icheck, DPoints* dpens, SPoints *spens, int npens)
{
	for (int i = 0; i < npens; i++)
		for (int j = 0; j < i; j++)
			this->DTW_Go(icheck, dpens, spens, i, j);
}

void SignChecker::DTW_Go(int icheck, DPoints* dpens, SPoints *spens, int i, int j)
{
	int Ni, Nj, Nt;
	Ni = spens[i].GetN();
	Nj = spens[j].GetN();
	Matrix<double> D(Ni, Nj);
	SPoints TPeni, TPenj;

	if (this->DTW_InitMatrix(icheck, D, dpens[i], spens[i], dpens[j], spens[j])) 
	{
		this->DTW_TraceBack(D, TPeni, TPenj);
		this->DTWch[icheck].SetElem(i, j, D(Ni - 1, Nj - 1));
		this->DTWch[icheck].SetElem(j, i, D(Ni - 1, Nj - 1));
		
		Nt = TPeni.GetN();
		this->Ncg[icheck].SetElem(i, j, Nt);
		this->Ncg[icheck].SetElem(j, i, Nt);

		if (Nt == 0) Nt = 1;
		this->CGch[icheck].SetElem(i, j, this->DTW_CaclGlobalDeformation(D, TPeni, TPenj)/double(Nt));
		this->CGch[icheck].SetElem(j, i, this->DTW_CaclGlobalDeformation(D, TPeni, TPenj) / double(Nt));

		this->CVch[icheck].SetElem(i, j, this->DTW_CalcDetermination(TPeni, dpens[i], TPenj, dpens[j]));
		this->CVch[icheck].SetElem(j, i, this->DTW_CalcDetermination(TPeni, dpens[i], TPenj, dpens[j]));
	}
}

bool SignChecker::DTW_InitMatrix(int icheck, Matrix<double>& D, DPoints& dpeni, SPoints& speni, DPoints& dpenj, SPoints& spenj)
{
	int Ni, Nj;
	int ik, jk;
	Ni = speni.GetN();
	Nj = spenj.GetN();

	double Dij, Dikj, Dijk, Dikjk, minDij;
	std::vector<double> xpars(4);
	std::vector<double> ypars(4);

	xpars[0] = dpeni.GetX(speni[0]);
	xpars[2] = dpeni.GetX(speni[0]);
	ypars[0] = dpeni.GetY(speni[0]);
	ypars[2] = dpeni.GetY(speni[0]);

	xpars[1] = dpenj.GetX(spenj[0]);
	xpars[3] = dpenj.GetX(spenj[0]);
	ypars[1] = dpenj.GetY(spenj[0]);
	ypars[3] = dpenj.GetY(spenj[0]);
	
	Dij = (*myDTWCheckFunctions[icheck])(xpars, ypars);
	D.SetElem(0, 0, Dij);

	xpars[0] = dpeni.GetX(speni[0]);
	xpars[2] = dpeni.GetX(speni[0]);
	ypars[0] = dpeni.GetY(speni[0]);
	ypars[2] = dpeni.GetY(speni[0]);

	for (int j = 1; j < Nj; j++)
	{
		xpars[1] = dpenj.GetX(spenj[j]);
		xpars[3] = dpenj.GetX(spenj[j - 1]);
		ypars[1] = dpenj.GetY(spenj[j]);
		ypars[3] = dpenj.GetY(spenj[j - 1]);
		Dij = (*myDTWCheckFunctions[icheck])(xpars, ypars);
		D.SetElem(0, j, Dij);
	}

	xpars[1] = dpenj.GetX(spenj[0]);
	xpars[3] = dpenj.GetX(spenj[0]);
	ypars[1] = dpenj.GetY(spenj[0]);
	ypars[3] = dpenj.GetY(spenj[0]);

	for (int i = 1; i < Ni; i++)
	{
		xpars[0] = dpeni.GetX(speni[i]);
		xpars[2] = dpeni.GetX(speni[i - 1]);
		ypars[0] = dpeni.GetY(speni[i]);
		ypars[2] = dpeni.GetY(speni[i - 1]);

		Dij = (*myDTWCheckFunctions[icheck])(xpars, ypars);
		D.SetElem(i, 0, Dij);
	}

	for (int i = 1; i < Ni; i++)
	{ 
		xpars[0] = dpeni.GetX(speni[i]);
		xpars[2] = dpeni.GetX(speni[i - 1]);
		ypars[0] = dpeni.GetY(speni[i]);
		ypars[2] = dpeni.GetY(speni[i - 1]);
		for (int j = 1; j < Nj; j++)
		{
			xpars[1] = dpenj.GetX(spenj[j]);
			xpars[3] = dpenj.GetX(spenj[j - 1]);
			ypars[1] = dpenj.GetY(spenj[j]);
			ypars[3] = dpenj.GetY(spenj[j - 1]);

			Dij = (*myDTWCheckFunctions[icheck])(xpars, ypars);
			D.SetElem(i, j, Dij);
		}
	}
	for (int j = 1; j < Nj; j++)
	{
		minDij = D(0, j-1);
		for (int k = 2; k < DTWwindow + 1; k++)
		{
			jk = std::max(j - k, 0);
			Dijk = D(0, jk);
			minDij = std::min(minDij, Dijk);
		}
		Dij = D(0, j);
		D.SetElem(0, j, Dij + minDij);
	}

	for (int i = 1; i < Ni; i++)
	{
		minDij = D(i - 1, 0);
		for (int k = 2; k < DTWwindow + 1; k++)
		{
			ik = std::max(i - k, 0);
			Dikj = D(ik, 0);
			minDij = std::min(minDij, Dikj);
		}
		Dij = D(i, 0);
		D.SetElem(i, 0, Dij + minDij);

		for (int j = 1; j < Nj; j++)
		{
			minDij = D(i - 1, j - 1);
			for (int k = 1; k < DTWwindow + 1; k++)
			{
				ik = std::max(i - k, 0);
				jk = std::max(j - k, 0);
				Dikj = D(ik, j);
				Dijk = D(i, jk);
				Dikjk = D(ik, jk);
				minDij = std::min(minDij, std::min(Dikjk, std::min(Dikj, Dijk)));
			}
			Dij = D(i, j);
			D.SetElem(i, j, Dij + minDij);
		}
	}
	return true;
}

void SignChecker::DTW_TraceBack(Matrix<double>& D, SPoints& tranformpeni, SPoints& tranformpenj)
{
	int Ni, Nj, i, j;
	Ni = D.GetN();
	Nj = D.GetM();
	double Dij, Dikj, Dijk, Dikjk, minDij;
	int mini, minj;
	int ik, jk;
	double GC = D(Ni - 1, Nj - 1);
	//std::cout << " DTW_TraceBack\n " ;
	i = Ni - 1; j = Nj - 1;
	do {
		minDij = D(i, j);
		mini = -1;
		minj = -1;
		for (int k = 1; k < DTWwindow + 1; k++)
		{
			ik = std::max(i - k, 0);
			jk = std::max(j - k, 0);
			Dikj = D(ik, j);
			Dijk = D(i, jk);
			Dikjk = D(ik, jk);
			if (Dikj < minDij)
				if (Dikj < Dijk)
					if (Dikj < Dikjk)
					{
						minDij = Dikj;
						mini = ik;
						i = ik;
						minj = j;
					}
					else
					{
						minDij = Dikjk;
						mini = ik;
						i = ik;
						minj = jk;
						j = jk;
					}
				else
					if (Dijk < Dikjk)
					{
						minDij = Dijk;
						mini = i;
						minj = jk;
						j = jk;
					}
					else
					{
						minDij = Dikjk;
						mini = ik;
						i = ik;
						minj = jk;
						j = jk;
					}
			else
				if (Dijk < minDij)
					if (Dijk < Dikjk)
					{
						minDij = Dikj;
						mini = i;
						minj = jk;
						j = jk;
					}
					else
					{
						minDij = Dikjk;
						mini = ik;
						i = ik;
						minj = jk;
						j = jk;
					}
				else
					if (Dikjk < minDij)
					{
						minDij = Dikjk;
						mini = i;
						minj = jk;
						j = jk;
					}
		}
		if ((mini != -1) && (minj != -1))
		{
			//std::cout << mini << " - " << minj << "\n";
			tranformpeni.PushElem(mini);
			tranformpenj.PushElem(minj);
		}
		else break;

	} while (i != 0 || j != 0);
}

double SignChecker::DTW_CaclGlobalDeformation(Matrix<double>& D, SPoints& tranformpeni, SPoints& tranformpenj)
{
	int Nk = tranformpeni.GetN();
	double CG = 0;

	for (int i = 1; i < Nk - 1; i++)
		CG += abs(D(tranformpeni[Nk - i], tranformpenj[Nk - i]) - D(tranformpeni[Nk - (i+1)], tranformpenj[Nk - (i+1)]));

	return CG;
}

double SignChecker::DTW_CalcDetermination(SPoints& tpeni, DPoints& dpeni, SPoints& tpenj, DPoints& dpenj)
{
	int N = tpeni.GetN();

	double xci, yci, xcj, ycj;
	xci = tpeni.AvgDX(dpeni);
	yci = tpeni.AvgDY(dpeni);
	xcj = tpenj.AvgDX(dpenj);
	ycj = tpenj.AvgDY(dpenj);

	double xi, yi, xj, yj;
	double sumi = 0;
	double sumj = 0;
	double sumsqx = 0;
	double sumsqy = 0;
	
	for (int i = 0; i < N; i++)
	{
		xi = dpeni.GetX(tpeni[i]);
		yi = dpeni.GetY(tpeni[i]);
		xj = dpenj.GetX(tpenj[i]);
		yj = dpenj.GetY(tpenj[i]);
		sumi += MultDif(xi, yi, xci, yci);
		sumj += MultDif(xj, yj, xcj, ycj);
		sumsqx += SumSq(xi, xj, xci, xcj);
		sumsqy += SumSq(yi, yj, yci, ycj);
	}
	if ((sumi + sumj) == 0)
	{
		sumi = Eps;
		sumj = 0;
	}
	if (sumsqx == 0) sumsqx = Eps;
	if (sumsqy == 0) sumsqy = Eps;
	return (sumi + sumj)*(sumi + sumj) / (sumsqx*sumsqy);
}