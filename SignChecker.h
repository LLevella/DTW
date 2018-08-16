#pragma once

extern double balance;

class SignChecker
{
	// Колличество простых и DTW проверок
	int SimpleChecks; 
	int DTWChecks;

	// окно для построения оптимального пути
	int DTWwindow;

	// результаты проверок
	double DTWTestResult;
	double SimpleTestResult; 

	// функции простых проверок
	typedef double(*SimpleCheckFunctions)(CTab *ipens, int i);
	std::unique_ptr <SimpleCheckFunctions[]> mySimpleCheckFunctions;

	// функции DTW проверок
	typedef double(*DTWCheckFunctions)(std::vector<double> &x, std::vector<double> &y);
	std::unique_ptr <DTWCheckFunctions[]> myDTWCheckFunctions;

	// результаты вычисления оптимального пути 
	std::unique_ptr <Matrix<double>[]> DTWch;
	std::unique_ptr <Matrix<double>[]> DTW;

	// результаты глобальный путь трансформации
	std::unique_ptr <Matrix<double>[]> CGch;
	std::unique_ptr <Matrix<double>[]> CG;

	// оценки корреляции
	std::unique_ptr <Matrix<double>[]> CVch;

	// колличество точек в глобальном пути трансформации
	std::unique_ptr <Matrix<int>[]> Ncg;

	// результаты простых проверок
	std::unique_ptr <Matrix<double>[]> SimpleM;
	std::unique_ptr <Matrix<double>[]> SimpleMch;

public:
	SignChecker();
	
	inline double GetTestResult() 
	{
		return balance * this->SimpleTestResult + (1.0 - balance)*this->DTWTestResult;
	}
	
	// Инициализация проверок
	bool SimpleCheckForRandomForge(CTab *ipens, int npens);
	bool DTWCheckForSimpleForge(DPoints* dpens, SPoints* spens, int npen);
	
	//оценка близости подписи
	double CalcDifference(Matrix<double>& ratioM_withCheckedSign, Matrix<double>& ratioM_withoutCheckedSign, int npens);
	double DTWCalcDifference(int ncheck, int npens);

	void GenerateMatixeByIPen(int icheck, CTab *ipens, int npens);
	void GenerateMatrixByDPen(int icheck, DPoints* dpens, SPoints *spens, int npens);
	
	void DTW_Go(int icheck, DPoints* dpens, SPoints *spens, int i, int j);
	bool DTW_InitMatrix(int icheck, Matrix<double>& d, DPoints& dpeni, SPoints& speni, DPoints& dpenj, SPoints& spenj);
	void DTW_TraceBack(Matrix<double>& D, SPoints& tranformpeni, SPoints& tranformpenj);
	double DTW_CaclGlobalDeformation(Matrix<double>& D, SPoints& tranformpeni, SPoints& tranformpenj);
	double DTW_CalcDetermination(SPoints& tpeni, DPoints& dpeni, SPoints& tpenj, DPoints& dpenj);

};