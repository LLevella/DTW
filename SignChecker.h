#pragma once

extern double balance;

class SignChecker
{
	// ����������� ������� � DTW ��������
	int SimpleChecks; 
	int DTWChecks;

	// ���� ��� ���������� ������������ ����
	int DTWwindow;

	// ���������� ��������
	double DTWTestResult;
	double SimpleTestResult; 

	// ������� ������� ��������
	typedef double(*SimpleCheckFunctions)(CTab *ipens, int i);
	std::unique_ptr <SimpleCheckFunctions[]> mySimpleCheckFunctions;

	// ������� DTW ��������
	typedef double(*DTWCheckFunctions)(std::vector<double> &x, std::vector<double> &y);
	std::unique_ptr <DTWCheckFunctions[]> myDTWCheckFunctions;

	// ���������� ���������� ������������ ���� 
	std::unique_ptr <Matrix<double>[]> DTWch;
	std::unique_ptr <Matrix<double>[]> DTW;

	// ���������� ���������� ���� �������������
	std::unique_ptr <Matrix<double>[]> CGch;
	std::unique_ptr <Matrix<double>[]> CG;

	// ������ ����������
	std::unique_ptr <Matrix<double>[]> CVch;

	// ����������� ����� � ���������� ���� �������������
	std::unique_ptr <Matrix<int>[]> Ncg;

	// ���������� ������� ��������
	std::unique_ptr <Matrix<double>[]> SimpleM;
	std::unique_ptr <Matrix<double>[]> SimpleMch;

public:
	SignChecker();
	
	inline double GetTestResult() 
	{
		return balance * this->SimpleTestResult + (1.0 - balance)*this->DTWTestResult;
	}
	
	// ������������� ��������
	bool SimpleCheckForRandomForge(CTab *ipens, int npens);
	bool DTWCheckForSimpleForge(DPoints* dpens, SPoints* spens, int npen);
	
	//������ �������� �������
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