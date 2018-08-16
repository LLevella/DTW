#pragma once

// шаблон матрицы

template <class TypeVect> class Matrix
{
protected:

	int N; // колличесвто строк
	int M; // колличество столбцов
	
	TypeVect** A;
	
public:

	Matrix() {};
	Matrix(int n, int m)
	{
		Init(n, m);
	}
	~Matrix()
	{
		for (int i = 0; i < this->N; i++)
			delete[] this->A[i];
		delete[] A;
	}
	
	void Init(int n, int m)
	{
		this->N = n;
		this->M = m;
		this->A = new TypeVect*[this->N];
		for (int i = 0; i < this->N; i++)
			this->A[i] = new TypeVect[this->M];

		for (int i = 0; i < this->N; i++)
			for (int j = 0; j < this->M; j++)
				this->A[i][j] = 0;
	}
	void Init(Matrix<TypeVect>& Ldef, int defi, int defj)
	{
		this->N = Ldef.GetN() - std::min(Ldef.GetN(), abs(defi));
		this->M = Ldef.GetM() - std::min(Ldef.GetM(), abs(defj));
		
		if (this->N == 1)
			this->N = Ldef.GetN();
		if (this->M == 1)
			this->M = Ldef.GetM();

		this->A = new TypeVect*[this->N];
		for (int i = 0; i < this->N; i++)
			this->A[i] = new TypeVect[this->M];

		for (int i = 0; i < this->N; i++)
			for (int j = 0; j < this->M; j++)
				this->A[i][j] = Ldef(i, j);
	}
	
	inline int GetN()
	{
		return N;
	}
	inline int GetM()
	{
		return M;
	}
	
	inline TypeVect &operator() (int i, int j)
	{
		return this->A[i][j];
	}
	inline void SetElem(int i, int j, TypeVect x) {
		this->A[i][j] = x;
	}
	
	// вычисление средних по строке, по столбцам, и в целом по матрице 
	TypeVect AvgByJ(int j)
	{
		double sum = 0.0;
		for (int i = 0; i< this->N; i++)
		{
			sum += this->A[i][j];
		}

		return sum / double(this->N);
	}
	TypeVect AvgByI(int i)
    {
        double sum = 0.0;
        for (int j = 0; j < this->M; j++)
        {
            sum += this->A[i][j];
        }
                
        return sum / double(this->M);
    }
    TypeVect AvgAll()
    {
        double sum = 0.0;
        for (int i = 0; i < this->N; i++)
            for (int j = 0; j < this->M; j++)
                sum += this->A[i][j];
    
        return sum / double(this->M * this->N);
    }

	// вычисление оценки среднеквадратичного отклонения 
	TypeVect SumSqByI(int i, TypeVect a)
	{
		double sum = 0.0;
		for (int j = 0; j < this->M; j++)
		{
			sum += (this->A[i][j] - a)*(this->A[i][j] - a);
		}
		return sum / double(this->M);
	}

	TypeVect GetIndMinElem(int& mini, int& minj)
	{
		TypeVect minel = this->A[0][0];
		for (int i = 1; i < this->N; i++)
			for (int j = 0; j < this->M; j++)
				if ((this->A[i][j] < minel) && (this->A[i][j] != 0))
				{
					minel = this->A[i][j];
					mini = i;
					minj = j;
				}
		return minel;
	}
	TypeVect GetIndMaxElem(int& maxi, int& maxj)
	{
		TypeVect maxel = this->A[0][0];
		for (int i = 1; i < this->N; i++)
			for(int j = 0; j < this->M; j++)
				if (this->A[i][j] > maxel)
				{
					maxel = this->A[i][j];
					maxi = i;
					maxj = j;
				}
		return maxel;
	}

	// индекс минимального/максимального элемента в сечениях по столбцам или строкам
	int GetJMinElInI(int i)
	{
		int jmin = 0;
		TypeVect minel = this->A[i][0];
		for (int j = 1; j < this->M; j++)
			if ((this->A[i][j] < minel) && (this->A[i][j] != 0))
			{
				minel = this->A[i][j];
				jmin = j;
			}
		return jmin;
	}
	int GetIMinElInJ(int j)
	{
		int imin = 0;
		TypeVect minel = this->A[0][j];
		for (int i = 1; i < this->N; i++)
			if ((this->A[i][j] < minel) && (this->A[i][j] != 0))
			{
				minel = this->A[i][j];
				imin = i;
			}
		return imin;
	}
	int GetJMaxElInI(int i)
	{
		int jmax = 0;
		TypeVect maxel = this->A[i][0];
		for (int j = 1; j < this->M; j++)
			if (this->A[i][j] > maxel)
			{
				maxel = this->A[i][j];
				jmax = j;
			}
		return jmax;
	}
	int GetIMaxElInJ(int j)
	{
		int imax = 0;
		TypeVect maxel = this->A[0][j];
		for (int i = 1; i < this->N; i++)
			if (this->A[i][j] > maxel)
			{
				maxel = this->A[i][j];
				imax = i;
			}
		return imax;
	}
};