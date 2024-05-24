#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <stdexcept>
#include <vector>

// шаблон матрицы

template <class TypeVect> class Matrix
{
protected:

	int N = 0; // колличесвто строк
	int M = 0; // колличество столбцов
	std::vector<TypeVect> data_;

public:

	Matrix() = default;
	Matrix(int n, int m)
	{
		Init(n, m);
	}
	~Matrix() = default;

	void Init(int n, int m)
	{
		if (n < 0 || m < 0)
			throw std::invalid_argument("Matrix dimensions must be non-negative");

		this->N = n;
		this->M = m;
		this->data_.assign(static_cast<std::size_t>(this->N * this->M), TypeVect{});
	}
	void Init(const Matrix<TypeVect>& Ldef, int defi, int defj)
	{
		int rows = Ldef.GetN() - std::min(Ldef.GetN(), std::abs(defi));
		int cols = Ldef.GetM() - std::min(Ldef.GetM(), std::abs(defj));

		if (rows == 1)
			rows = Ldef.GetN();
		if (cols == 1)
			cols = Ldef.GetM();

		Init(rows, cols);

		for (int i = 0; i < this->N; i++)
			for (int j = 0; j < this->M; j++)
				(*this)(i, j) = Ldef(i, j);
	}

	inline int GetN() const
	{
		return N;
	}
	inline int GetM() const
	{
		return M;
	}

	inline TypeVect &operator() (int i, int j)
	{
		return this->data_[Index(i, j)];
	}
	inline const TypeVect& operator() (int i, int j) const
	{
		return this->data_[Index(i, j)];
	}
	inline void SetElem(int i, int j, TypeVect x) {
		(*this)(i, j) = x;
	}

	// вычисление средних по строке, по столбцам, и в целом по матрице
	TypeVect AvgByJ(int j) const
	{
		EnsureColumn(j);
		double sum = 0.0;
		for (int i = 0; i< this->N; i++)
		{
			sum += (*this)(i, j);
		}

		return sum / double(this->N);
	}
	TypeVect AvgByI(int i) const
    {
		EnsureRow(i);
        double sum = 0.0;
        for (int j = 0; j < this->M; j++)
        {
            sum += (*this)(i, j);
        }

        return sum / double(this->M);
    }
    TypeVect AvgAll() const
    {
		EnsureNotEmpty();
        double sum = 0.0;
        for (int i = 0; i < this->N; i++)
            for (int j = 0; j < this->M; j++)
                sum += (*this)(i, j);

        return sum / double(this->M * this->N);
    }

	// вычисление оценки среднеквадратичного отклонения
	TypeVect SumSqByI(int i, TypeVect a) const
	{
		EnsureRow(i);
		double sum = 0.0;
		for (int j = 0; j < this->M; j++)
		{
			sum += ((*this)(i, j) - a)*((*this)(i, j) - a);
		}
		return sum / double(this->M);
	}

	TypeVect GetIndMinElem(int& mini, int& minj) const
	{
		EnsureNotEmpty();
		bool found = false;
		TypeVect minel{};
		for (int i = 0; i < this->N; i++)
			for (int j = 0; j < this->M; j++)
				if (((*this)(i, j) != TypeVect{}) && (!found || (*this)(i, j) < minel))
				{
					minel = (*this)(i, j);
					mini = i;
					minj = j;
					found = true;
				}
		if (!found)
		{
			mini = 0;
			minj = 0;
		}
		return minel;
	}
	TypeVect GetIndMaxElem(int& maxi, int& maxj) const
	{
		EnsureNotEmpty();
		TypeVect maxel = (*this)(0, 0);
		maxi = 0;
		maxj = 0;
		for (int i = 0; i < this->N; i++)
			for(int j = 0; j < this->M; j++)
				if ((*this)(i, j) > maxel)
				{
					maxel = (*this)(i, j);
					maxi = i;
					maxj = j;
				}
		return maxel;
	}

	// индекс минимального/максимального элемента в сечениях по столбцам или строкам
	int GetJMinElInI(int i) const
	{
		EnsureRow(i);
		int jmin = 0;
		bool found = false;
		TypeVect minel{};
		for (int j = 0; j < this->M; j++)
			if (((*this)(i, j) != TypeVect{}) && (!found || (*this)(i, j) < minel))
			{
				minel = (*this)(i, j);
				jmin = j;
				found = true;
			}
		return jmin;
	}
	int GetIMinElInJ(int j) const
	{
		EnsureColumn(j);
		int imin = 0;
		bool found = false;
		TypeVect minel{};
		for (int i = 0; i < this->N; i++)
			if (((*this)(i, j) != TypeVect{}) && (!found || (*this)(i, j) < minel))
			{
				minel = (*this)(i, j);
				imin = i;
				found = true;
			}
		return imin;
	}
	int GetJMaxElInI(int i) const
	{
		EnsureRow(i);
		int jmax = 0;
		TypeVect maxel = (*this)(i, 0);
		for (int j = 1; j < this->M; j++)
			if ((*this)(i, j) > maxel)
			{
				maxel = (*this)(i, j);
				jmax = j;
			}
		return jmax;
	}
	int GetIMaxElInJ(int j) const
	{
		EnsureColumn(j);
		int imax = 0;
		TypeVect maxel = (*this)(0, j);
		for (int i = 1; i < this->N; i++)
			if ((*this)(i, j) > maxel)
			{
				maxel = (*this)(i, j);
				imax = i;
			}
		return imax;
	}

private:
	std::size_t Index(int i, int j) const
	{
		if (i < 0 || i >= this->N || j < 0 || j >= this->M)
			throw std::out_of_range("Matrix index out of range");
		return static_cast<std::size_t>(i * this->M + j);
	}

	void EnsureNotEmpty() const
	{
		if (this->N == 0 || this->M == 0)
			throw std::out_of_range("Matrix is empty");
	}

	void EnsureRow(int i) const
	{
		if (i < 0 || i >= this->N || this->M == 0)
			throw std::out_of_range("Matrix row is out of range");
	}

	void EnsureColumn(int j) const
	{
		if (j < 0 || j >= this->M || this->N == 0)
			throw std::out_of_range("Matrix column is out of range");
	}
};
