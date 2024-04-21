#pragma once

#include <cstddef>
#include <initializer_list>
#include <utility>
#include <vector>

extern double agree_param;
extern double Eps;
extern double balance;

struct CTab
{
    double ratio_xy = 0.0;
    double avg_xy = 0.0;
    double sin_xy = 0.0;
    double ratio_touches = 0.0;
};

class DPoints
{
public:
    using Point = std::pair<double, double>;

    DPoints() = default;
    DPoints(std::initializer_list<Point> points) : points_(points) {}

    void PushElem(double x, double y)
    {
        points_.emplace_back(x, y);
    }

    void AddPoint(double x, double y)
    {
        PushElem(x, y);
    }

    int GetN() const
    {
        return static_cast<int>(points_.size());
    }

    double GetX(int i) const
    {
        return points_.at(static_cast<std::size_t>(i)).first;
    }

    double GetY(int i) const
    {
        return points_.at(static_cast<std::size_t>(i)).second;
    }

private:
    std::vector<Point> points_;
};

class SPoints
{
public:
    SPoints() = default;
    SPoints(std::initializer_list<int> indexes) : indexes_(indexes) {}

    static SPoints FullRange(int size)
    {
        SPoints result;
        for (int i = 0; i < size; ++i)
            result.PushElem(i);
        return result;
    }

    void PushElem(int index)
    {
        indexes_.push_back(index);
    }

    int GetN() const
    {
        return static_cast<int>(indexes_.size());
    }

    int operator[](int i) const
    {
        return indexes_.at(static_cast<std::size_t>(i));
    }

    int& operator[](int i)
    {
        return indexes_.at(static_cast<std::size_t>(i));
    }

    double AvgDX(const DPoints& points) const
    {
        if (indexes_.empty())
            return 0.0;

        double sum = 0.0;
        for (const int index : indexes_)
            sum += points.GetX(index);
        return sum / static_cast<double>(indexes_.size());
    }

    double AvgDY(const DPoints& points) const
    {
        if (indexes_.empty())
            return 0.0;

        double sum = 0.0;
        for (const int index : indexes_)
            sum += points.GetY(index);
        return sum / static_cast<double>(indexes_.size());
    }

private:
    std::vector<int> indexes_;
};

double RatioXY(CTab* ipens, int i);
double AvgXY(CTab* ipens, int i);
double SinXY(CTab* ipens, int i);
double RatioTouches(CTab* ipens, int i);

double CityBlockMetric(std::vector<double>& x, std::vector<double>& y);
double EuclideMetric(std::vector<double>& x, std::vector<double>& y);
double VelocityMetric(std::vector<double>& x, std::vector<double>& y);
double SinMetric(std::vector<double>& x, std::vector<double>& y);

double MultDif(double x, double y, double xc, double yc);
double SumSq(double x, double y, double xc, double yc);
