#pragma once

#include <cstddef>
#include <initializer_list>
#include <utility>
#include <vector>

extern double Eps;
extern double balance;

struct CTab
{
    double ratio_xy = 0.0;
    double avg_xy = 0.0;
    double sin_xy = 0.0;
    double ratio_touches = 0.0;
};

struct DPoint
{
    double x = 0.0;
    double y = 0.0;
    double pressure = 0.0;
    double timestamp = 0.0;
    bool   pen_down = true;
    bool   has_pressure = false;
    bool   has_timestamp = false;
};

class DPoints
{
public:
    using PairPoint = std::pair<double, double>;

    DPoints() = default;

    DPoints(std::initializer_list<PairPoint> points)
    {
        points_.reserve(points.size());
        for (const auto& point : points)
        {
            DPoint dp;
            dp.x = point.first;
            dp.y = point.second;
            points_.push_back(dp);
        }
    }

    static DPoints FromExtended(std::initializer_list<DPoint> points)
    {
        DPoints result;
        result.points_.assign(points.begin(), points.end());
        return result;
    }

    void PushElem(double x, double y)
    {
        DPoint dp;
        dp.x = x;
        dp.y = y;
        points_.push_back(dp);
    }

    void AddPoint(double x, double y) { PushElem(x, y); }

    void AddPoint(double x, double y, double pressure)
    {
        DPoint dp;
        dp.x = x;
        dp.y = y;
        dp.pressure = pressure;
        dp.has_pressure = true;
        points_.push_back(dp);
    }

    void AddPoint(double x, double y, double pressure, double timestamp, bool pen_down = true)
    {
        DPoint dp;
        dp.x = x;
        dp.y = y;
        dp.pressure = pressure;
        dp.has_pressure = true;
        dp.timestamp = timestamp;
        dp.has_timestamp = true;
        dp.pen_down = pen_down;
        points_.push_back(dp);
    }

    void AddPoint(const DPoint& point) { points_.push_back(point); }

    int GetN() const { return static_cast<int>(points_.size()); }

    double GetX(int i) const { return points_.at(static_cast<std::size_t>(i)).x; }
    double GetY(int i) const { return points_.at(static_cast<std::size_t>(i)).y; }
    double GetPressure(int i) const { return points_.at(static_cast<std::size_t>(i)).pressure; }
    double GetTimestamp(int i) const { return points_.at(static_cast<std::size_t>(i)).timestamp; }
    bool   GetPenDown(int i) const { return points_.at(static_cast<std::size_t>(i)).pen_down; }
    bool   HasPressure(int i) const { return points_.at(static_cast<std::size_t>(i)).has_pressure; }
    bool   HasTimestamp(int i) const { return points_.at(static_cast<std::size_t>(i)).has_timestamp; }
    const DPoint& GetPoint(int i) const { return points_.at(static_cast<std::size_t>(i)); }

    bool AllHavePressure() const
    {
        for (const auto& p : points_)
            if (!p.has_pressure)
                return false;
        return !points_.empty();
    }

    bool AllHaveTimestamp() const
    {
        for (const auto& p : points_)
            if (!p.has_timestamp)
                return false;
        return !points_.empty();
    }

private:
    std::vector<DPoint> points_;
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
