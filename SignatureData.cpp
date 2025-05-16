#include "SignatureData.h"

#include <cmath>

double Eps = 1e-12;
double balance = 0.5;

double RatioXY(CTab* ipens, int i)
{
    return ipens[i].ratio_xy;
}

double AvgXY(CTab* ipens, int i)
{
    return ipens[i].avg_xy;
}

double SinXY(CTab* ipens, int i)
{
    return ipens[i].sin_xy;
}

double RatioTouches(CTab* ipens, int i)
{
    return ipens[i].ratio_touches;
}

double CityBlockMetric(std::vector<double>& x, std::vector<double>& y)
{
    return std::abs(x[0] - x[1]) + std::abs(y[0] - y[1]);
}

double EuclideMetric(std::vector<double>& x, std::vector<double>& y)
{
    const double dx = x[0] - x[1];
    const double dy = y[0] - y[1];
    return std::sqrt(dx * dx + dy * dy);
}

double VelocityMetric(std::vector<double>& x, std::vector<double>& y)
{
    const double vx0 = x[0] - x[2];
    const double vy0 = y[0] - y[2];
    const double vx1 = x[1] - x[3];
    const double vy1 = y[1] - y[3];
    const double v0 = std::sqrt(vx0 * vx0 + vy0 * vy0);
    const double v1 = std::sqrt(vx1 * vx1 + vy1 * vy1);
    return std::abs(v0 - v1);
}

double SinMetric(std::vector<double>& x, std::vector<double>& y)
{
    const double vx0 = x[0] - x[2];
    const double vy0 = y[0] - y[2];
    const double vx1 = x[1] - x[3];
    const double vy1 = y[1] - y[3];
    const double len0 = std::sqrt(vx0 * vx0 + vy0 * vy0);
    const double len1 = std::sqrt(vx1 * vx1 + vy1 * vy1);

    if (len0 <= Eps || len1 <= Eps)
        return 0.0;

    const double cross = vx0 * vy1 - vy0 * vx1;
    return std::abs(cross) / (len0 * len1);
}

double MultDif(double x, double y, double xc, double yc)
{
    return (x - xc) * (y - yc);
}

double SumSq(double x, double y, double xc, double yc)
{
    const double dx = x - xc;
    const double dy = y - yc;
    return dx * dx + dy * dy;
}
