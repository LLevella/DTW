#include "SignChecker.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <vector>

extern double Eps;

namespace
{
constexpr double Pi = 3.141592653589793238462643383279502884;

struct SampledPen
{
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> pressure;
    std::vector<double> timestamp;
    std::vector<double> tangent;        // углы касательной [-π, π]
    std::vector<double> curvature;
    std::vector<double> seglen;         // длина сегмента к предыдущей точке
    std::vector<double> velocity;       // seglen / Δt (если есть timestamp)
    std::vector<double> pseudopressure; // оценка скорости через плотность точек
    bool has_pressure = false;
    bool has_timestamp = false;
};

struct FeatureNeeds
{
    bool pressure = false;
    bool timestamp = false;
    bool tangent = false;
    bool curvature = false;
    bool seglen = false;
    bool velocity = false;
    bool pseudopressure = false;
};

FeatureNeeds DetermineFeatureNeeds(const SignCheckerConfig& cfg, bool includeShapeFeatures = false)
{
    FeatureNeeds needs;
    using DC = SignCheckerConfig::DtwChannel;

    for (const DC channel : cfg.dtw_channels)
    {
        switch (channel)
        {
        case DC::CityBlock:
        case DC::Direction:
            break;
        case DC::SegmentLength:
            needs.seglen = true;
            break;
        case DC::TangentAngle:
            needs.tangent = true;
            break;
        case DC::Curvature:
            needs.curvature = true;
            break;
        case DC::Pseudopressure:
            needs.seglen = true;
            needs.pseudopressure = true;
            break;
        case DC::Pressure:
            needs.pressure = true;
            break;
        case DC::Velocity:
            needs.timestamp = true;
            needs.seglen = true;
            needs.velocity = true;
            break;
        }
    }

    if (includeShapeFeatures)
        needs.curvature = true;

    return needs;
}

double RelativeDifference(double left, double right)
{
    const double denominator = std::max(std::max(std::abs(left), std::abs(right)), Eps);
    if (denominator <= Eps)
        return 0.0;
    return std::abs(left - right) / denominator;
}

double Median(std::vector<double> values)
{
    if (values.empty())
        return 0.0;

    const std::size_t n = values.size();
    const std::size_t mid = n / 2;
    std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(mid), values.end());
    if (n % 2 == 1)
        return values[mid];

    const double upper = values[mid];
    std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(mid - 1), values.begin() + static_cast<std::ptrdiff_t>(mid));
    return 0.5 * (values[mid - 1] + upper);
}

double RowRmsPrefix(const Matrix<double>& matrix, int row, int count)
{
    if (count <= 0)
        return 0.0;
    double sumSquares = 0.0;
    for (int j = 0; j < count; ++j)
    {
        const double value = matrix(row, j);
        sumSquares += value * value;
    }
    return std::sqrt(sumSquares / static_cast<double>(count));
}

double RowMedianPrefix(const Matrix<double>& matrix, int row, int count)
{
    if (count <= 0)
        return 0.0;
    std::vector<double> values;
    values.reserve(static_cast<std::size_t>(count));
    for (int j = 0; j < count; ++j)
        values.push_back(std::abs(matrix(row, j)));
    return Median(std::move(values));
}

double PairwiseRms(const Matrix<double>& matrix, int count)
{
    if (count < 2)
        return 0.0;
    double sumSquares = 0.0;
    int pairs = 0;
    for (int i = 1; i < count; ++i)
        for (int j = 0; j < i; ++j)
        {
            const double value = matrix(i, j);
            sumSquares += value * value;
            ++pairs;
        }
    return pairs == 0 ? 0.0 : std::sqrt(sumSquares / static_cast<double>(pairs));
}

double PairwiseMedian(const Matrix<double>& matrix, int count)
{
    if (count < 2)
        return 0.0;
    std::vector<double> values;
    values.reserve(static_cast<std::size_t>(count) * static_cast<std::size_t>(count - 1) / 2);
    for (int i = 1; i < count; ++i)
        for (int j = 0; j < i; ++j)
            values.push_back(std::abs(matrix(i, j)));
    return Median(std::move(values));
}

void SmoothInplace(std::vector<double>& v, int window)
{
    if (window <= 1 || v.size() < 2)
        return;

    const int half = window / 2;
    std::vector<double> out(v.size());
    double sum = 0.0;
    int begin = 0;
    int end = 0;

    for (std::size_t i = 0; i < v.size(); ++i)
    {
        const int wantedBegin = std::max(0, static_cast<int>(i) - half);
        const int wantedEnd = std::min(static_cast<int>(v.size()) - 1, static_cast<int>(i) + half);

        while (begin < wantedBegin)
        {
            sum -= v[static_cast<std::size_t>(begin)];
            ++begin;
        }
        while (end <= wantedEnd)
        {
            sum += v[static_cast<std::size_t>(end)];
            ++end;
        }

        const int count = end - begin;
        out[i] = count > 0 ? sum / static_cast<double>(count) : v[i];
    }
    v.swap(out);
}

void ArclengthResample(std::vector<double>& x,
                       std::vector<double>& y,
                       std::vector<double>* pressure,
                       std::vector<double>* timestamp,
                       int target)
{
    if (target < 2 || x.size() < 2 || x.size() != y.size())
        return;

    const std::size_t n = x.size();
    std::vector<double> s(n, 0.0);
    for (std::size_t i = 1; i < n; ++i)
    {
        const double dx = x[i] - x[i - 1];
        const double dy = y[i] - y[i - 1];
        s[i] = s[i - 1] + std::sqrt(dx * dx + dy * dy);
    }
    const double total = s.back();
    if (total <= Eps)
        return;

    const std::size_t targetSize = static_cast<std::size_t>(target);
    std::vector<double> nx(targetSize);
    std::vector<double> ny(targetSize);
    std::vector<double> np;
    std::vector<double> nt;
    if (pressure)
        np.resize(targetSize);
    if (timestamp)
        nt.resize(targetSize);

    std::size_t idx = 1;
    for (std::size_t k = 0; k < targetSize; ++k)
    {
        const double target_s = total * static_cast<double>(k) / static_cast<double>(target - 1);
        while (idx < n - 1 && s[idx] < target_s)
            ++idx;
        const double seg = s[idx] - s[idx - 1];
        const double a = seg <= Eps ? 0.0 : (target_s - s[idx - 1]) / seg;
        nx[k] = x[idx - 1] + a * (x[idx] - x[idx - 1]);
        ny[k] = y[idx - 1] + a * (y[idx] - y[idx - 1]);
        if (pressure)
            np[k] = (*pressure)[idx - 1] + a * ((*pressure)[idx] - (*pressure)[idx - 1]);
        if (timestamp)
            nt[k] = (*timestamp)[idx - 1] + a * ((*timestamp)[idx] - (*timestamp)[idx - 1]);
    }

    x.swap(nx);
    y.swap(ny);
    if (pressure)
        pressure->swap(np);
    if (timestamp)
        timestamp->swap(nt);
}

void PcaRotate(std::vector<double>& x, std::vector<double>& y)
{
    if (x.size() < 2)
        return;
    const double n = static_cast<double>(x.size());
    double mx = 0.0;
    double my = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        mx += x[i];
        my += y[i];
    }
    mx /= n;
    my /= n;

    double cxx = 0.0;
    double cyy = 0.0;
    double cxy = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        const double dx = x[i] - mx;
        const double dy = y[i] - my;
        cxx += dx * dx;
        cyy += dy * dy;
        cxy += dx * dy;
    }
    cxx /= n;
    cyy /= n;
    cxy /= n;

    const double trace = cxx + cyy;
    const double det = cxx * cyy - cxy * cxy;
    const double disc = std::max(0.0, trace * trace / 4.0 - det);
    const double lam = trace / 2.0 + std::sqrt(disc);

    double vx;
    double vy;
    if (std::abs(cxy) > Eps)
    {
        vx = lam - cyy;
        vy = cxy;
    }
    else
    {
        vx = (cxx >= cyy) ? 1.0 : 0.0;
        vy = (cxx >= cyy) ? 0.0 : 1.0;
    }
    const double len = std::sqrt(vx * vx + vy * vy);
    if (len <= Eps)
        return;
    vx /= len;
    vy /= len;

    for (std::size_t i = 0; i < x.size(); ++i)
    {
        const double dx = x[i] - mx;
        const double dy = y[i] - my;
        x[i] = mx + vx * dx + vy * dy;
        y[i] = my - vy * dx + vx * dy;
    }
}

void ZScoreNormalize(std::vector<double>& x, std::vector<double>& y)
{
    if (x.empty())
        return;
    const double n = static_cast<double>(x.size());
    double mx = 0.0;
    double my = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        mx += x[i];
        my += y[i];
    }
    mx /= n;
    my /= n;

    double vx = 0.0;
    double vy = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        vx += (x[i] - mx) * (x[i] - mx);
        vy += (y[i] - my) * (y[i] - my);
    }
    const double sx = std::sqrt(std::max(vx / n, Eps));
    const double sy = std::sqrt(std::max(vy / n, Eps));

    for (std::size_t i = 0; i < x.size(); ++i)
    {
        x[i] = (x[i] - mx) / sx;
        y[i] = (y[i] - my) / sy;
    }
}

void CenterAndDiagScale(std::vector<double>& x, std::vector<double>& y)
{
    if (x.empty())
        return;
    double minX = std::numeric_limits<double>::infinity();
    double minY = std::numeric_limits<double>::infinity();
    double maxX = -std::numeric_limits<double>::infinity();
    double maxY = -std::numeric_limits<double>::infinity();
    double sumX = 0.0;
    double sumY = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        sumX += x[i];
        sumY += y[i];
        minX = std::min(minX, x[i]);
        minY = std::min(minY, y[i]);
        maxX = std::max(maxX, x[i]);
        maxY = std::max(maxY, y[i]);
    }
    const double cx = sumX / static_cast<double>(x.size());
    const double cy = sumY / static_cast<double>(y.size());
    const double w = maxX - minX;
    const double h = maxY - minY;
    const double scale = std::max(std::sqrt(w * w + h * h), Eps);

    for (std::size_t i = 0; i < x.size(); ++i)
    {
        x[i] = (x[i] - cx) / scale;
        y[i] = (y[i] - cy) / scale;
    }
}

void ComputeDerivedFeatures(SampledPen& pen, const FeatureNeeds& needs)
{
    const std::size_t n = pen.x.size();
    const bool needsSeglen = needs.seglen || needs.velocity || needs.pseudopressure;
    if (needs.tangent)
        pen.tangent.assign(n, 0.0);
    if (needs.curvature)
        pen.curvature.assign(n, 0.0);
    if (needsSeglen)
        pen.seglen.assign(n, 0.0);
    if (needs.velocity && pen.has_timestamp)
        pen.velocity.assign(n, 0.0);
    if (needs.pseudopressure)
        pen.pseudopressure.assign(n, 0.0);

    if (n < 2)
        return;

    if (needsSeglen)
    {
        for (std::size_t i = 1; i < n; ++i)
        {
            const double dx = pen.x[i] - pen.x[i - 1];
            const double dy = pen.y[i] - pen.y[i - 1];
            pen.seglen[i] = std::sqrt(dx * dx + dy * dy);
        }
        pen.seglen[0] = pen.seglen[1];
    }

    if (needs.velocity && pen.has_timestamp)
    {
        for (std::size_t i = 1; i < n; ++i)
        {
            const double dt = pen.timestamp[i] - pen.timestamp[i - 1];
            pen.velocity[i] = dt > Eps ? pen.seglen[i] / dt : 0.0;
        }
        pen.velocity[0] = pen.velocity[1];
    }

    if (needs.pseudopressure)
    {
        std::vector<double> seglenSorted(pen.seglen);
        const double medSeg = Median(std::move(seglenSorted));
        const double scale = std::max(medSeg, Eps);
        for (std::size_t i = 0; i < n; ++i)
            pen.pseudopressure[i] = scale / std::max(pen.seglen[i], Eps);
    }

    if (needs.tangent)
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            const std::size_t prev = (i == 0) ? 0 : i - 1;
            const std::size_t next = (i + 1 < n) ? i + 1 : i;
            const double dx = pen.x[next] - pen.x[prev];
            const double dy = pen.y[next] - pen.y[prev];
            pen.tangent[i] = std::atan2(dy, dx);
        }
    }

    if (needs.curvature)
    {
        for (std::size_t i = 1; i + 1 < n; ++i)
        {
            const double xp = pen.x[i + 1] - pen.x[i - 1];
            const double yp = pen.y[i + 1] - pen.y[i - 1];
            const double xpp = pen.x[i + 1] - 2.0 * pen.x[i] + pen.x[i - 1];
            const double ypp = pen.y[i + 1] - 2.0 * pen.y[i] + pen.y[i - 1];
            const double speedSq = xp * xp + yp * yp;
            if (speedSq <= Eps)
            {
                pen.curvature[i] = 0.0;
                continue;
            }
            const double denom = speedSq * std::sqrt(speedSq);
            pen.curvature[i] = (xp * ypp - yp * xpp) / denom;
        }
    }
}

SampledPen BuildSampledPen(DPoints& points, SPoints& samples, const SignCheckerConfig& cfg, const FeatureNeeds& needs)
{
    SampledPen pen;
    const int count = samples.GetN();
    if (count <= 0)
        return pen;

    pen.x.reserve(static_cast<std::size_t>(count));
    pen.y.reserve(static_cast<std::size_t>(count));

    bool hasP = true;
    bool hasT = true;
    for (int i = 0; i < count; ++i)
    {
        const int idx = samples[i];
        pen.x.push_back(points.GetX(idx));
        pen.y.push_back(points.GetY(idx));
        if (needs.pressure && !points.HasPressure(idx))
            hasP = false;
        if (needs.timestamp && !points.HasTimestamp(idx))
            hasT = false;
    }
    pen.has_pressure = needs.pressure && hasP;
    pen.has_timestamp = needs.timestamp && hasT;
    if (pen.has_pressure)
    {
        pen.pressure.reserve(static_cast<std::size_t>(count));
        for (int i = 0; i < count; ++i)
            pen.pressure.push_back(points.GetPressure(samples[i]));
    }
    if (pen.has_timestamp)
    {
        pen.timestamp.reserve(static_cast<std::size_t>(count));
        for (int i = 0; i < count; ++i)
            pen.timestamp.push_back(points.GetTimestamp(samples[i]));
    }

    if (cfg.smooth)
    {
        SmoothInplace(pen.x, cfg.smooth_window);
        SmoothInplace(pen.y, cfg.smooth_window);
        if (pen.has_pressure)
            SmoothInplace(pen.pressure, cfg.smooth_window);
    }

    if (cfg.arclength_resample && cfg.resample_points >= 2)
    {
        std::vector<double>* pp = pen.has_pressure ? &pen.pressure : nullptr;
        std::vector<double>* tp = pen.has_timestamp ? &pen.timestamp : nullptr;
        ArclengthResample(pen.x, pen.y, pp, tp, cfg.resample_points);
    }

    if (cfg.pca_rotate)
        PcaRotate(pen.x, pen.y);

    if (cfg.zscore_normalize)
        ZScoreNormalize(pen.x, pen.y);
    else
        CenterAndDiagScale(pen.x, pen.y);

    ComputeDerivedFeatures(pen, needs);
    return pen;
}

double AngularDistance(double a, double b)
{
    double diff = std::abs(a - b);
    while (diff > Pi)
        diff = std::abs(diff - 2.0 * Pi);
    return diff / Pi; // [0, 1]
}

double LocalDTWCost(SignCheckerConfig::DtwChannel channel,
                    const SampledPen& left, int i,
                    const SampledPen& right, int j)
{
    const std::size_t li = static_cast<std::size_t>(i);
    const std::size_t lj = static_cast<std::size_t>(j);

    using DC = SignCheckerConfig::DtwChannel;
    switch (channel)
    {
    case DC::CityBlock:
        return std::abs(left.x[li] - right.x[lj]) + std::abs(left.y[li] - right.y[lj]);

    case DC::Direction:
    {
        if (i <= 0 || j <= 0)
            return 0.0;
        const double lvx = left.x[li] - left.x[li - 1];
        const double lvy = left.y[li] - left.y[li - 1];
        const double rvx = right.x[lj] - right.x[lj - 1];
        const double rvy = right.y[lj] - right.y[lj - 1];
        const double leftLen = std::sqrt(lvx * lvx + lvy * lvy);
        const double rightLen = std::sqrt(rvx * rvx + rvy * rvy);
        if (leftLen <= Eps || rightLen <= Eps)
            return 0.0;
        const double cross = lvx * rvy - lvy * rvx;
        return std::abs(cross) / (leftLen * rightLen);
    }

    case DC::SegmentLength:
        if (left.seglen.empty() || right.seglen.empty())
            return 0.0;
        return std::abs(left.seglen[li] - right.seglen[lj]);

    case DC::TangentAngle:
        if (left.tangent.empty() || right.tangent.empty())
            return 0.0;
        return AngularDistance(left.tangent[li], right.tangent[lj]);

    case DC::Curvature:
        if (left.curvature.empty() || right.curvature.empty())
            return 0.0;
        return std::abs(left.curvature[li] - right.curvature[lj]);

    case DC::Pseudopressure:
        if (left.pseudopressure.empty() || right.pseudopressure.empty())
            return 0.0;
        return std::abs(left.pseudopressure[li] - right.pseudopressure[lj]);

    case DC::Pressure:
        if (!left.has_pressure || !right.has_pressure)
            return 0.0;
        return std::abs(left.pressure[li] - right.pressure[lj]);

    case DC::Velocity:
        if (!left.has_timestamp || !right.has_timestamp || left.velocity.empty() || right.velocity.empty())
            return 0.0;
        return std::abs(left.velocity[li] - right.velocity[lj]);
    }
    return 0.0;
}

bool OutsideSakoeChiba(int i, int j, int Ni, int Nj, int band)
{
    if (band < 0 || Ni <= 1 || Nj <= 1)
        return false;
    const double iScaled = static_cast<double>(i) * static_cast<double>(Nj - 1) / static_cast<double>(Ni - 1);
    return std::abs(iScaled - static_cast<double>(j)) > static_cast<double>(band);
}

bool FillDTWMatrix(Matrix<double>& D,
                   SignCheckerConfig::DtwChannel channel,
                   const SampledPen& left,
                   const SampledPen& right,
                   int window,
                   int band)
{
    const int Ni = static_cast<int>(left.x.size());
    const int Nj = static_cast<int>(right.x.size());
    if (Ni == 0 || Nj == 0 || D.GetN() != Ni || D.GetM() != Nj)
        return false;

    for (int i = 0; i < Ni; i++)
    {
        for (int j = 0; j < Nj; j++)
        {
            if (OutsideSakoeChiba(i, j, Ni, Nj, band))
            {
                D.SetElem(i, j, std::numeric_limits<double>::infinity());
                continue;
            }

            double bestPrevious = std::numeric_limits<double>::infinity();
            for (int k = 1; k <= window; k++)
            {
                const int ik = std::max(i - k, 0);
                const int jk = std::max(j - k, 0);

                if (ik != i || jk != j)
                    bestPrevious = std::min(bestPrevious, D(ik, jk));
                if (ik != i)
                    bestPrevious = std::min(bestPrevious, D(ik, j));
                if (jk != j)
                    bestPrevious = std::min(bestPrevious, D(i, jk));
            }

            if (i == 0 && j == 0)
                bestPrevious = 0.0;

            D.SetElem(i, j, LocalDTWCost(channel, left, i, right, j) + bestPrevious);
        }
    }

    return true;
}

void TraceBackMatrix(const Matrix<double>& D, int window, SPoints& transformedI, SPoints& transformedJ)
{
    const int Ni = D.GetN();
    const int Nj = D.GetM();
    if (Ni == 0 || Nj == 0)
        return;

    int i = Ni - 1;
    int j = Nj - 1;
    window = std::max(1, window);
    transformedI.PushElem(i);
    transformedJ.PushElem(j);

    while (i != 0 || j != 0)
    {
        double minDij = std::numeric_limits<double>::infinity();
        int mini = i;
        int minj = j;

        const auto tryCandidate = [&](int candidateI, int candidateJ)
        {
            if (candidateI == i && candidateJ == j)
                return;

            const double value = D(candidateI, candidateJ);
            if (value < minDij)
            {
                minDij = value;
                mini = candidateI;
                minj = candidateJ;
            }
        };

        for (int k = 1; k <= window; k++)
        {
            const int ik = std::max(i - k, 0);
            const int jk = std::max(j - k, 0);

            tryCandidate(ik, jk);
            tryCandidate(ik, j);
            tryCandidate(i, jk);
        }

        if (mini == i && minj == j)
            break;

        i = mini;
        j = minj;
        transformedI.PushElem(i);
        transformedJ.PushElem(j);
    }
}

double GlobalDeformation(const Matrix<double>& D, SPoints& transformedI, SPoints& transformedJ)
{
    const int count = transformedI.GetN();
    double result = 0.0;
    for (int i = 1; i < count; i++)
        result += std::abs(D(transformedI[i - 1], transformedJ[i - 1]) - D(transformedI[i], transformedJ[i]));
    return result;
}

double DeterminationFromSamples(SPoints& pathI, const SampledPen& left, SPoints& pathJ, const SampledPen& right)
{
    const int count = pathI.GetN();
    if (count == 0 || pathJ.GetN() != count)
        return 0.0;

    double xci = 0.0;
    double yci = 0.0;
    double xcj = 0.0;
    double ycj = 0.0;
    for (int k = 0; k < count; ++k)
    {
        const std::size_t i = static_cast<std::size_t>(pathI[k]);
        const std::size_t j = static_cast<std::size_t>(pathJ[k]);
        xci += left.x[i];
        yci += left.y[i];
        xcj += right.x[j];
        ycj += right.y[j];
    }
    xci /= static_cast<double>(count);
    yci /= static_cast<double>(count);
    xcj /= static_cast<double>(count);
    ycj /= static_cast<double>(count);

    double sumi = 0.0;
    double sumj = 0.0;
    double sumsqx = 0.0;
    double sumsqy = 0.0;
    for (int k = 0; k < count; ++k)
    {
        const std::size_t i = static_cast<std::size_t>(pathI[k]);
        const std::size_t j = static_cast<std::size_t>(pathJ[k]);
        const double xi = left.x[i];
        const double yi = left.y[i];
        const double xj = right.x[j];
        const double yj = right.y[j];

        sumi += (xi - xci) * (yi - yci);
        sumj += (xj - xcj) * (yj - ycj);
        sumsqx += (xi - xci) * (xi - xci) + (xj - xcj) * (xj - xcj);
        sumsqy += (yi - yci) * (yi - yci) + (yj - ycj) * (yj - ycj);
    }

    if ((sumi + sumj) == 0)
    {
        sumi = Eps;
        sumj = 0;
    }
    if (sumsqx == 0)
        sumsqx = Eps;
    if (sumsqy == 0)
        sumsqy = Eps;

    const double determination = (sumi + sumj) * (sumi + sumj) / (sumsqx * sumsqy);
    return std::max(0.0, std::min(1.0, determination));
}

void StoreDTWPair(const SampledPen& left,
                  const SampledPen& right,
                  Matrix<double>& dtw,
                  Matrix<double>& cg,
                  Matrix<double>& cv,
                  Matrix<int>& ncg,
                  SignCheckerConfig::DtwChannel channel,
                  int window,
                  int band,
                  int i,
                  int j)
{
    const int Ni = static_cast<int>(left.x.size());
    const int Nj = static_cast<int>(right.x.size());
    if (Ni == 0 || Nj == 0)
        return;

    Matrix<double> D(Ni, Nj);
    if (!FillDTWMatrix(D, channel, left, right, window, band))
        return;

    SPoints pathI;
    SPoints pathJ;
    TraceBackMatrix(D, window, pathI, pathJ);

    const int pathLength = std::max(pathI.GetN(), 1);
    double finalCost = D(Ni - 1, Nj - 1);
    if (!std::isfinite(finalCost))
        finalCost = 1.0;

    const double normalizedCost = finalCost / static_cast<double>(pathLength);
    dtw.SetElem(i, j, normalizedCost);
    dtw.SetElem(j, i, normalizedCost);
    ncg.SetElem(i, j, pathI.GetN());
    ncg.SetElem(j, i, pathI.GetN());

    const double globalDeformation = GlobalDeformation(D, pathI, pathJ) / static_cast<double>(pathLength);
    cg.SetElem(i, j, globalDeformation);
    cg.SetElem(j, i, globalDeformation);

    const double determination = DeterminationFromSamples(pathI, left, pathJ, right);
    cv.SetElem(i, j, determination);
    cv.SetElem(j, i, determination);
}

struct ShapeFeatures
{
    double path_over_diag = 0.0;
    double aspect = 0.0;
    double curvature_mean = 0.0;
};

ShapeFeatures ComputeShape(DPoints& points, SPoints& samples, const SignCheckerConfig& cfg)
{
    ShapeFeatures features;
    const int n = samples.GetN();
    if (n <= 0)
        return features;

    double minX = std::numeric_limits<double>::infinity();
    double minY = std::numeric_limits<double>::infinity();
    double maxX = -std::numeric_limits<double>::infinity();
    double maxY = -std::numeric_limits<double>::infinity();
    double pathLen = 0.0;
    double prevX = points.GetX(samples[0]);
    double prevY = points.GetY(samples[0]);
    for (int i = 0; i < n; ++i)
    {
        const double x = points.GetX(samples[i]);
        const double y = points.GetY(samples[i]);
        minX = std::min(minX, x);
        minY = std::min(minY, y);
        maxX = std::max(maxX, x);
        maxY = std::max(maxY, y);
        if (i > 0)
        {
            const double dx = x - prevX;
            const double dy = y - prevY;
            pathLen += std::sqrt(dx * dx + dy * dy);
        }
        prevX = x;
        prevY = y;
    }

    const double w = std::max(maxX - minX, Eps);
    const double h = std::max(maxY - minY, Eps);
    const double diag = std::sqrt(w * w + h * h);
    features.path_over_diag = diag <= Eps ? 0.0 : pathLen / diag;
    features.aspect = w / h;

    SampledPen sampled = BuildSampledPen(points, samples, cfg, DetermineFeatureNeeds(cfg, true));
    if (!sampled.curvature.empty())
    {
        double sum = 0.0;
        for (const double k : sampled.curvature)
            sum += std::abs(k);
        features.curvature_mean = sum / static_cast<double>(sampled.curvature.size());
    }
    return features;
}
} // namespace

SignChecker::SignChecker()
{
    DTWwindow = 1;
    SimpleChecks = 4;
    DTWChecks = static_cast<int>(config_.dtw_channels.size());
    ShapeChecks = 3;

    SimpleTestResult = 0.0;
    DTWTestResult = 0.0;
    ShapeTestResult = 0.0;

    mySimpleCheckFunctions = std::unique_ptr<SimpleCheckFunctions[]>(new SimpleCheckFunctions[SimpleChecks]);
    mySimpleCheckFunctions[0] = &RatioXY;
    mySimpleCheckFunctions[1] = &AvgXY;
    mySimpleCheckFunctions[2] = &SinXY;
    mySimpleCheckFunctions[3] = &RatioTouches;

    SimpleM = std::unique_ptr<Matrix<double>[]>(new Matrix<double>[SimpleChecks]);
    SimpleMch = std::unique_ptr<Matrix<double>[]>(new Matrix<double>[SimpleChecks]);

    DTWch = std::unique_ptr<Matrix<double>[]>(new Matrix<double>[DTWChecks]);
    DTW   = std::unique_ptr<Matrix<double>[]>(new Matrix<double>[DTWChecks]);
    CGch  = std::unique_ptr<Matrix<double>[]>(new Matrix<double>[DTWChecks]);
    CG    = std::unique_ptr<Matrix<double>[]>(new Matrix<double>[DTWChecks]);
    CVch  = std::unique_ptr<Matrix<double>[]>(new Matrix<double>[DTWChecks]);
    Ncg   = std::unique_ptr<Matrix<int>[]>(new Matrix<int>[DTWChecks]);

    ShapeM   = std::unique_ptr<Matrix<double>[]>(new Matrix<double>[ShapeChecks]);
    ShapeMch = std::unique_ptr<Matrix<double>[]>(new Matrix<double>[ShapeChecks]);

    ChannelWeights.assign(static_cast<std::size_t>(DTWChecks),
                          1.0 / static_cast<double>(std::max(DTWChecks, 1)));
}

void SignChecker::SetConfig(const SignCheckerConfig& config)
{
    config_ = config;
    if (config_.dtw_channels.empty())
        config_.dtw_channels = {SignCheckerConfig::DtwChannel::CityBlock};

    DTWwindow = std::max(1, config_.dtw_window);
    DTWChecks = static_cast<int>(config_.dtw_channels.size());

    DTWch = std::unique_ptr<Matrix<double>[]>(new Matrix<double>[DTWChecks]);
    DTW   = std::unique_ptr<Matrix<double>[]>(new Matrix<double>[DTWChecks]);
    CGch  = std::unique_ptr<Matrix<double>[]>(new Matrix<double>[DTWChecks]);
    CG    = std::unique_ptr<Matrix<double>[]>(new Matrix<double>[DTWChecks]);
    CVch  = std::unique_ptr<Matrix<double>[]>(new Matrix<double>[DTWChecks]);
    Ncg   = std::unique_ptr<Matrix<int>[]>(new Matrix<int>[DTWChecks]);

    ChannelWeights.assign(static_cast<std::size_t>(DTWChecks),
                          1.0 / static_cast<double>(DTWChecks));
}

void SignChecker::GenerateMatixeByIPen(int icheck, CTab *ipens, int npens)
{
    std::vector<double> values(static_cast<std::size_t>(npens));
    for (int i = 0; i < npens; i++)
        values[static_cast<std::size_t>(i)] = (*mySimpleCheckFunctions[icheck])(ipens, i);

    for (int i = 0; i < npens; i++)
    {
        for (int j = 0; j < i; j++)
        {
            double argval = std::abs(values[static_cast<std::size_t>(i)] - values[static_cast<std::size_t>(j)]);
            if (npens < 3)
            {
                const double denominator = std::max(
                    std::max(std::abs(values[static_cast<std::size_t>(i)]), std::abs(values[static_cast<std::size_t>(j)])),
                    Eps);
                argval /= denominator;
            }
            this->SimpleMch[icheck].SetElem(i, j, argval);
            this->SimpleMch[icheck].SetElem(j, i, argval);
        }
    }
}

bool SignChecker::SimpleCheckForRandomForge(CTab *ipens, int npens)
{
    this->SimpleTestResult = 0.0;
    if (ipens == nullptr || npens < 2)
        return false;

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
    this->DTWTestResult = 0.0;
    if (dpens == nullptr || spens == nullptr || npens < 2)
        return false;

    for (int i = 0; i < npens; ++i)
        if (spens[i].GetN() == 0)
            return false;

    const FeatureNeeds needs = DetermineFeatureNeeds(config_);
    std::vector<SampledPen> preparedPens;
    preparedPens.reserve(static_cast<std::size_t>(npens));
    for (int i = 0; i < npens; ++i)
        preparedPens.push_back(BuildSampledPen(dpens[i], spens[i], config_, needs));

    std::vector<double> rawScores(static_cast<std::size_t>(DTWChecks), 0.0);
    for (int i = 0; i < this->DTWChecks; i++)
    {
        this->DTWch[i].Init(npens, npens);
        this->CGch[i].Init(npens, npens);
        this->CVch[i].Init(npens, npens);
        this->Ncg[i].Init(npens, npens);
        const auto channel = config_.dtw_channels[static_cast<std::size_t>(i)];
        for (int ipen = 0; ipen < npens; ipen++)
        {
            for (int jpen = 0; jpen < ipen; jpen++)
            {
                StoreDTWPair(preparedPens[static_cast<std::size_t>(ipen)],
                             preparedPens[static_cast<std::size_t>(jpen)],
                             this->DTWch[i],
                             this->CGch[i],
                             this->CVch[i],
                             this->Ncg[i],
                             channel,
                             std::max(1, this->DTWwindow),
                             config_.sakoe_chiba_band,
                             ipen,
                             jpen);
            }
        }
        this->DTW[i].Init(this->DTWch[i], 1, 1);
        this->CG[i].Init(this->CGch[i], 1, 1);

        rawScores[static_cast<std::size_t>(i)] = 1.0 - this->DTWCalcDifference(i, npens);
    }

    if (config_.auto_weight_channels && DTWChecks > 1 && npens >= 3)
    {
        std::vector<double> weights(static_cast<std::size_t>(DTWChecks), 0.0);
        double weightSum = 0.0;
        for (int c = 0; c < DTWChecks; ++c)
        {
            double mean = 0.0;
            int pairs = 0;
            for (int i = 1; i < npens - 1; ++i)
                for (int j = 0; j < i; ++j)
                {
                    mean += this->DTWch[c](i, j);
                    ++pairs;
                }
            if (pairs == 0)
            {
                weights[static_cast<std::size_t>(c)] = 1.0;
                weightSum += 1.0;
                continue;
            }
            mean /= static_cast<double>(pairs);
            double var = 0.0;
            for (int i = 1; i < npens - 1; ++i)
                for (int j = 0; j < i; ++j)
                {
                    const double d = this->DTWch[c](i, j) - mean;
                    var += d * d;
                }
            var /= static_cast<double>(pairs);
            const double w = 1.0 / std::max(var, Eps);
            weights[static_cast<std::size_t>(c)] = w;
            weightSum += w;
        }
        if (weightSum > Eps)
            for (auto& w : weights)
                w /= weightSum;
        ChannelWeights = weights;
    }
    else
    {
        ChannelWeights.assign(static_cast<std::size_t>(DTWChecks),
                              1.0 / static_cast<double>(DTWChecks));
    }

    for (int i = 0; i < DTWChecks; ++i)
        this->DTWTestResult += ChannelWeights[static_cast<std::size_t>(i)] * rawScores[static_cast<std::size_t>(i)];

    return true;
}

bool SignChecker::ShapeCheckFromDPoints(DPoints* dpens, SPoints* spens, int npens)
{
    this->ShapeTestResult = 0.0;
    if (dpens == nullptr || spens == nullptr || npens < 2)
        return false;

    std::vector<ShapeFeatures> features;
    features.reserve(static_cast<std::size_t>(npens));
    for (int i = 0; i < npens; ++i)
        features.push_back(ComputeShape(dpens[i], spens[i], config_));

    auto fillMatrix = [&](int icheck, double ShapeFeatures::*field)
    {
        this->ShapeMch[icheck].Init(npens, npens);
        for (int i = 0; i < npens; ++i)
            for (int j = 0; j < i; ++j)
            {
                const double diff = std::abs(features[static_cast<std::size_t>(i)].*field
                                             - features[static_cast<std::size_t>(j)].*field);
                this->ShapeMch[icheck].SetElem(i, j, diff);
                this->ShapeMch[icheck].SetElem(j, i, diff);
            }
        this->ShapeM[icheck].Init(this->ShapeMch[icheck], 1, 1);
    };

    fillMatrix(0, &ShapeFeatures::path_over_diag);
    fillMatrix(1, &ShapeFeatures::aspect);
    fillMatrix(2, &ShapeFeatures::curvature_mean);

    for (int i = 0; i < ShapeChecks; ++i)
        this->ShapeTestResult += 1.0 - CalcDifference(this->ShapeMch[i], this->ShapeM[i], npens);
    this->ShapeTestResult /= static_cast<double>(ShapeChecks);
    return true;
}

double SignChecker::CalcDifference(Matrix<double>& ratioM_withCheckedSign,
                                   Matrix<double>& ratioM_withoutCheckedSign,
                                   int npens)
{
    if (npens < 2)
        return 1.0;

    const double checkedRms = config_.trim_outliers
        ? RowMedianPrefix(ratioM_withCheckedSign, npens - 1, npens - 1)
        : RowRmsPrefix(ratioM_withCheckedSign, npens - 1, npens - 1);

    if (npens < 3)
    {
        if (checkedRms > 1)
            return 1.0 - (1.0 / checkedRms);
        return checkedRms;
    }

    const double referenceRms = config_.trim_outliers
        ? PairwiseMedian(ratioM_withoutCheckedSign, ratioM_withoutCheckedSign.GetN())
        : PairwiseRms(ratioM_withoutCheckedSign, ratioM_withoutCheckedSign.GetN());
    return RelativeDifference(checkedRms, referenceRms);
}

double SignChecker::DTWCalcDifference(int icheck, int npens)
{
    if (npens <= 0)
        return 1.0;

    int jCVcheckmax = this->CVch[icheck].GetJMaxElInI(npens - 1);
    jCVcheckmax = std::min(jCVcheckmax, this->DTW[icheck].GetN() - 1);

    double DTWavg = this->DTW[icheck].AvgByI(jCVcheckmax);
    double sigma = std::sqrt(this->DTW[icheck].SumSqByI(jCVcheckmax, DTWavg));
    double sigmacheck = std::abs(this->DTWch[icheck](npens - 1, jCVcheckmax) - DTWavg);

    auto fold = [&](double s, double sc)
    {
        if (config_.zscore_scoring)
        {
            const double sc2 = sc * sc;
            const double s2 = s * s;
            const double denom = sc2 + s2;
            return denom <= Eps ? 0.0 : sc2 / denom;
        }
        return RelativeDifference(s, sc);
    };

    double ratval = fold(sigma, sigmacheck);

    double CGavg = this->CG[icheck].AvgByI(jCVcheckmax);
    sigma = std::sqrt(this->CG[icheck].SumSqByI(jCVcheckmax, CGavg));
    sigmacheck = std::abs(this->CGch[icheck](npens - 1, jCVcheckmax) - CGavg);
    ratval += fold(sigma, sigmacheck);

    return ratval / 2.0;
}

void SignChecker::GenerateMatrixByDPen(int icheck, DPoints* dpens, SPoints *spens, int npens)
{
    if (icheck < 0 || icheck >= static_cast<int>(config_.dtw_channels.size()) || dpens == nullptr || spens == nullptr)
        return;

    const FeatureNeeds needs = DetermineFeatureNeeds(config_);
    std::vector<SampledPen> preparedPens;
    preparedPens.reserve(static_cast<std::size_t>(npens));
    for (int i = 0; i < npens; ++i)
        preparedPens.push_back(BuildSampledPen(dpens[i], spens[i], config_, needs));

    const auto channel = config_.dtw_channels[static_cast<std::size_t>(icheck)];
    for (int i = 0; i < npens; i++)
        for (int j = 0; j < i; j++)
            StoreDTWPair(preparedPens[static_cast<std::size_t>(i)],
                         preparedPens[static_cast<std::size_t>(j)],
                         this->DTWch[icheck],
                         this->CGch[icheck],
                         this->CVch[icheck],
                         this->Ncg[icheck],
                         channel,
                         std::max(1, this->DTWwindow),
                         config_.sakoe_chiba_band,
                         i,
                         j);
}

void SignChecker::DTW_Go(int icheck, DPoints* dpens, SPoints *spens, int i, int j)
{
    if (spens[i].GetN() == 0 || spens[j].GetN() == 0)
        return;

    if (icheck < 0 || icheck >= static_cast<int>(config_.dtw_channels.size()))
        return;

    const FeatureNeeds needs = DetermineFeatureNeeds(config_);
    const SampledPen left = BuildSampledPen(dpens[i], spens[i], config_, needs);
    const SampledPen right = BuildSampledPen(dpens[j], spens[j], config_, needs);
    StoreDTWPair(left,
                 right,
                 this->DTWch[icheck],
                 this->CGch[icheck],
                 this->CVch[icheck],
                 this->Ncg[icheck],
                 config_.dtw_channels[static_cast<std::size_t>(icheck)],
                 std::max(1, this->DTWwindow),
                 config_.sakoe_chiba_band,
                 i,
                 j);
}

bool SignChecker::DTW_InitMatrix(int icheck, Matrix<double>& D, DPoints& dpeni, SPoints& speni, DPoints& dpenj, SPoints& spenj)
{
    const FeatureNeeds needs = DetermineFeatureNeeds(config_);
    const SampledPen left = BuildSampledPen(dpeni, speni, config_, needs);
    const SampledPen right = BuildSampledPen(dpenj, spenj, config_, needs);
    const int Ni = static_cast<int>(left.x.size());
    const int Nj = static_cast<int>(right.x.size());
    if (Ni == 0 || Nj == 0 || D.GetN() != Ni || D.GetM() != Nj)
        return false;

    using DC = SignCheckerConfig::DtwChannel;
    if (icheck < 0 || icheck >= static_cast<int>(config_.dtw_channels.size()))
        return false;
    const DC channel = config_.dtw_channels[static_cast<std::size_t>(icheck)];
    return FillDTWMatrix(D, channel, left, right, std::max(1, this->DTWwindow), config_.sakoe_chiba_band);
}

void SignChecker::DTW_TraceBack(Matrix<double>& D, SPoints& tranformpeni, SPoints& tranformpenj)
{
    TraceBackMatrix(D, std::max(1, this->DTWwindow), tranformpeni, tranformpenj);
}

double SignChecker::DTW_CaclGlobalDeformation(Matrix<double>& D, SPoints& tranformpeni, SPoints& tranformpenj)
{
    return GlobalDeformation(D, tranformpeni, tranformpenj);
}

double SignChecker::DTW_CalcDetermination(SPoints& tpeni, DPoints& dpeni, SPoints& tpenj, DPoints& dpenj)
{
    int N = tpeni.GetN();
    if (N == 0 || tpenj.GetN() != N)
        return 0.0;

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
    const double determination = (sumi + sumj) * (sumi + sumj) / (sumsqx * sumsqy);
    return std::max(0.0, std::min(1.0, determination));
}
