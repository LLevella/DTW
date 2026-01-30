#include "SignChecker.h"

#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace
{
DPoints MakeArc(int n, double phaseShift, double scale)
{
    DPoints points;
    for (int i = 0; i < n; ++i)
    {
        const double t = static_cast<double>(i) / static_cast<double>(n - 1);
        const double angle = phaseShift + t * 3.14159265358979323846;
        const double pressure = 0.5 + 0.1 * std::sin(2.0 * angle);
        const double timestamp = static_cast<double>(i) * 0.01;
        points.AddPoint(scale * std::cos(angle), scale * std::sin(angle), pressure, timestamp, true);
    }
    return points;
}

std::vector<DPoints> MakeDataset(int signatures, int points)
{
    std::vector<DPoints> result;
    result.reserve(static_cast<std::size_t>(signatures));
    for (int i = 0; i < signatures; ++i)
    {
        const double phase = 0.002 * static_cast<double>(i % 5);
        const double scale = 1.0 + 0.01 * static_cast<double>(i % 7);
        result.push_back(MakeArc(points, phase, scale));
    }
    return result;
}

std::vector<SPoints> FullRanges(const std::vector<DPoints>& points)
{
    std::vector<SPoints> result;
    result.reserve(points.size());
    for (const auto& pen : points)
        result.push_back(SPoints::FullRange(pen.GetN()));
    return result;
}

void RunCase(const std::string& name, SignCheckerConfig config, int signatures, int points, int iterations)
{
    auto pens = MakeDataset(signatures, points);
    auto samples = FullRanges(pens);
    SignChecker checker;
    checker.SetConfig(config);
    balance = 0.0;

    const auto started = std::chrono::steady_clock::now();
    double score = 0.0;
    for (int i = 0; i < iterations; ++i)
    {
        if (!checker.DTWCheckForSimpleForge(pens.data(), samples.data(), static_cast<int>(pens.size())))
        {
            std::cerr << name << ": DTW check failed\n";
            return;
        }
        score += checker.GetDTWResult();
    }
    const auto finished = std::chrono::steady_clock::now();
    const auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finished - started).count();

    std::cout << name
              << " signatures=" << signatures
              << " points=" << points
              << " iterations=" << iterations
              << " total_us=" << elapsed
              << " avg_score=" << (score / static_cast<double>(iterations))
              << '\n';
}
}

int main()
{
    SignCheckerConfig baseline;
    RunCase("baseline", baseline, 24, 96, 10);

    SignCheckerConfig full;
    full.arclength_resample = true;
    full.resample_points = 96;
    full.smooth = true;
    full.smooth_window = 5;
    full.sakoe_chiba_band = 12;
    full.auto_weight_channels = true;
    full.dtw_channels = {
        SignCheckerConfig::DtwChannel::CityBlock,
        SignCheckerConfig::DtwChannel::TangentAngle,
        SignCheckerConfig::DtwChannel::Curvature,
        SignCheckerConfig::DtwChannel::Pressure,
        SignCheckerConfig::DtwChannel::Velocity,
    };
    RunCase("full", full, 24, 96, 10);
}
