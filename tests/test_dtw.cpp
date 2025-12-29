#include "Matrix.h"
#include "SignChecker.h"

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace
{
void ExpectNear(double actual, double expected)
{
    assert(std::abs(actual - expected) < 1e-9);
}

DPoints MakeLine()
{
    return DPoints{{0.0, 0.0}, {1.0, 1.0}, {2.0, 0.0}};
}

DPoints MakeShiftedScaledLine()
{
    return DPoints{{10.0, -3.0}, {12.0, -1.0}, {14.0, -3.0}};
}

DPoints MakePressureTimedLine(double pressureOffset = 0.0, double timeScale = 1.0)
{
    DPoints points;
    points.AddPoint(0.0, 0.0, 0.4 + pressureOffset, 0.00 * timeScale, true);
    points.AddPoint(1.0, 0.0, 0.6 + pressureOffset, 0.10 * timeScale, true);
    points.AddPoint(2.0, 0.0, 0.5 + pressureOffset, 0.20 * timeScale, false);
    return points;
}

DPoints MakeArc(int n, double phaseShift = 0.0, double scale = 1.0)
{
    DPoints points;
    for (int i = 0; i < n; ++i)
    {
        const double t = static_cast<double>(i) / static_cast<double>(n - 1);
        const double angle = phaseShift + t * 3.14159265358979;
        points.AddPoint(scale * std::cos(angle), scale * std::sin(angle));
    }
    return points;
}

void TestMatrixStorageAndExtremes()
{
    Matrix<double> matrix;
    matrix.Init(2, 3);
    matrix.SetElem(0, 2, 5.0);
    matrix.SetElem(1, 0, 2.0);

    int i = -1;
    int j = -1;
    ExpectNear(matrix.GetIndMaxElem(i, j), 5.0);
    assert(i == 0);
    assert(j == 2);

    ExpectNear(matrix.GetIndMinElem(i, j), 2.0);
    assert(i == 1);
    assert(j == 0);
    assert(matrix.GetJMinElInI(0) == 2);

    Matrix<double> copy = matrix;
    copy.SetElem(0, 2, 9.0);
    ExpectNear(matrix(0, 2), 5.0);
    ExpectNear(copy(0, 2), 9.0);

    bool threw = false;
    try
    {
        Matrix<double> empty;
        empty.AvgAll();
    }
    catch (const std::out_of_range&)
    {
        threw = true;
    }
    assert(threw);
}

void TestTraceBackHandlesEqualCosts()
{
    Matrix<double> costs(3, 3);
    SignChecker checker;
    SPoints transformedI;
    SPoints transformedJ;

    checker.DTW_TraceBack(costs, transformedI, transformedJ);

    assert(transformedI.GetN() == 3);
    assert(transformedJ.GetN() == 3);
    assert(transformedI[0] == 2);
    assert(transformedJ[0] == 2);
    assert(transformedI[1] == 1);
    assert(transformedJ[1] == 1);
    assert(transformedI[2] == 0);
    assert(transformedJ[2] == 0);
}

void TestSimpleChecksResetAndAvoidZeroDivision()
{
    balance = 1.0;
    CTab pens[2] = {};
    SignChecker checker;

    assert(checker.SimpleCheckForRandomForge(pens, 2));
    ExpectNear(checker.GetTestResult(), 1.0);

    assert(checker.SimpleCheckForRandomForge(pens, 2));
    ExpectNear(checker.GetTestResult(), 1.0);

    assert(!checker.SimpleCheckForRandomForge(nullptr, 2));
}

void TestSimpleDifferenceUsesReferenceRms()
{
    SignChecker checker;
    Matrix<double> withChecked(3, 3);
    Matrix<double> withoutChecked(2, 2);

    withChecked.SetElem(0, 1, 2.0);
    withChecked.SetElem(1, 0, 2.0);
    withChecked.SetElem(2, 0, 2.0);
    withChecked.SetElem(0, 2, 2.0);
    withChecked.SetElem(2, 1, 2.0);
    withChecked.SetElem(1, 2, 2.0);
    withoutChecked.SetElem(0, 1, 2.0);
    withoutChecked.SetElem(1, 0, 2.0);

    ExpectNear(checker.CalcDifference(withChecked, withoutChecked, 3), 0.0);

    withChecked.SetElem(2, 0, 8.0);
    withChecked.SetElem(0, 2, 8.0);
    withChecked.SetElem(2, 1, 8.0);
    withChecked.SetElem(1, 2, 8.0);
    ExpectNear(checker.CalcDifference(withChecked, withoutChecked, 3), 0.75);
}

void TestDTWIdenticalSignatures()
{
    balance = 0.0;
    DPoints dpens[2] = {MakeLine(), MakeLine()};
    SPoints spens[2] = {SPoints::FullRange(dpens[0].GetN()), SPoints::FullRange(dpens[1].GetN())};
    SignChecker checker;

    assert(checker.DTWCheckForSimpleForge(dpens, spens, 2));
    ExpectNear(checker.GetTestResult(), 1.0);

    assert(checker.DTWCheckForSimpleForge(dpens, spens, 2));
    ExpectNear(checker.GetTestResult(), 1.0);

    Matrix<double> costs(3, 3);
    assert(checker.DTW_InitMatrix(0, costs, dpens[0], spens[0], dpens[1], spens[1]));
    ExpectNear(costs(2, 2), 0.0);
}

void TestDTWNormalizesTranslationAndScale()
{
    DPoints dpens[2] = {MakeLine(), MakeShiftedScaledLine()};
    SPoints spens[2] = {SPoints::FullRange(dpens[0].GetN()), SPoints::FullRange(dpens[1].GetN())};
    SignChecker checker;

    Matrix<double> cityBlockCosts(3, 3);
    assert(checker.DTW_InitMatrix(0, cityBlockCosts, dpens[0], spens[0], dpens[1], spens[1]));
    ExpectNear(cityBlockCosts(2, 2), 0.0);

    Matrix<double> velocityCosts(3, 3);
    assert(checker.DTW_InitMatrix(2, velocityCosts, dpens[0], spens[0], dpens[1], spens[1]));
    ExpectNear(velocityCosts(2, 2), 0.0);
}

void TestExtendedDPointFieldsPreserved()
{
    DPoints points;
    points.AddPoint(0.0, 0.0, 0.5, 0.0, true);
    points.AddPoint(1.0, 0.0, 0.7, 0.01, true);
    points.AddPoint(2.0, 0.0, 0.6, 0.02, false);

    assert(points.AllHavePressure());
    assert(points.AllHaveTimestamp());
    assert(points.GetPressure(1) == 0.7);
    assert(points.GetTimestamp(2) == 0.02);
    assert(points.GetPenDown(2) == false);

    DPoints implicit{{0.0, 0.0}, {1.0, 1.0}};
    assert(!implicit.AllHavePressure());
    assert(implicit.GetPenDown(0) == true);
}

void TestCustomChannelsViaConfig()
{
    DPoints arc1 = MakeArc(20);
    DPoints arc2 = MakeArc(20, 0.0, 2.0);
    DPoints dpens[2] = {arc1, arc2};
    SPoints spens[2] = {SPoints::FullRange(arc1.GetN()), SPoints::FullRange(arc2.GetN())};

    SignChecker checker;
    SignCheckerConfig cfg;
    cfg.dtw_channels = {
        SignCheckerConfig::DtwChannel::TangentAngle,
        SignCheckerConfig::DtwChannel::Curvature,
    };
    checker.SetConfig(cfg);

    balance = 0.0;
    assert(checker.DTWCheckForSimpleForge(dpens, spens, 2));
    assert(checker.GetDTWResult() > 0.99); // одинаковая кривизна и направления — должны совпадать
}

void TestPressureAndVelocityChannels()
{
    DPoints dpens[2] = {MakePressureTimedLine(), MakePressureTimedLine()};
    SPoints spens[2] = {SPoints::FullRange(dpens[0].GetN()), SPoints::FullRange(dpens[1].GetN())};
    SignChecker checker;
    SignCheckerConfig cfg;
    cfg.dtw_channels = {
        SignCheckerConfig::DtwChannel::Pressure,
        SignCheckerConfig::DtwChannel::Velocity,
        SignCheckerConfig::DtwChannel::Pseudopressure,
    };
    checker.SetConfig(cfg);

    Matrix<double> pressureCosts(3, 3);
    assert(checker.DTW_InitMatrix(0, pressureCosts, dpens[0], spens[0], dpens[1], spens[1]));
    ExpectNear(pressureCosts(2, 2), 0.0);

    Matrix<double> velocityCosts(3, 3);
    assert(checker.DTW_InitMatrix(1, velocityCosts, dpens[0], spens[0], dpens[1], spens[1]));
    ExpectNear(velocityCosts(2, 2), 0.0);

    Matrix<double> pseudoPressureCosts(3, 3);
    assert(checker.DTW_InitMatrix(2, pseudoPressureCosts, dpens[0], spens[0], dpens[1], spens[1]));
    ExpectNear(pseudoPressureCosts(2, 2), 0.0);
}

void TestEmptyDtwChannelsFallsBackToCityBlock()
{
    DPoints dpens[2] = {MakeLine(), MakeLine()};
    SPoints spens[2] = {SPoints::FullRange(dpens[0].GetN()), SPoints::FullRange(dpens[1].GetN())};
    SignChecker checker;
    SignCheckerConfig cfg;
    cfg.dtw_channels.clear();
    checker.SetConfig(cfg);

    balance = 0.0;
    assert(checker.GetConfig().dtw_channels.size() == 1);
    assert(checker.DTWCheckForSimpleForge(dpens, spens, 2));
    ExpectNear(checker.GetDTWResult(), 1.0);
}

void TestSakoeChibaBandRejectsFarPaths()
{
    DPoints dpens[2] = {MakeArc(16), MakeArc(16)};
    SPoints spens[2] = {SPoints::FullRange(16), SPoints::FullRange(16)};

    SignChecker checker;
    SignCheckerConfig cfg;
    cfg.sakoe_chiba_band = 0; // только диагональ
    checker.SetConfig(cfg);

    Matrix<double> costs(16, 16);
    assert(checker.DTW_InitMatrix(0, costs, dpens[0], spens[0], dpens[1], spens[1]));
    assert(!std::isfinite(costs(0, 15))); // дальний угол вне полосы
    assert(std::isfinite(costs(8, 8)));
}

void TestArclengthResampleStableOnIdenticalSignatures()
{
    DPoints dpens[2] = {MakeArc(15), MakeArc(15)};
    SPoints spens[2] = {SPoints::FullRange(15), SPoints::FullRange(15)};

    SignChecker checker;
    SignCheckerConfig cfg;
    cfg.arclength_resample = true;
    cfg.resample_points = 32;
    checker.SetConfig(cfg);

    balance = 0.0;
    assert(checker.DTWCheckForSimpleForge(dpens, spens, 2));
    ExpectNear(checker.GetDTWResult(), 1.0);
}

void TestShapeCheckIdenticalGivesPerfectScore()
{
    DPoints dpens[2] = {MakeArc(20), MakeArc(20)};
    SPoints spens[2] = {SPoints::FullRange(20), SPoints::FullRange(20)};

    SignChecker checker;
    assert(checker.ShapeCheckFromDPoints(dpens, spens, 2));
    ExpectNear(checker.GetShapeResult(), 1.0);
}

void TestShapeWeightContributesToTotal()
{
    DPoints dpens[2] = {MakeArc(20), MakeArc(20)};
    SPoints spens[2] = {SPoints::FullRange(20), SPoints::FullRange(20)};

    SignChecker checker;
    SignCheckerConfig cfg;
    cfg.shape_weight = 0.5;
    checker.SetConfig(cfg);

    balance = 0.0;
    assert(checker.DTWCheckForSimpleForge(dpens, spens, 2));
    assert(checker.ShapeCheckFromDPoints(dpens, spens, 2));
    ExpectNear(checker.GetTestResult(), 1.0);
}

void TestZscoreScoringSimpleCase()
{
    DPoints dpens[2] = {MakeArc(10), MakeArc(10)};
    SPoints spens[2] = {SPoints::FullRange(10), SPoints::FullRange(10)};

    SignChecker checker;
    SignCheckerConfig cfg;
    cfg.zscore_scoring = true;
    checker.SetConfig(cfg);

    balance = 0.0;
    assert(checker.DTWCheckForSimpleForge(dpens, spens, 2));
    ExpectNear(checker.GetDTWResult(), 1.0);
}

void TestSmoothingPcaAndZscoreKeepIdenticalSignatures()
{
    DPoints dpens[2] = {MakeArc(24), MakeArc(24)};
    SPoints spens[2] = {SPoints::FullRange(24), SPoints::FullRange(24)};

    SignChecker checker;
    SignCheckerConfig cfg;
    cfg.smooth = true;
    cfg.smooth_window = 5;
    cfg.pca_rotate = true;
    cfg.zscore_normalize = true;
    cfg.dtw_channels = {
        SignCheckerConfig::DtwChannel::CityBlock,
        SignCheckerConfig::DtwChannel::TangentAngle,
        SignCheckerConfig::DtwChannel::Curvature,
    };
    checker.SetConfig(cfg);

    balance = 0.0;
    assert(checker.DTWCheckForSimpleForge(dpens, spens, 2));
    ExpectNear(checker.GetDTWResult(), 1.0);
}

void TestTrimOutliersUsesMedianDistance()
{
    SignChecker checker;
    SignCheckerConfig cfg;
    cfg.trim_outliers = true;
    checker.SetConfig(cfg);

    Matrix<double> withChecked(4, 4);
    Matrix<double> withoutChecked(3, 3);

    withoutChecked.SetElem(1, 0, 2.0);
    withoutChecked.SetElem(0, 1, 2.0);
    withoutChecked.SetElem(2, 0, 2.0);
    withoutChecked.SetElem(0, 2, 2.0);
    withoutChecked.SetElem(2, 1, 100.0);
    withoutChecked.SetElem(1, 2, 100.0);

    withChecked.SetElem(3, 0, 2.0);
    withChecked.SetElem(0, 3, 2.0);
    withChecked.SetElem(3, 1, 2.0);
    withChecked.SetElem(1, 3, 2.0);
    withChecked.SetElem(3, 2, 2.0);
    withChecked.SetElem(2, 3, 2.0);

    ExpectNear(checker.CalcDifference(withChecked, withoutChecked, 4), 0.0);
}

void TestAutoWeightChannelsNormalizes()
{
    DPoints dpens[3] = {MakeArc(12), MakeArc(12), MakeArc(12)};
    SPoints spens[3] = {SPoints::FullRange(12), SPoints::FullRange(12), SPoints::FullRange(12)};

    SignChecker checker;
    SignCheckerConfig cfg;
    cfg.auto_weight_channels = true;
    checker.SetConfig(cfg);

    balance = 0.0;
    assert(checker.DTWCheckForSimpleForge(dpens, spens, 3));
    // Все три подписи одинаковы, так что результат всё равно близок к 1
    assert(checker.GetDTWResult() > 0.95);
}
}

int main()
{
    TestMatrixStorageAndExtremes();
    TestTraceBackHandlesEqualCosts();
    TestSimpleChecksResetAndAvoidZeroDivision();
    TestSimpleDifferenceUsesReferenceRms();
    TestDTWIdenticalSignatures();
    TestDTWNormalizesTranslationAndScale();
    TestExtendedDPointFieldsPreserved();
    TestCustomChannelsViaConfig();
    TestPressureAndVelocityChannels();
    TestEmptyDtwChannelsFallsBackToCityBlock();
    TestSakoeChibaBandRejectsFarPaths();
    TestArclengthResampleStableOnIdenticalSignatures();
    TestShapeCheckIdenticalGivesPerfectScore();
    TestShapeWeightContributesToTotal();
    TestZscoreScoringSimpleCase();
    TestSmoothingPcaAndZscoreKeepIdenticalSignatures();
    TestTrimOutliersUsesMedianDistance();
    TestAutoWeightChannelsNormalizes();
}
