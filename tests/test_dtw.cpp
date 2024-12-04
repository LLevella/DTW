#include "Matrix.h"
#include "SignChecker.h"

#include <cassert>
#include <cmath>
#include <stdexcept>

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
}

int main()
{
    TestMatrixStorageAndExtremes();
    TestTraceBackHandlesEqualCosts();
    TestSimpleChecksResetAndAvoidZeroDivision();
    TestDTWIdenticalSignatures();
    TestDTWNormalizesTranslationAndScale();
}
