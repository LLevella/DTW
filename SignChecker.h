#pragma once

#include <memory>
#include <vector>

#include "Matrix.h"
#include "SignatureData.h"

extern double balance;

struct SignCheckerConfig
{
    enum class DtwChannel : int
    {
        CityBlock      = 0, // |dx| + |dy| между нормализованными координатами
        Direction      = 1, // |cross| / |v_i||v_j| (текущая sin-метрика направлений)
        SegmentLength  = 2, // |seglen_i - seglen_j|
        TangentAngle   = 3, // угол между касательными по модулю π
        Curvature      = 4, // |κ_i - κ_j|
        Pseudopressure = 5, // 1 / seglen — оценка скорости через плотность точек
        Pressure       = 6, // |p_i - p_j| (требует has_pressure)
        Velocity       = 7, // |v_i - v_j| (требует has_timestamp)
    };

    // --- Препроцессинг ---
    bool smooth = false;             // скользящее среднее по x,y,p
    int  smooth_window = 3;          // нечётный размер окна
    bool arclength_resample = false; // равномерный ресемпл по длине дуги
    int  resample_points = 64;
    bool pca_rotate = false;         // выровнять по главной оси
    bool zscore_normalize = false;   // z-score per axis вместо bbox-diag scale

    // --- DTW шаблон шага и глобальные ограничения ---
    int dtw_window = 1;              // радиус локального шага
    int sakoe_chiba_band = -1;       // -1 = выключено; иначе радиус полосы вокруг диагонали

    // --- Каналы DTW ---
    std::vector<DtwChannel> dtw_channels = {
        DtwChannel::CityBlock,
        DtwChannel::Direction,
        DtwChannel::SegmentLength,
    };

    // --- Решающее правило ---
    bool auto_weight_channels = false; // взвешивать каналы по 1/Var(intra-ref)
    bool trim_outliers = false;        // медиана/MAD вместо среднего/RMS
    bool zscore_scoring = false;       // sigmacheck^2 / (sigmacheck^2 + sigma^2)
    double shape_weight = 0.0;         // вклад shape-проверок в GetTestResult, 0..1
};

class SignChecker
{
    int SimpleChecks;
    int DTWChecks;
    int ShapeChecks;

    int DTWwindow;

    double DTWTestResult;
    double SimpleTestResult;
    double ShapeTestResult;

    typedef double(*SimpleCheckFunctions)(CTab *ipens, int i);
    std::unique_ptr<SimpleCheckFunctions[]> mySimpleCheckFunctions;

    std::unique_ptr<Matrix<double>[]> DTWch;
    std::unique_ptr<Matrix<double>[]> DTW;

    std::unique_ptr<Matrix<double>[]> CGch;
    std::unique_ptr<Matrix<double>[]> CG;

    std::unique_ptr<Matrix<double>[]> CVch;
    std::unique_ptr<Matrix<int>[]>    Ncg;

    std::unique_ptr<Matrix<double>[]> SimpleM;
    std::unique_ptr<Matrix<double>[]> SimpleMch;

    std::unique_ptr<Matrix<double>[]> ShapeM;
    std::unique_ptr<Matrix<double>[]> ShapeMch;

    std::vector<double> ChannelWeights;

    SignCheckerConfig config_;

public:
    SignChecker();

    void SetConfig(const SignCheckerConfig& config);
    const SignCheckerConfig& GetConfig() const { return config_; }

    inline double GetTestResult()
    {
        const double w_shape = std::max(0.0, std::min(1.0, config_.shape_weight));
        const double w_rest  = 1.0 - w_shape;
        const double simple  = balance * SimpleTestResult + (1.0 - balance) * DTWTestResult;
        return w_rest * simple + w_shape * ShapeTestResult;
    }

    inline double GetSimpleResult() const { return SimpleTestResult; }
    inline double GetDTWResult() const    { return DTWTestResult; }
    inline double GetShapeResult() const  { return ShapeTestResult; }

    bool SimpleCheckForRandomForge(CTab *ipens, int npens);
    bool DTWCheckForSimpleForge(DPoints* dpens, SPoints* spens, int npen);
    bool ShapeCheckFromDPoints(DPoints* dpens, SPoints* spens, int npens);

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
