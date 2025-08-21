<!-- markdownlint-disable MD024 -->

# DTW

[English](#english) | [Русский](#русский)

---

## English

DTW is a small C++17 library for handwritten signature comparison. It combines simple signature features, configurable Dynamic Time Warping channels, shape-level checks, and optional preprocessing for noisy point streams.

### Features

- simple feature comparison via `CTab`;
- DTW over normalized signature trajectories;
- configurable DTW channels:
  - city-block coordinate distance;
  - direction;
  - segment length;
  - tangent angle;
  - curvature;
  - pseudo-pressure from point density;
  - real pressure when available;
  - velocity when timestamps are available;
- optional smoothing, arclength resampling, PCA rotation, and z-score normalization;
- optional Sakoe-Chiba band constraint for DTW;
- shape checks based on path length, aspect ratio, and curvature;
- optional channel auto-weighting, outlier trimming, and z-score scoring;
- CMake package install support for `find_package(DTW)`;
- CI matrix plus AddressSanitizer/UndefinedBehaviorSanitizer checks.

### Layout

```text
.
├── Matrix.h                 # matrix template
├── SignatureData.h/.cpp     # signature points, feature structs, legacy metrics
├── SignChecker.h/.cpp       # verification logic and SignCheckerConfig
├── tests/                   # regression tests
├── cmake/                   # CMake package config template
└── .github/workflows/ci.yml # CI, sanitizer job, release artifact packaging
```

### Requirements

- CMake 3.16 or newer;
- a C++17-capable compiler, for example `g++` or `clang++`.

### Quick Start

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Release build:

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --parallel
ctest --test-dir build-release --output-on-failure
```

Sanitizer build:

```bash
cmake -S . -B build-sanitize \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined" \
  -DCMAKE_SHARED_LINKER_FLAGS="-fsanitize=address,undefined"
cmake --build build-sanitize --parallel
ctest --test-dir build-sanitize --output-on-failure
```

### Installation

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --parallel
cmake --install build-release --prefix dist
```

Use the installed package from another CMake project:

```cmake
find_package(DTW REQUIRED)

add_executable(app main.cpp)
target_link_libraries(app PRIVATE DTW::dtw)
```

For a local install prefix:

```bash
cmake -S your-app -B your-app/build -DCMAKE_PREFIX_PATH=/path/to/DTW/dist
```

### Usage Example

```cpp
#include "SignChecker.h"

int main()
{
    DPoints reference{{0.0, 0.0}, {1.0, 1.0}, {2.0, 0.0}};
    DPoints checked{{10.0, -3.0}, {12.0, -1.0}, {14.0, -3.0}};

    DPoints signatures[2] = {reference, checked};
    SPoints samples[2] = {
        SPoints::FullRange(signatures[0].GetN()),
        SPoints::FullRange(signatures[1].GetN()),
    };

    SignCheckerConfig config;
    config.arclength_resample = true;
    config.resample_points = 64;
    config.sakoe_chiba_band = 8;
    config.dtw_channels = {
        SignCheckerConfig::DtwChannel::CityBlock,
        SignCheckerConfig::DtwChannel::TangentAngle,
        SignCheckerConfig::DtwChannel::Curvature,
    };
    config.shape_weight = 0.25;

    SignChecker checker;
    checker.SetConfig(config);

    balance = 0.0; // 0.0: DTW only, 1.0: simple checks only

    if (!checker.DTWCheckForSimpleForge(signatures, samples, 2))
        return 1;
    if (!checker.ShapeCheckFromDPoints(signatures, samples, 2))
        return 1;

    const double score = checker.GetTestResult();
    return score >= 0.9 ? 0 : 2;
}
```

`SignChecker` compares a set of signatures. The last signature in the array is treated as the checked signature; earlier signatures form the reference base.

### Pressure and Timestamp Data

`DPoints` can store only coordinates or extended point data:

```cpp
DPoints points;
points.AddPoint(0.0, 0.0, 0.45, 0.00, true);
points.AddPoint(1.0, 0.2, 0.60, 0.01, true);
points.AddPoint(2.0, 0.0, 0.50, 0.02, false);
```

Use `DtwChannel::Pressure` only when all compared points have pressure. Use `DtwChannel::Velocity` only when all compared points have timestamps. Missing data is handled safely and contributes zero cost for that channel.

### Configuration

`SignCheckerConfig` controls preprocessing, DTW channels, and scoring:

- `smooth`, `smooth_window` - moving average over coordinates and pressure;
- `arclength_resample`, `resample_points` - uniform sampling by path length;
- `pca_rotate` - align the signature by its principal axis;
- `zscore_normalize` - normalize coordinates per axis instead of bounding-box diagonal scale;
- `dtw_window` - local DTW step radius;
- `sakoe_chiba_band` - global band around the diagonal, `-1` disables it;
- `dtw_channels` - list of DTW channels to combine;
- `auto_weight_channels` - weight channels by inverse intra-reference variance;
- `trim_outliers` - use median-based scoring for simple/shape differences;
- `zscore_scoring` - alternative DTW scoring rule;
- `shape_weight` - contribution of shape checks to `GetTestResult()`.

### Core Types

- `CTab` stores simple features: `ratio_xy`, `avg_xy`, `sin_xy`, `ratio_touches`.
- `DPoint` stores coordinates plus optional pressure, timestamp, and pen-down state.
- `DPoints` stores a signature point sequence.
- `SPoints` stores selected point indices in `DPoints`.
- `SignCheckerConfig` describes preprocessing, DTW channels, and scoring.
- `SignChecker` runs simple, DTW, and shape checks.
- `Matrix<T>` is used for distance matrices and intermediate computations.

### Tests

The tests live in `tests/test_dtw.cpp` and cover:

- safe copying and indexing of `Matrix`;
- matrix extrema search;
- stable traceback under equal DTW costs;
- result reset across repeated checks;
- translation and scale normalization;
- pressure, velocity, and pseudo-pressure channels;
- tangent and curvature channels;
- empty DTW channel fallback;
- Sakoe-Chiba band behavior;
- arclength resampling;
- shape score contribution;
- z-score scoring;
- auto channel weighting.

Run:

```bash
ctest --test-dir build --output-on-failure
```

### CI/CD

GitHub Actions runs on pull requests, pushes to `master`/`main`, manual dispatch, and tags matching `v*`.

CI jobs:

- `g++` and `clang++` in `Debug` and `Release`;
- AddressSanitizer and UndefinedBehaviorSanitizer;
- release packaging on `v*` tags.

Release artifacts contain the static library, headers, and CMake package config.

### Notes

The public `SignChecker` methods from the original project are still available. New behavior is opt-in through `SignCheckerConfig` unless it is a safe default, such as normalized DTW coordinates.

---

## Русский

DTW - небольшая C++17-библиотека для сравнения рукописных подписей. Она сочетает простые признаки, настраиваемые DTW-каналы, shape-проверки и опциональный препроцессинг для шумных потоков точек.

### Возможности

- сравнение простых признаков через `CTab`;
- DTW по нормализованным траекториям подписи;
- настраиваемые DTW-каналы:
  - city-block расстояние по координатам;
  - направление;
  - длина сегмента;
  - угол касательной;
  - кривизна;
  - pseudo-pressure по плотности точек;
  - реальное давление, если оно есть;
  - скорость, если есть timestamps;
- smoothing, arclength resampling, PCA rotation и z-score normalization;
- ограничение пути DTW полосой Sakoe-Chiba;
- shape-проверки по длине пути, aspect ratio и кривизне;
- auto-weight каналов, trim-outliers и z-score scoring;
- установка через CMake package для `find_package(DTW)`;
- CI-матрица плюс проверки AddressSanitizer/UndefinedBehaviorSanitizer.

### Структура

```text
.
├── Matrix.h                 # шаблон матрицы
├── SignatureData.h/.cpp     # точки подписи, признаки, legacy-метрики
├── SignChecker.h/.cpp       # логика проверки и SignCheckerConfig
├── tests/                   # регрессионные тесты
├── cmake/                   # шаблон CMake package config
└── .github/workflows/ci.yml # CI, sanitizer job, упаковка релиза
```

### Требования

- CMake 3.16 или новее;
- компилятор с поддержкой C++17, например `g++` или `clang++`.

### Быстрый Старт

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Release-сборка:

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --parallel
ctest --test-dir build-release --output-on-failure
```

Sanitizer-сборка:

```bash
cmake -S . -B build-sanitize \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined" \
  -DCMAKE_SHARED_LINKER_FLAGS="-fsanitize=address,undefined"
cmake --build build-sanitize --parallel
ctest --test-dir build-sanitize --output-on-failure
```

### Установка

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --parallel
cmake --install build-release --prefix dist
```

Подключение из другого CMake-проекта:

```cmake
find_package(DTW REQUIRED)

add_executable(app main.cpp)
target_link_libraries(app PRIVATE DTW::dtw)
```

Если установка лежит локально:

```bash
cmake -S your-app -B your-app/build -DCMAKE_PREFIX_PATH=/path/to/DTW/dist
```

### Пример Использования

```cpp
#include "SignChecker.h"

int main()
{
    DPoints reference{{0.0, 0.0}, {1.0, 1.0}, {2.0, 0.0}};
    DPoints checked{{10.0, -3.0}, {12.0, -1.0}, {14.0, -3.0}};

    DPoints signatures[2] = {reference, checked};
    SPoints samples[2] = {
        SPoints::FullRange(signatures[0].GetN()),
        SPoints::FullRange(signatures[1].GetN()),
    };

    SignCheckerConfig config;
    config.arclength_resample = true;
    config.resample_points = 64;
    config.sakoe_chiba_band = 8;
    config.dtw_channels = {
        SignCheckerConfig::DtwChannel::CityBlock,
        SignCheckerConfig::DtwChannel::TangentAngle,
        SignCheckerConfig::DtwChannel::Curvature,
    };
    config.shape_weight = 0.25;

    SignChecker checker;
    checker.SetConfig(config);

    balance = 0.0; // 0.0: только DTW, 1.0: только простые проверки

    if (!checker.DTWCheckForSimpleForge(signatures, samples, 2))
        return 1;
    if (!checker.ShapeCheckFromDPoints(signatures, samples, 2))
        return 1;

    const double score = checker.GetTestResult();
    return score >= 0.9 ? 0 : 2;
}
```

`SignChecker` сравнивает набор подписей. Последняя подпись в массиве считается проверяемой, предыдущие подписи используются как эталонная база.

### Давление И Время

`DPoints` может хранить только координаты или расширенные данные точки:

```cpp
DPoints points;
points.AddPoint(0.0, 0.0, 0.45, 0.00, true);
points.AddPoint(1.0, 0.2, 0.60, 0.01, true);
points.AddPoint(2.0, 0.0, 0.50, 0.02, false);
```

`DtwChannel::Pressure` стоит включать, когда у всех сравниваемых точек есть pressure. `DtwChannel::Velocity` стоит включать, когда у всех точек есть timestamps. Если данных нет, канал обрабатывается безопасно и дает нулевую стоимость.

### Конфигурация

`SignCheckerConfig` управляет препроцессингом, DTW-каналами и scoring:

- `smooth`, `smooth_window` - скользящее среднее по координатам и pressure;
- `arclength_resample`, `resample_points` - равномерный ресемплинг по длине пути;
- `pca_rotate` - выравнивание подписи по главной оси;
- `zscore_normalize` - нормализация координат по осям вместо bbox-diag scale;
- `dtw_window` - радиус локального шага DTW;
- `sakoe_chiba_band` - глобальная полоса вокруг диагонали, `-1` выключает;
- `dtw_channels` - список DTW-каналов;
- `auto_weight_channels` - веса каналов по обратной внутриклассовой дисперсии;
- `trim_outliers` - медианное сравнение для simple/shape различий;
- `zscore_scoring` - альтернативное DTW scoring rule;
- `shape_weight` - вклад shape-проверок в `GetTestResult()`.

### Основные Типы

- `CTab` хранит простые признаки: `ratio_xy`, `avg_xy`, `sin_xy`, `ratio_touches`.
- `DPoint` хранит координаты, optional pressure, timestamp и pen-down state.
- `DPoints` хранит последовательность точек подписи.
- `SPoints` хранит индексы выбранных точек в `DPoints`.
- `SignCheckerConfig` задает препроцессинг, DTW-каналы и scoring.
- `SignChecker` запускает simple, DTW и shape проверки.
- `Matrix<T>` используется для матриц расстояний и промежуточных расчетов.

### Тесты

Тесты находятся в `tests/test_dtw.cpp` и проверяют:

- безопасное копирование и индексацию `Matrix`;
- поиск экстремумов в матрице;
- стабильный traceback при равных DTW-стоимостях;
- сброс результата между повторными вызовами;
- нормализацию переноса и масштаба;
- pressure, velocity и pseudo-pressure каналы;
- tangent и curvature каналы;
- fallback при пустом списке DTW-каналов;
- поведение Sakoe-Chiba band;
- arclength resampling;
- вклад shape score;
- z-score scoring;
- auto channel weighting.

Запуск:

```bash
ctest --test-dir build --output-on-failure
```

### CI/CD

GitHub Actions запускается на pull request, push в `master`/`main`, ручной запуск и теги `v*`.

CI jobs:

- `g++` и `clang++` в `Debug` и `Release`;
- AddressSanitizer и UndefinedBehaviorSanitizer;
- упаковка релизного артефакта для тегов `v*`.

Релизный артефакт содержит статическую библиотеку, заголовки и CMake package config.

### Примечания

Публичные методы `SignChecker` из исходного проекта сохранены. Новое поведение включается через `SignCheckerConfig`, кроме безопасных дефолтов вроде нормализации координат в DTW.
