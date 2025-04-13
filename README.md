<!-- markdownlint-disable MD024 -->

# DTW

[English](#english) | [Русский](#русский)

---

## English

A small C++17 library for handwritten signature comparison. The project combines simple signature features with checks based on Dynamic Time Warping (DTW).

### Features

- computation of simple signature features via `CTab`;
- comparison of signature point trajectories via DTW;
- distance, speed, and direction metrics;
- safe `Matrix<T>` template for internal computations;
- CMake build, tests via CTest and GitHub Actions;
- library installation with CMake package config for `find_package(DTW)`.

### Layout

```text
.
├── Matrix.h                 # matrix template
├── SignatureData.h/.cpp     # signature points, simple features and metrics
├── SignChecker.h/.cpp       # main signature verification logic
├── tests/                   # regression tests
├── cmake/                   # CMake package config template
└── .github/workflows/ci.yml # CI and release artifact packaging
```

### Requirements

- CMake 3.16 or newer;
- a C++17-capable compiler, e.g. `g++` or `clang++`.

### Quick start

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

For a Release build:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

### Installation

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
cmake --install build --prefix dist
```

After installation, the library can be picked up from another CMake project:

```cmake
find_package(DTW REQUIRED)

add_executable(app main.cpp)
target_link_libraries(app PRIVATE DTW::dtw)
```

If the install lives in a local directory:

```bash
cmake -S your-app -B your-app/build -DCMAKE_PREFIX_PATH=/path/to/DTW/dist
```

### Usage example

```cpp
#include "SignChecker.h"

int main()
{
    balance = 0.0; // 0.0: DTW only, 1.0: simple checks only

    DPoints signatures[2] = {
        DPoints{{0.0, 0.0}, {1.0, 1.0}, {2.0, 0.0}},
        DPoints{{0.0, 0.0}, {1.0, 1.0}, {2.0, 0.0}},
    };

    SPoints sampled[2] = {
        SPoints::FullRange(signatures[0].GetN()),
        SPoints::FullRange(signatures[1].GetN()),
    };

    SignChecker checker;
    if (!checker.DTWCheckForSimpleForge(signatures, sampled, 2))
        return 1;

    const double score = checker.GetTestResult();
    return score >= 0.9 ? 0 : 2;
}
```

`SignChecker` compares a set of signatures, where the last signature in the array is the one being verified and the previous ones are used as the reference base.

### Core types

- `CTab` stores simple features: `ratio_xy`, `avg_xy`, `sin_xy`, `ratio_touches`.
- `DPoints` stores signature point coordinates.
- `SPoints` stores indices of selected points in `DPoints`.
- `SignChecker` runs the simple and DTW checks.
- `Matrix<T>` is used for distance matrices and intermediate computations.

### Tests

The tests live in `tests/test_dtw.cpp` and cover:

- safe copying and indexing of `Matrix`;
- locating extrema in a matrix;
- stable traceback under equal DTW costs;
- absence of result accumulation across repeated calls;
- the baseline scenario of identical signatures.

To run:

```bash
ctest --test-dir build --output-on-failure
```

### CI/CD

The GitHub Actions workflow runs on pull requests, pushes to `master`/`main`, manual dispatch, and `v*` tags.

The CI matrix builds and tests the project in `Debug` and `Release` on:

- `g++`;
- `clang++`.

For `v*`-style tags, a release artifact is additionally produced via:

```bash
cmake --install build --prefix dist
```

The artifact contains the static library, headers, and CMake package config.

### Current API status

The project preserves the original public `SignChecker` interface but is now built as a standalone library. The global parameters `balance`, `Eps`, and `agree_param` are declared in `SignatureData.h` and defined in `SignatureData.cpp`.

---

## Русский

Небольшая C++17-библиотека для сравнения рукописных подписей. Проект сочетает простые признаки подписи с проверками на основе Dynamic Time Warping (DTW).

### Возможности

- расчет простых признаков подписи через `CTab`;
- сравнение траекторий точек подписи через DTW;
- метрики расстояния, скорости и направления;
- безопасная матрица `Matrix<T>` для внутренних расчетов;
- CMake-сборка, тесты через CTest и GitHub Actions;
- установка библиотеки с CMake package config для `find_package(DTW)`.

### Структура

```text
.
├── Matrix.h                 # шаблон матрицы
├── SignatureData.h/.cpp     # точки подписи, простые признаки и метрики
├── SignChecker.h/.cpp       # основная логика проверки подписи
├── tests/                   # регрессионные тесты
├── cmake/                   # CMake package config template
└── .github/workflows/ci.yml # CI и упаковка релизного артефакта
```

### Требования

- CMake 3.16 или новее;
- компилятор с поддержкой C++17, например `g++` или `clang++`.

### Быстрый старт

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Для Release-сборки:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

### Установка

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
cmake --install build --prefix dist
```

После установки библиотеку можно найти из другого CMake-проекта:

```cmake
find_package(DTW REQUIRED)

add_executable(app main.cpp)
target_link_libraries(app PRIVATE DTW::dtw)
```

Если установка лежит в локальной директории:

```bash
cmake -S your-app -B your-app/build -DCMAKE_PREFIX_PATH=/path/to/DTW/dist
```

### Пример использования

```cpp
#include "SignChecker.h"

int main()
{
    balance = 0.0; // 0.0: только DTW, 1.0: только простые проверки

    DPoints signatures[2] = {
        DPoints{{0.0, 0.0}, {1.0, 1.0}, {2.0, 0.0}},
        DPoints{{0.0, 0.0}, {1.0, 1.0}, {2.0, 0.0}},
    };

    SPoints sampled[2] = {
        SPoints::FullRange(signatures[0].GetN()),
        SPoints::FullRange(signatures[1].GetN()),
    };

    SignChecker checker;
    if (!checker.DTWCheckForSimpleForge(signatures, sampled, 2))
        return 1;

    const double score = checker.GetTestResult();
    return score >= 0.9 ? 0 : 2;
}
```

`SignChecker` сравнивает набор подписей, где последняя подпись в массиве считается проверяемой, а предыдущие используются как эталонная база.

### Основные типы

- `CTab` хранит простые признаки: `ratio_xy`, `avg_xy`, `sin_xy`, `ratio_touches`.
- `DPoints` хранит координаты точек подписи.
- `SPoints` хранит индексы выбранных точек в `DPoints`.
- `SignChecker` запускает простые проверки и DTW-проверки.
- `Matrix<T>` используется для матриц расстояний и промежуточных расчетов.

### Тесты

Тесты находятся в `tests/test_dtw.cpp` и проверяют:

- безопасное копирование и индексацию `Matrix`;
- поиск экстремумов в матрице;
- стабильный traceback при равных DTW-стоимостях;
- отсутствие накопления результата между повторными вызовами;
- базовый сценарий идентичных подписей.

Запуск:

```bash
ctest --test-dir build --output-on-failure
```

### CI/CD

GitHub Actions workflow запускается на pull request, push в `master`/`main`, ручной запуск и теги `v*`.

Матрица CI собирает и тестирует проект в `Debug` и `Release` на:

- `g++`;
- `clang++`.

Для тегов вида `v*` дополнительно собирается релизный артефакт через:

```bash
cmake --install build --prefix dist
```

Артефакт содержит статическую библиотеку, заголовки и CMake package config.

### Текущее состояние API

Проект сохраняет исходный публичный интерфейс `SignChecker`, но теперь собирается как самостоятельная библиотека. Глобальные параметры `balance`, `Eps` и `agree_param` объявлены в `SignatureData.h` и определены в `SignatureData.cpp`.
