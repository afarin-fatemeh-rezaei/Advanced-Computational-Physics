# Advanced Computational Physics

**Computational physics exercises and final project: Random walks, Monte Carlo methods, and chaotic systems**

---

## 📋 Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Exercises](#exercises)
  - [1. Probability Distribution](#1-probability-distribution)
  - [2. Random Walk (Single Particle)](#2-random-walk-single-particle)
  - [3. Random Walk Ensemble](#3-random-walk-ensemble)
  - [4. Acceptance-Rejection Sampling](#4-acceptance-rejection-sampling)
  - [5. Monte Carlo π Estimate](#5-monte-carlo-π-estimate)
- [Final Project: Duffing System](#final-project-duffing-system)
- [Dependencies](#dependencies)
- [How to Run](#how-to-run)
- [References](#references)
- [Acknowledgments](#acknowledgments)

---

## 📖 Overview

This repository contains coursework for an **Advanced Computational Physics** course, covering fundamental computational methods in physics:

- **Random walks** and statistical ensembles
- **Monte Carlo methods** for integration and sampling
- **Acceptance-rejection sampling** for probability distributions
- **Chaos theory** and numerical integration of nonlinear systems

The materials are organized into two main sections:
1. **Exercises** – Standalone Python scripts exploring core computational methods
2. **Final Project** – Numerical integration and chaos analysis of the Duffing oscillator

---

## 📁 Repository Structure

```
Advanced_Computational_Physics/
│
├── Probability Distribution.py          # Uniform random distribution analysis
├── Random Walk.py                       # Single particle random walk
├── Random Walk Ensemble.py              # Ensemble comparison (single vs. N-particle)
├── Acceptance-Rejection Sampling.py     # Sampling from arbitrary distributions
├── Monte Carlo π Estimate.py            # π estimation via Monte Carlo integration
│
├── duffing_system/                      # Final project folder
│   ├── Duffing.py                       # Python (Cash-Karp method)
│   ├── Duffing(version2).py             # Python (Dormand-Prince method)
│   ├── Duffing.cpp                      # C++ (Numerical Recipes)
│   ├── Report.pdf                       # Persian project report
│   └── README.md                        # Project documentation
│
└── README.md                            # This file
```

---

## 🧪 Exercises

### 1. Probability Distribution

**File**: `Probability Distribution.py`

**Objective**: Generate and analyze a uniform random distribution.

**Methodology**:
- Generate 1000 random numbers uniformly distributed between 0 and 1
- Bin the data into 10 equal-width bins
- Calculate probabilities for each bin
- Compute the variance of the probability distribution
- Compare variance with \(1/\sqrt{N}\)

**Key Results**:
```python
Probabilities for each bin:
Bin 1 (0.0-0.1): 0.1020
Bin 2 (0.1-0.2): 0.1010
...
Variance: 0.0095
1/sqrt(1000): 0.0316
```

**Concepts Explored**:
- Uniform random number generation
- Histogram binning and probability distributions
- Statistical variance and convergence with sample size

**How to Run**:
```bash
python "Probability Distribution.py"
```

---

### 2. Random Walk (Single Particle)

**File**: `Random Walk.py`

**Objective**: Simulate a single-particle random walk in one dimension.

**Methodology**:
- Generate \(M = 1000\) steps, each either +1 or -1 (equal probability)
- Calculate the time average: \(\langle s \rangle = \frac{1}{M}\sum_{i=1}^M s_i\)
- Plot the cumulative position as a function of step number

**Also includes**:
- Ensemble of 100 particles, each taking 10 steps
- Bar plot showing final positions of all 100 particles
- Comparison with the single-particle time average

**Key Results**:
```
Single Particle M-step (M=1000) Time Average: -0.002
Average of Final Positions for N=100 Particles (10 steps each): -0.020
```

**Concepts Explored**:
- Random walk theory and Brownian motion
- Time averages vs. ensemble averages
- Law of large numbers
- Central limit theorem

**How to Run**:
```bash
python "Random Walk.py"
```

---

### 3. Random Walk Ensemble

**File**: `Random Walk Ensemble.py`

**Objective**: Compare single-particle time averages with ensemble averages for increasing system sizes.

**Methodology**:
- For each total number of steps \(M\) [100, 1000, 10000, 100000, 1000000, 10000000]:
  - **Single particle**: Generate \(M\) steps, calculate time average
  - **Ensemble**: Partition into \(N = M/P\) particles, each taking \(P = 10\) steps, average the final positions
- Calculate absolute error and RMSE between the two methods
- Plot results on log scale

**Key Results**:
| M | Time Average | Ensemble Average | Absolute Error | RMSE |
|---|--------------|------------------|----------------|------|
| 100 | 0.040000 | 0.000000 | 0.040000 | 0.040000 |
| 1000 | -0.060000 | 0.100000 | 0.160000 | 0.160000 |
| 10000 | 0.013000 | -0.010000 | 0.023000 | 0.023000 |
| ... | ... | ... | ... | ... |

**Concepts Explored**:
- Ergodic hypothesis: time averages vs. ensemble averages
- Statistical convergence with system size
- Error scaling with \(1/\sqrt{N}\)

**How to Run**:
```bash
python "Random Walk Ensemble.py"
```

---

### 4. Acceptance-Rejection Sampling

**File**: `Acceptance-Rejection Sampling.py`

**Objective**: Sample from a non-standard probability distribution using the acceptance-rejection method.

**Methodology**:
- Target distribution: \(p(x) = e^{-(x-1)^2}\) (unnormalized Gaussian centered at \(x=1\))
- Normalize \(p(x)\) numerically using trapezoidal integration
- Sample from \(p(x)\) using acceptance-rejection:
  - Generate \(y \sim \text{Uniform}(-2, 4)\)
  - Generate \(u \sim \text{Uniform}(0, 1)\)
  - Accept \(y\) if \(u < p(y)\)
- Generate 100, 10,000, and 1,000,000 samples
- Compare histograms with true normalized PDF
- Calculate Mean Absolute Error (MAE) and Root Mean Squared Error (RMSE)

**Key Results**:
```
Error Analysis:
Sample Size   | MAE      | RMSE
----------------------------------------
100           | 0.0154   | 0.0224
10,000        | 0.0016   | 0.0023
1,000,000     | 0.0002   | 0.0002
```

**Concepts Explored**:
- Acceptance-rejection sampling algorithm
- Numerical integration (trapezoidal rule)
- Convergence of sampling with increasing samples
- Error metrics (MAE, RMSE)

**How to Run**:
```bash
python "Acceptance-Rejection Sampling.py"
```

---

### 5. Monte Carlo π Estimate

**File**: `Monte Carlo π Estimate.py`

**Objective**: Estimate the value of π using Monte Carlo integration.

**Methodology**:
- Generate random points uniformly in the square \([-1, 1] \times [-1, 1]\)
- Count points inside the unit circle: \(x^2 + y^2 \leq 1\)
- Area estimate: \(A = 4 \times \frac{\text{points inside}}{\text{total points}}\)
- Compare to true value \(\pi\)
- Test sample sizes: 10, 100, 1000, 10⁴, 10⁵, 10⁶

**Key Results**:
```
Samples     Estimate    Error
10          3.60000     0.45841
100         3.12000     0.02159
1000        3.15600     0.01441
10000       3.13800     0.00359
100000      3.14124     0.00035
1000000     3.14114     0.00045
```

**Concepts Explored**:
- Monte Carlo integration
- Error scaling with sample size (\(\propto 1/\sqrt{N}\))
- Geometric probability

**How to Run**:
```bash
python "Monte Carlo π Estimate.py"
```

---

## 🚀 Final Project: Duffing System

**Folder**: `duffing_system/`

**Objective**: Numerically integrate and analyze the chaotic Duffing oscillator using adaptive Runge-Kutta methods.

### The Duffing System

\[
\begin{align*}
\frac{dx}{dt} &= v \\
\frac{dv}{dt} &= -0.2v + x - 0.1x^3 + \cos(t) \\
\frac{ds}{dt} &= \cos\left(\frac{t}{2}\right)
\end{align*}
\]

### Implementations

| File | Method | Language |
|------|--------|----------|
| `Duffing.py` | Cash-Karp (5th order) | Python |
| `Duffing(version2).py` | Dormand-Prince (5th order) | Python |
| `Duffing.cpp` | Cash-Karp (Numerical Recipes) | C++ |

### Key Features

- **Adaptive step-size control** for accuracy and efficiency
- **Poincaré sections** at \(t = 2\pi n\)
- **Phase space analysis** (x vs. v)
- **Error analysis** and step-size evolution
- **3D trajectory visualization** (t, x, v)

### Documentation

See `duffing_system/README.md` for:
- Detailed methodology
- Algorithm comparison
- Sensitivity analysis
- How to run each implementation

---

## 📦 Dependencies

### Python
- Python 3.x
- NumPy
- Matplotlib

### C++
- C++ compiler (g++, clang++)
- Numerical Recipes library (`nr.h`)

### Installation

```bash
# Python dependencies
pip install numpy matplotlib
```

---

## 🚀 How to Run

### Exercises

```bash
# Probability distribution
python "Probability Distribution.py"

# Random walk (single particle)
python "Random Walk.py"

# Random walk ensemble comparison
python "Random Walk Ensemble.py"

# Acceptance-rejection sampling
python "Acceptance-Rejection Sampling.py"

# Monte Carlo π estimate
python "Monte Carlo π Estimate.py"
```

### Final Project

```bash
# Python (Cash-Karp)
cd duffing_system
python Duffing.py

# Python (Dormand-Prince)
python "Duffing(version2).py"

# C++
g++ Duffing.cpp -o duffing -lnr
./duffing
```

---

## 📚 References

1. **Metropolis, N., & Ulam, S.** (1949). "The Monte Carlo Method." *Journal of the American Statistical Association*, 44(247), 335-341.

2. **Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P.** (1992). *Numerical Recipes in C: The Art of Scientific Computing*, 2nd Edition. Cambridge University Press.

3. **Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P.** (2007). *Numerical Recipes: The Art of Scientific Computing*, 3rd Edition. Cambridge University Press.

4. **Strogatz, S. H.** (2018). *Nonlinear Dynamics and Chaos: With Applications to Physics, Biology, Chemistry, and Engineering*, 2nd Edition. CRC Press.

5. **Landau, R. H., Páez, M. J., & Bordeianu, C. C.** (2015). *Computational Physics: Problem Solving with Python*, 3rd Edition. Wiley-VCH.

---

## 🙏 Acknowledgments

This repository contains exercises and the final project for the **Advanced Computational Physics** course. The exercises explore fundamental computational methods, while the final project applies these methods to a chaotic dynamical system.

---

**Happy Computing!** ⚡🧮
