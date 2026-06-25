# Duffing System: Advanced Computational Physics Final Project

**Numerical integration and chaos analysis of the Duffing oscillator using adaptive Runge-Kutta methods**

---

## 📋 Table of Contents

- [Overview](#overview)
- [The Duffing System](#the-duffing-system)
- [Repository Structure](#repository-structure)
- [Implementation Details](#implementation-details)
  - [Cash-Karp Method (Duffing.py)](#cash-karp-method-duffingpy)
  - [Dormand-Prince Method (Duffing(version2).py)](#dormand-prince-method-duffingversion2py)
  - [C++ Implementation (Duffing.cpp)](#c-implementation-duffingcpp)
- [Results](#results)
- [Algorithm Comparison](#algorithm-comparison)
- [Dependencies](#dependencies)
- [How to Run](#how-to-run)
- [References](#references)
- [Acknowledgments](#acknowledgments)

---

## 📖 Overview

This project, completed as the **final project for the Advanced Computational Physics** course, involves the numerical integration of the **Duffing oscillator**—a nonlinear, chaotic system with no analytical solution. Three implementations are provided:

1. **Python (Duffing.py)**: Cash-Karp method (Numerical Recipes, 2nd Edition)
2. **Python (Duffing(version2).py)**: Dormand-Prince method (Numerical Recipes, 3rd Edition)
3. **C++ (Duffing.cpp)**: Cash-Karp method using the Numerical Recipes library

The project analyzes the system's chaotic behavior through:
- **Phase space trajectories** (x vs. v)
- **Poincaré sections** (stroboscopic map at t = 2πn)
- **3D trajectories** (t, x, v)
- **Error analysis** and step-size evolution

---

## 🔬 The Duffing System

The Duffing oscillator is a nonlinear, driven, damped oscillator described by:

\[
\begin{align*}
\frac{dx}{dt} &= v \\
\frac{dv}{dt} &= -0.2v + x - 0.1x^3 + \cos(t) \\
\frac{ds}{dt} &= \cos\left(\frac{t}{2}\right)
\end{align*}
\]

where:
- \(x\) = position
- \(v\) = velocity
- \(t\) = time
- The third equation is an **auxiliary variable** for constructing the Poincaré section

**Initial conditions**: \(x(0) = 0\), \(v(0) = 0\), \(s(0) = 0\)  
**Integration range**: \(t \in [0, 400]\)  
**Accuracy parameter**: \(\epsilon = 10^{-4}\)

The system exhibits **chaotic behavior**, making it sensitive to initial conditions, integration parameters, and algorithm details.

---

## 📁 Repository Structure

```
duffing_system/
│
├── Duffing.py                          # Python (Cash-Karp method)
├── Duffing(version2).py                # Python (Dormand-Prince method)
├── Duffing.cpp                         # C++ (Cash-Karp, Numerical Recipes)
├── Report.pdf                          # Persian project report (detailed analysis)
└── README.md                           # This file
```

---

## 💻 Implementation Details

### Cash-Karp Method (Duffing.py)

**File**: `Duffing.py`

This implementation follows the **Cash-Karp** embedded Runge-Kutta method from the **2nd Edition of Numerical Recipes**.

**Key Features**:
- Fifth-order method with fourth-order error estimation
- Adaptive step-size control via `rkck()` and `rkqs()`
- Auxiliary variable for Poincaré section
- Extensive print statements for debugging

**Algorithm Parameters**:
```python
SAFETY = 0.9       # Safety factor for step control
PGRDW = -0.2       # Power for step growth
PSHRNK = -0.25     # Power for step shrinkage
MAXSTP = 10000     # Maximum number of steps
TINY = 1.0e-30     # Safety factor for scaling
```

**Plotting**: Generates six subplots:
1. x(t) - position vs. time
2. Phase space (x vs. v)
3. Error evolution (semilogy)
4. Step size evolution (semilogy)
5. Poincaré section (scatter plot)
6. 3D trajectory (t, x, v)

**Output**: Saves detailed output to `main.txt`

---

### Dormand-Prince Method (Duffing(version2).py)

**File**: `Duffing(version2).py`

This implementation follows the **Dormand-Prince** embedded Runge-Kutta method from the **3rd Edition of Numerical Recipes**.

**Key Features**:
- Fifth-order method with fourth-order error estimation
- More sophisticated step-size controller (PI-type)
- Dense output capability (interpolation between steps)
- Object-oriented design (Output, StepperBase, StepperDopr5, Odeint classes)

**Algorithm Parameters**:
```python
SAFE = 0.9         # Safety factor
MINSKALE = 0.2     # Minimum scale factor
MAXSKALE = 10.0    # Maximum scale factor
BETA = 0.0         # Error correction parameter
ALPHA = 0.2        # Step size adjustment exponent
MAXSTP = 50000     # Maximum number of steps
```

**Poincaré Section Method**: 
The auxiliary variable method is used to locate \(t = 2\pi n\) points through interpolation.

---

### C++ Implementation (Duffing.cpp)

**File**: `Duffing.cpp`

This implementation uses the **Numerical Recipes** library (NR) with:
- `rkck()`: Cash-Karp stepper
- `rkqs()`: Adaptive step-size control
- `odeint()`: ODE integration driver

**Key Differences from Python versions**:
- Two variables only (x and v) - no auxiliary variable for Poincaré
- Output saved to CSV for external plotting
- Uses `nr.h` header from Numerical Recipes library

**Integration Parameters**:
```cpp
kmax = 1000;        // Maximum stored points
dxsav = 0.1;        // Save points every 0.1 time units
eps = 1.0e-4;       // Accuracy parameter
h1 = 0.1;           // Initial step size
```

**Output**: Saves results to `results.csv` with columns: `t, x, v`

---

## 📊 Results

### Sample Output from `Duffing.cpp`

```
Integration complete. Good steps: 3803, Bad steps: 252
Final x: -2.0746167, v: 2.127398
Results saved to results.csv
```

### Sample Output from `Duffing.py`

```
Integration completed!
Steps taken: 10640
Successful steps: 0, Failed steps: 0
Final position: -2.074582, Final velocity: 2.127327, Final s: 0.000117
Poincaré points: 63
```

---

## 📈 Algorithm Comparison

### Cash-Karp (Duffing.py) vs. Dormand-Prince (Duffing(version2).py)

| Feature | Cash-Karp (v1) | Dormand-Prince (v2) |
|---------|----------------|---------------------|
| **Method Order** | 5(4) | 5(4) |
| **Error Estimation** | Embedded RK | Embedded RK |
| **Step Controller** | Simple | Advanced (PI controller) |
| **Dense Output** | No | Yes (interpolation) |
| **Error Scaling** | `yscal[i] = |y[i]| + |dydx[i]*h| + TINY` | `scale = atol + rtol * max(|y|, |yout|)` |
| **Step Growth** | `hnext = SAFETY * h * errmax^PGRDW` | `scale = SAFE * err^(-ALPHA) * errold^BETA` |
| **Code Structure** | Procedural | Object-oriented |

### Key Differences in Error Calculation

**Version 1 (Cash-Karp)**:
```python
yscal[i] = abs(y[i]) + abs(dydx[i] * h) + TINY
errmax = max(errmax, abs(yerr[i] / yscal[i]))
errmax /= eps
```

**Version 2 (Dormand-Prince)**:
```python
scale = atol + rtol * max(abs(y[i]), abs(yout[i]))
err += (yerr[i] / scale)**2
err = sqrt(err / n)
```

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

### Python (Version 1 - Cash-Karp)

```bash
python Duffing.py
```

This will:
1. Integrate the Duffing system using the Cash-Karp method
2. Save output to `main.txt`
3. Display six plots (x(t), phase space, error, step size, Poincaré, 3D)

### Python (Version 2 - Dormand-Prince)

```bash
python "Duffing(version2).py"
```

This will:
1. Integrate the Duffing system using the Dormand-Prince method
2. Display six plots (x(t), phase space, error, step size, Poincaré, 3D)

### C++

```bash
# Compile
g++ Duffing.cpp -o duffing -lnr

# Run
./duffing
```

This will:
1. Integrate the Duffing system using the Cash-Karp method
2. Save results to `results.csv`
3. Print integration summary to console

### Plotting C++ Results (Python)

```python
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('results.csv')
plt.plot(data['t'], data['x'])
plt.xlabel('t')
plt.ylabel('x')
plt.show()
```

---

## 📚 References

1. **Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P.** (1992). *Numerical Recipes in C: The Art of Scientific Computing*, 2nd Edition. Cambridge University Press.

2. **Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P.** (2007). *Numerical Recipes: The Art of Scientific Computing*, 3rd Edition. Cambridge University Press.

3. **Duffing, G.** (1918). *Erzwungene Schwingungen bei veränderlicher Eigenfrequenz*. Vieweg.

4. **Strogatz, S. H.** (2018). *Nonlinear Dynamics and Chaos: With Applications to Physics, Biology, Chemistry, and Engineering*, 2nd Edition. CRC Press.

---

## 🙏 Acknowledgments

This project was completed as the **final project for the Advanced Computational Physics** course. The goal was to implement and analyze adaptive Runge-Kutta methods for a chaotic system, comparing different algorithmic approaches and studying parameter sensitivity. The detailed sensitivity analysis is documented in `Report.pdf`.

---

**Happy Computing!** ⚡🌀
