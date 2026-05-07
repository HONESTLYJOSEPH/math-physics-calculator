# 🧮 University Math & Physics Calculator

Two command-line projects — one in **C**, one in **Python** — covering university-level math and physics problems.
✨😊

---

## 📚 Topics Covered

| Topic | Tools |
|---|---|
| **Algebra** | Quadratic equations, linear systems, polynomial roots |
| **Calculus** | Derivatives, integrals, root finding, ODE solver |
| **Mechanics** | Kinematics, projectile motion, circular motion, energy |
| **Thermodynamics** | Ideal gas law, heat transfer, Carnot efficiency |
| **Statistics** | Descriptive stats, hypothesis testing, linear regression |

---

## ⚙️ C Project

**File:** `math_physics_calc.c`

Zero external dependencies — just a C compiler.

**Compile & run:**
```bash
make && ./math_physics_calc
```

Or manually:
```bash
gcc -Wall -O2 -o math_physics_calc math_physics_calc.c -lm
./math_physics_calc
```

**Requirements:** Any C11 compiler (gcc, clang, MSVC)

---

## 🐍 Python Project

**File:** `math_physics_toolkit.py`

Includes automatic plots for trajectories, histograms, regression lines, and ODE solutions.

**Install dependencies:**
```bash
pip install -r requirements.txt
```

**Run:**
```bash
python math_physics_toolkit.py
```

**Requirements:** Python 3.8+, numpy, scipy, matplotlib

---

## 📁 Files

```
├── math_physics_calc.c       # C project source code
├── Makefile                  # Build file for C project
├── math_physics_toolkit.py   # Python project
├── requirements.txt          # Python dependencies
└── README.md                 # This file
```

---

## 🚀 Quick Start (Python)

```bash
git clone https://github.com/YOUR_USERNAME/math-physics-calculator
cd math-physics-calculator
pip install -r requirements.txt
python math_physics_toolkit.py
```
