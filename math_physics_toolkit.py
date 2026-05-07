#!/usr/bin/env python3
"""
================================================================
  University Math & Physics Toolkit  —  Python Edition
  Topics: Algebra, Calculus, Mechanics, Thermodynamics,
          Statistics
  Dependencies: numpy, scipy, matplotlib (all standard)
================================================================
"""

import sys
import math
import cmath
import statistics as st

# ── Try importing scientific stack ──────────────────────────
try:
    import numpy as np
    import scipy.optimize as opt
    import scipy.integrate as integrate
    import scipy.stats as stats
    import matplotlib.pyplot as plt
    SCI = True
except ImportError:
    SCI = False
    print("⚠  numpy / scipy / matplotlib not found — running in "
          "pure-Python fallback mode (some features limited).\n"
          "   Install with:  pip install numpy scipy matplotlib\n")

# ── Console colours ─────────────────────────────────────────
CYAN  = "\033[36m"; GREEN = "\033[32m"; YELLOW = "\033[33m"
RED   = "\033[31m"; BOLD  = "\033[1m";  RESET  = "\033[0m"
GREY  = "\033[90m"

def section(title: str) -> None:
    print(f"\n{CYAN}{'─'*50}\n  {title}\n{'─'*50}{RESET}")

def ok(*args):  print(GREEN, *args, RESET)
def err(*args): print(RED,   *args, RESET)

# ============================================================
# 1.  ALGEBRA
# ============================================================

def quadratic_solver():
    section("Quadratic Solver  ax² + bx + c = 0")
    a = float(input("  a = "))
    b = float(input("  b = "))
    c = float(input("  c = "))

    if abs(a) < 1e-12:
        if abs(b) < 1e-12:
            err("  Degenerate — no variable terms.")
        else:
            ok(f"  Linear root: x = {-c/b:.6f}")
        return

    disc = b**2 - 4*a*c
    if disc >= 0:
        sq = math.sqrt(disc)
        x1, x2 = (-b + sq)/(2*a), (-b - sq)/(2*a)
        ok(f"  Real roots:\n    x₁ = {x1:.6f}\n    x₂ = {x2:.6f}")
    else:
        r = complex(-b, math.sqrt(-disc)) / (2*a)
        ok(f"  Complex roots:\n    x₁ = {r.real:.6f} + {r.imag:.6f}i\n"
           f"    x₂ = {r.real:.6f} - {r.imag:.6f}i")


def linear_system_solver():
    section("Linear System Solver  Ax = b")
    n = int(input("  Number of equations (1–8): "))
    print("  Enter augmented matrix [A|b], one row at a time:")
    rows = []
    for i in range(n):
        row = list(map(float, input(f"  Row {i+1}: ").split()))
        if len(row) != n + 1:
            err("  Wrong number of values."); return
        rows.append(row)

    if SCI:
        A = np.array([r[:n] for r in rows])
        b_vec = np.array([r[n] for r in rows])
        try:
            x = np.linalg.solve(A, b_vec)
            ok("  Solution:")
            for i, xi in enumerate(x):
                ok(f"    x{i+1} = {xi:.8f}")
            # Residual
            res = np.linalg.norm(A @ x - b_vec)
            print(GREY + f"    ‖Ax - b‖ = {res:.2e}" + RESET)
        except np.linalg.LinAlgError:
            err("  Singular matrix — no unique solution.")
    else:
        # Pure-Python Gaussian elimination
        M = [rows[i][:] for i in range(n)]
        for col in range(n):
            pivot = max(range(col, n), key=lambda r: abs(M[r][col]))
            M[col], M[pivot] = M[pivot], M[col]
            if abs(M[col][col]) < 1e-12:
                err("  Singular."); return
            for row in range(col+1, n):
                f = M[row][col] / M[col][col]
                for j in range(col, n+1): M[row][j] -= f * M[col][j]
        x = [0.0]*n
        for i in range(n-1, -1, -1):
            x[i] = M[i][n]
            for j in range(i+1, n): x[i] -= M[i][j]*x[j]
            x[i] /= M[i][i]
        ok("  Solution:")
        for i, xi in enumerate(x): ok(f"    x{i+1} = {xi:.8f}")


def polynomial_roots():
    """Find all roots of a polynomial using numpy."""
    if not SCI:
        err("  Requires numpy."); return
    section("Polynomial Root Finder  (numpy)")
    deg = int(input("  Polynomial degree: "))
    print(f"  Enter {deg+1} coefficients (highest to lowest power):")
    coeffs = list(map(float, input("  > ").split()))
    roots = np.roots(coeffs)
    ok("  Roots:")
    for i, r in enumerate(roots):
        if abs(r.imag) < 1e-8:
            ok(f"    r{i+1} = {r.real:.8f}")
        else:
            ok(f"    r{i+1} = {r.real:.6f} {'+' if r.imag>=0 else '-'} {abs(r.imag):.6f}i")

# ============================================================
# 2.  CALCULUS
# ============================================================

_FUNCTIONS = {
    "1": ("sin(x)",        lambda x: math.sin(x)),
    "2": ("cos(x)",        lambda x: math.cos(x)),
    "3": ("e^x",           lambda x: math.exp(x)),
    "4": ("ln|x|",         lambda x: math.log(abs(x))),
    "5": ("x³ - 3x² + 2x - 1", lambda x: x**3 - 3*x**2 + 2*x - 1),
    "6": ("sin(x)·e^(-x)", lambda x: math.sin(x)*math.exp(-x)),
    "7": ("x²·cos(x)",     lambda x: x**2 * math.cos(x)),
}

def pick_function():
    print("  Functions:")
    for k, (name, _) in _FUNCTIONS.items():
        print(f"    {k}) {name}")
    ch = input("  Choice: ").strip()
    if ch not in _FUNCTIONS:
        print(GREY + "  Defaulting to sin(x)" + RESET)
        return _FUNCTIONS["1"]
    return _FUNCTIONS[ch]


def numerical_derivative():
    section("Numerical Derivative — Central Difference & Exact (scipy)")
    name, f = pick_function()
    x = float(input("  Evaluate f'(x) at x = "))
    h = 1e-6
    cd = (f(x+h) - f(x-h)) / (2*h)
    ok(f"  f(x)  = {name}")
    ok(f"  f'({x}) ≈ {cd:.10f}  (central difference, h={h})")
    if SCI:
        from scipy.misc import derivative as spd
        exact = spd(f, x, dx=1e-6, n=1)
        ok(f"  f'({x}) ≈ {exact:.10f}  (scipy.misc.derivative)")


def numerical_integral():
    section("Numerical Integration — Simpson's Rule & Adaptive (scipy)")
    name, f = pick_function()
    a = float(input("  Lower limit a = "))
    b = float(input("  Upper limit b = "))

    # Pure Simpson
    n = 1000
    h = (b - a) / n
    S = f(a) + f(b)
    for i in range(1, n):
        S += (4 if i % 2 else 2) * f(a + i*h)
    simp = S * h / 3
    ok(f"  ∫[{a:.3f}→{b:.3f}] {name} dx")
    ok(f"    Simpson (n={n}): {simp:.10f}")
    if SCI:
        result, error = integrate.quad(f, a, b)
        ok(f"    scipy.integrate.quad: {result:.10f}  (est. error {error:.2e})")


def ode_solver():
    """Solve dy/dt = f(t, y) with scipy."""
    if not SCI:
        err("  Requires scipy."); return
    section("ODE Solver  dy/dt = f(t, y)  [scipy.solve_ivp]")
    print("  Built-in equations:")
    print("    1) Simple decay:  dy/dt = -ky")
    print("    2) Logistic:       dy/dt = r·y·(1 - y/K)")
    print("    3) Harmonic osc:  dy/dt = v,  dv/dt = -ω²y")
    ch = input("  Choice: ").strip()

    if ch == "1":
        k = float(input("  k = "))
        y0 = float(input("  y(0) = "))
        t_end = float(input("  t_end = "))
        sol = integrate.solve_ivp(lambda t,y: [-k*y[0]], [0, t_end],
                                  [y0], dense_output=True)
        t_pts = np.linspace(0, t_end, 400)
        y_pts = sol.sol(t_pts)[0]
        ok(f"  Final value y({t_end:.2f}) = {y_pts[-1]:.6f}")
        plt.figure(figsize=(8,4))
        plt.plot(t_pts, y_pts, color="#00b4d8", lw=2)
        plt.title(f"Exponential Decay  dy/dt = -{k}y"); plt.xlabel("t"); plt.ylabel("y")
        plt.grid(alpha=0.3); plt.tight_layout(); plt.show()

    elif ch == "2":
        r = float(input("  Growth rate r = "))
        K = float(input("  Carrying capacity K = "))
        y0 = float(input("  y(0) = "))
        t_end = float(input("  t_end = "))
        sol = integrate.solve_ivp(lambda t,y: [r*y[0]*(1-y[0]/K)],
                                  [0, t_end], [y0], dense_output=True)
        t_pts = np.linspace(0, t_end, 400)
        y_pts = sol.sol(t_pts)[0]
        ok(f"  y({t_end:.2f}) = {y_pts[-1]:.6f}")
        plt.figure(figsize=(8,4))
        plt.plot(t_pts, y_pts, color="#90e0ef", lw=2)
        plt.axhline(K, ls="--", color="#ff6b6b", label=f"K={K}")
        plt.title("Logistic Growth"); plt.xlabel("t"); plt.ylabel("y")
        plt.legend(); plt.grid(alpha=0.3); plt.tight_layout(); plt.show()

    elif ch == "3":
        omega = float(input("  Angular frequency ω = "))
        y0_disp = float(input("  y(0) displacement = "))
        y0_vel  = float(input("  y'(0) velocity    = "))
        t_end = float(input("  t_end = "))
        def harm(t, Y): return [Y[1], -omega**2 * Y[0]]
        sol = integrate.solve_ivp(harm, [0, t_end], [y0_disp, y0_vel], dense_output=True)
        t_pts = np.linspace(0, t_end, 800)
        y_pts = sol.sol(t_pts)
        ok(f"  y({t_end:.2f}) = {y_pts[0,-1]:.6f}")
        plt.figure(figsize=(8,5))
        plt.plot(t_pts, y_pts[0], label="displacement y", lw=2)
        plt.plot(t_pts, y_pts[1], label="velocity v", lw=2)
        plt.title("Harmonic Oscillator"); plt.xlabel("t")
        plt.legend(); plt.grid(alpha=0.3); plt.tight_layout(); plt.show()
    else:
        err("  Unknown choice.")


# ============================================================
# 3.  MECHANICS
# ============================================================

def kinematics_1d():
    section("1-D Kinematics  (constant acceleration)")
    u = float(input("  Initial velocity u (m/s)  = "))
    a_acc = float(input("  Acceleration     a (m/s²) = "))
    t = float(input("  Time             t (s)    = "))
    v = u + a_acc * t
    s = u*t + 0.5*a_acc*t**2
    ok(f"  v = {v:.4f} m/s")
    ok(f"  s = {s:.4f} m")
    ok(f"  v² = {v**2:.4f} m²/s²")


def projectile_motion():
    section("Projectile Motion  (flat Earth, no drag)")
    v0    = float(input("  Launch speed v₀ (m/s)    = "))
    angle = float(input("  Launch angle θ  (deg)    = "))
    h0    = float(input("  Initial height h₀ (m)    = "))
    g     = 9.80665

    theta = math.radians(angle)
    vx = v0 * math.cos(theta)
    vy = v0 * math.sin(theta)

    disc = vy**2 + 2*g*h0
    if disc < 0:
        err("  No real solution."); return
    t_flight = (vy + math.sqrt(disc)) / g
    t_peak   = vy / g
    h_peak   = h0 + vy*t_peak - 0.5*g*t_peak**2 if t_peak >= 0 else h0
    rng      = vx * t_flight

    ok(f"  Time of flight  = {t_flight:.4f} s")
    ok(f"  Range           = {rng:.4f} m")
    ok(f"  Peak height     = {h_peak:.4f} m")

    if SCI:
        t_pts = np.linspace(0, t_flight, 500)
        x_pts = vx * t_pts
        y_pts = h0 + vy*t_pts - 0.5*g*t_pts**2
        plt.figure(figsize=(8,5))
        plt.plot(x_pts, y_pts, color="#48cae4", lw=2)
        plt.axhline(0, color="#aaa", lw=0.8)
        plt.title("Projectile Trajectory")
        plt.xlabel("Horizontal distance (m)"); plt.ylabel("Height (m)")
        plt.grid(alpha=0.3); plt.tight_layout(); plt.show()


def energy_work():
    section("Energy & Work")
    print("  1) Kinetic energy   KE = ½mv²")
    print("  2) Gravitational PE = mgh")
    print("  3) Work-energy theorem")
    ch = input("  Choice: ").strip()
    g = 9.80665
    if ch == "1":
        m = float(input("  Mass m (kg) = "))
        v = float(input("  Speed v (m/s) = "))
        ok(f"  KE = {0.5*m*v**2:.4f} J")
    elif ch == "2":
        m = float(input("  Mass m (kg) = "))
        h = float(input("  Height h (m) = "))
        ok(f"  PE = {m*g*h:.4f} J")
    elif ch == "3":
        m  = float(input("  Mass m (kg) = "))
        v1 = float(input("  Initial speed v₁ (m/s) = "))
        v2 = float(input("  Final speed   v₂ (m/s) = "))
        W  = 0.5*m*(v2**2 - v1**2)
        ok(f"  Net work done W = {W:.4f} J")


# ============================================================
# 4.  THERMODYNAMICS
# ============================================================

R_GAS  = 8.314     # J/(mol·K)
R_LITATM = 0.082057  # L·atm/(mol·K)

def ideal_gas_law():
    section("Ideal Gas Law  PV = nRT")
    print("  Solve for: [P]ressure / [V]olume / [n]oles / [T]emperature")
    var = input("  Variable: ").strip().upper()
    if var == "P":
        n = float(input("  n (mol) = ")); V = float(input("  V (L) = "))
        T = float(input("  T (K) = "))
        ok(f"  P = {n*R_LITATM*T/V:.4f} atm   ({n*R_GAS*T/(V/1000):.2f} Pa)")
    elif var == "V":
        n = float(input("  n (mol) = ")); P = float(input("  P (atm) = "))
        T = float(input("  T (K) = "))
        ok(f"  V = {n*R_LITATM*T/P:.4f} L")
    elif var == "N":
        P = float(input("  P (atm) = ")); V = float(input("  V (L) = "))
        T = float(input("  T (K) = "))
        ok(f"  n = {P*V/(R_LITATM*T):.4f} mol")
    elif var == "T":
        P = float(input("  P (atm) = ")); V = float(input("  V (L) = "))
        n = float(input("  n (mol) = "))
        ok(f"  T = {P*V/(n*R_LITATM):.4f} K")
    else:
        err("  Unrecognised variable.")


def heat_transfer():
    section("Sensible Heat Transfer  Q = mcΔT")
    m  = float(input("  Mass m (kg) = "))
    c  = float(input("  Specific heat c (J/kg·K) = "))
    dT = float(input("  ΔT (K or °C) = "))
    Q  = m * c * dT
    ok(f"  Q = {Q:.4f} J  ({Q/1000:.4f} kJ)  ({Q/4184:.4f} kcal)")


def carnot_efficiency():
    section("Carnot Cycle Analysis")
    T_H = float(input("  Hot reservoir  T_H (K) = "))
    T_C = float(input("  Cold reservoir T_C (K) = "))
    if T_C >= T_H:
        err("  T_C must be less than T_H."); return
    eta   = 1 - T_C/T_H
    W_out = float(input("  Heat input Q_H (J) = "))
    Q_C   = W_out * T_C / T_H
    W     = W_out - Q_C
    ok(f"  Carnot efficiency η = {eta:.4f} ({eta*100:.2f}%)")
    ok(f"  Work output       W = {W:.4f} J")
    ok(f"  Heat rejected   Q_C = {Q_C:.4f} J")


# ============================================================
# 5.  STATISTICS
# ============================================================

def descriptive_stats():
    section("Descriptive Statistics")
    raw = input("  Enter data values (space-separated): ")
    data = list(map(float, raw.split()))
    n = len(data)
    if n == 0:
        err("  No data."); return

    mean = st.mean(data)
    med  = st.median(data)
    mode_val = None
    try: mode_val = st.mode(data)
    except: pass
    var  = st.variance(data)      # sample variance
    sd   = st.stdev(data)
    mn, mx = min(data), max(data)

    ok(f"  n        = {n}")
    ok(f"  Mean     = {mean:.6f}")
    ok(f"  Median   = {med:.6f}")
    if mode_val is not None: ok(f"  Mode     = {mode_val:.6f}")
    ok(f"  Std Dev  = {sd:.6f}")
    ok(f"  Variance = {var:.6f}")
    ok(f"  Min      = {mn:.6f}")
    ok(f"  Max      = {mx:.6f}")
    ok(f"  Range    = {mx - mn:.6f}")

    if SCI and n >= 2:
        # Histogram
        plt.figure(figsize=(7,4))
        plt.hist(data, bins="auto", color="#48cae4", edgecolor="#023e8a", alpha=0.8)
        plt.axvline(mean, color="#f72585", lw=2, label=f"Mean={mean:.3f}")
        plt.axvline(med,  color="#ffd60a", lw=2, ls="--", label=f"Median={med:.3f}")
        plt.title("Data Distribution"); plt.xlabel("Value"); plt.ylabel("Frequency")
        plt.legend(); plt.tight_layout(); plt.show()


def hypothesis_test():
    if not SCI:
        err("  Requires scipy."); return
    section("Hypothesis Testing")
    print("  1) One-sample t-test")
    print("  2) Two-sample t-test (independent)")
    print("  3) Chi-squared goodness-of-fit")
    ch = input("  Choice: ").strip()
    if ch == "1":
        raw = input("  Sample data (space-separated): ")
        data = list(map(float, raw.split()))
        mu0 = float(input("  Null hypothesis mean μ₀ = "))
        t_stat, p_val = stats.ttest_1samp(data, mu0)
        ok(f"  t-statistic = {t_stat:.4f}")
        ok(f"  p-value     = {p_val:.6f}")
        ok(f"  {'REJECT' if p_val < 0.05 else 'FAIL TO REJECT'} H₀ at α=0.05")
    elif ch == "2":
        r1 = input("  Group 1 data: ")
        r2 = input("  Group 2 data: ")
        d1 = list(map(float, r1.split()))
        d2 = list(map(float, r2.split()))
        t_stat, p_val = stats.ttest_ind(d1, d2)
        ok(f"  t-statistic = {t_stat:.4f}")
        ok(f"  p-value     = {p_val:.6f}")
        ok(f"  {'REJECT' if p_val < 0.05 else 'FAIL TO REJECT'} H₀ at α=0.05")
    elif ch == "3":
        obs_raw  = input("  Observed counts: ")
        exp_raw  = input("  Expected counts: ")
        observed = list(map(float, obs_raw.split()))
        expected = list(map(float, exp_raw.split()))
        chi2, p_val = stats.chisquare(observed, expected)
        ok(f"  χ² = {chi2:.4f}")
        ok(f"  p-value = {p_val:.6f}")
        ok(f"  {'REJECT' if p_val < 0.05 else 'FAIL TO REJECT'} H₀ at α=0.05")


def regression():
    if not SCI:
        err("  Requires scipy / numpy."); return
    section("Linear Regression  y = mx + b")
    raw_x = input("  x values (space-separated): ")
    raw_y = input("  y values (space-separated): ")
    x = np.array(list(map(float, raw_x.split())))
    y = np.array(list(map(float, raw_y.split())))
    if len(x) != len(y):
        err("  x and y must have the same length."); return
    slope, intercept, r, p_val, se = stats.linregress(x, y)
    ok(f"  Slope      m = {slope:.6f}")
    ok(f"  Intercept  b = {intercept:.6f}")
    ok(f"  R-value    r = {r:.6f}")
    ok(f"  R²           = {r**2:.6f}")
    ok(f"  Std error    = {se:.6f}")
    ok(f"  p-value      = {p_val:.6e}")

    x_fit = np.linspace(x.min(), x.max(), 200)
    y_fit = slope * x_fit + intercept
    plt.figure(figsize=(7,5))
    plt.scatter(x, y, color="#f72585", zorder=5, label="Data")
    plt.plot(x_fit, y_fit, color="#4cc9f0", lw=2, label=f"y={slope:.3f}x+{intercept:.3f}")
    plt.title("Linear Regression"); plt.xlabel("x"); plt.ylabel("y")
    plt.legend(); plt.grid(alpha=0.3); plt.tight_layout(); plt.show()


# ============================================================
# MAIN MENU
# ============================================================

MENU = [
    # (label, function)
    ("─── ALGEBRA ─────────────────────────────────", None),
    (" 1) Quadratic Equation Solver",               quadratic_solver),
    (" 2) Linear System Solver  (numpy linalg)",    linear_system_solver),
    (" 3) Polynomial Root Finder  (numpy)",         polynomial_roots),
    ("─── CALCULUS ─────────────────────────────────", None),
    (" 4) Numerical Derivative",                    numerical_derivative),
    (" 5) Numerical Integral  (Simpson + scipy)",   numerical_integral),
    (" 6) ODE Solver  (scipy.solve_ivp + plot)",    ode_solver),
    ("─── MECHANICS ────────────────────────────────", None),
    (" 7) 1-D Kinematics",                          kinematics_1d),
    (" 8) Projectile Motion  (+ plot)",             projectile_motion),
    (" 9) Energy & Work",                           energy_work),
    ("─── THERMODYNAMICS ───────────────────────────", None),
    ("10) Ideal Gas Law",                           ideal_gas_law),
    ("11) Heat Transfer  Q = mcΔT",                 heat_transfer),
    ("12) Carnot Cycle",                            carnot_efficiency),
    ("─── STATISTICS ───────────────────────────────", None),
    ("13) Descriptive Statistics  (+ histogram)",   descriptive_stats),
    ("14) Hypothesis Testing  (t-test, chi²)",      hypothesis_test),
    ("15) Linear Regression  (+ plot)",             regression),
]

ACTION_MAP = {
    str(i+1): fn for i, (_, fn)
    in enumerate([(l, f) for l, f in MENU if f is not None])
}

def main():
    banner = (
        f"\n{BOLD}{CYAN}"
        "╔══════════════════════════════════════════════════════╗\n"
        "║  University Math & Physics Toolkit  [Python Edition] ║\n"
        "╚══════════════════════════════════════════════════════╝"
        f"{RESET}"
    )
    while True:
        print(banner)
        k = 0
        for label, fn in MENU:
            if fn is None:
                print(f"\n{BOLD}{label}{RESET}")
            else:
                k += 1
                print(f"  {label}")
        print(f"\n{BOLD}  {'─'*44}{RESET}")
        print("   0) Exit")
        choice = input(f"\n{CYAN}  Your choice: {RESET}").strip()
        if choice == "0":
            print(CYAN + "\n  Goodbye!\n" + RESET)
            sys.exit(0)
        fn = ACTION_MAP.get(choice)
        if fn:
            try:
                fn()
            except ValueError as e:
                err(f"  Input error: {e}")
            except KeyboardInterrupt:
                print()
        else:
            err("  Unknown option.")
        input(GREY + "\n  [Press Enter to continue]" + RESET)

if __name__ == "__main__":
    main()
