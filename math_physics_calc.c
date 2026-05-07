/*
 * ============================================================
 *  University Math & Physics Calculator  —  C Project
 *  Topics: Algebra, Calculus, Mechanics, Thermodynamics,
 *          Statistics
 * ============================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846
#define G  9.80665   /* standard gravity m/s^2 */
#define R  8.314     /* universal gas constant J/(mol·K) */

/* ─── ANSI colours ─────────────────────────────────────────── */
#define BOLD   "\033[1m"
#define CYAN   "\033[36m"
#define GREEN  "\033[32m"
#define YELLOW "\033[33m"
#define RED    "\033[31m"
#define RESET  "\033[0m"

/* ============================================================
   1. ALGEBRA
   ============================================================ */

/* Solve ax^2 + bx + c = 0 */
void quadratic_solver(void) {
    double a, b, c;
    printf(CYAN "\n── Quadratic Solver: ax² + bx + c = 0 ──\n" RESET);
    printf("  a = "); scanf("%lf", &a);
    printf("  b = "); scanf("%lf", &b);
    printf("  c = "); scanf("%lf", &c);

    if (fabs(a) < 1e-12) {
        if (fabs(b) < 1e-12) {
            printf(RED "  Degenerate equation (a=0, b=0).\n" RESET);
        } else {
            printf(GREEN "  Linear root: x = %.6f\n" RESET, -c/b);
        }
        return;
    }

    double disc = b*b - 4*a*c;
    if (disc > 0) {
        double sq = sqrt(disc);
        printf(GREEN "  Two real roots:\n"
               "    x₁ = %.6f\n    x₂ = %.6f\n" RESET,
               (-b + sq)/(2*a), (-b - sq)/(2*a));
    } else if (fabs(disc) < 1e-12) {
        printf(GREEN "  One repeated root: x = %.6f\n" RESET, -b/(2*a));
    } else {
        double re = -b / (2*a);
        double im = sqrt(-disc) / (2*a);
        printf(YELLOW "  Complex roots:\n"
               "    x₁ = %.6f + %.6fi\n"
               "    x₂ = %.6f - %.6fi\n" RESET,
               re, im, re, im);
    }
}

/* Gaussian elimination for Ax = b (up to 4×4) */
void linear_system_solver(void) {
    int n;
    printf(CYAN "\n── Linear System Solver (Gaussian Elimination) ──\n" RESET);
    printf("  Number of equations (1-4): ");
    scanf("%d", &n);
    if (n < 1 || n > 4) { printf(RED "  Out of range.\n" RESET); return; }

    double a[4][5];
    printf("  Enter augmented matrix [A|b] row by row:\n");
    for (int i = 0; i < n; i++) {
        printf("  Row %d: ", i+1);
        for (int j = 0; j <= n; j++) scanf("%lf", &a[i][j]);
    }

    /* Forward elimination */
    for (int col = 0; col < n; col++) {
        /* Partial pivot */
        int pivot = col;
        for (int row = col+1; row < n; row++)
            if (fabs(a[row][col]) > fabs(a[pivot][col])) pivot = row;
        if (pivot != col)
            for (int j = 0; j <= n; j++) {
                double tmp = a[col][j]; a[col][j] = a[pivot][j]; a[pivot][j] = tmp;
            }
        if (fabs(a[col][col]) < 1e-12) { printf(RED "  Singular matrix.\n" RESET); return; }
        for (int row = col+1; row < n; row++) {
            double factor = a[row][col] / a[col][col];
            for (int j = col; j <= n; j++) a[row][j] -= factor * a[col][j];
        }
    }

    /* Back substitution */
    double x[4];
    for (int i = n-1; i >= 0; i--) {
        x[i] = a[i][n];
        for (int j = i+1; j < n; j++) x[i] -= a[i][j] * x[j];
        x[i] /= a[i][i];
    }

    printf(GREEN "\n  Solution:\n" RESET);
    for (int i = 0; i < n; i++) printf(GREEN "    x%d = %.6f\n" RESET, i+1, x[i]);
}

/* ============================================================
   2. CALCULUS
   ============================================================ */

/* Numerical derivative (central difference) of user-defined f */
typedef double (*Func)(double);

/* We store a simple expression via a menu */
static int g_func_choice = 0;

double eval_func(double x) {
    switch (g_func_choice) {
        case 1: return sin(x);
        case 2: return cos(x);
        case 3: return exp(x);
        case 4: return log(fabs(x));
        case 5: return x*x*x - 3*x*x + 2*x - 1;
        case 6: return sin(x)*exp(-x);
        default: return x*x;
    }
}

const char *func_name(void) {
    switch (g_func_choice) {
        case 1: return "sin(x)";
        case 2: return "cos(x)";
        case 3: return "e^x";
        case 4: return "ln|x|";
        case 5: return "x³ - 3x² + 2x - 1";
        case 6: return "sin(x)·e^(-x)";
        default: return "x²";
    }
}

void pick_function(void) {
    printf("  Select function:\n"
           "    1) sin(x)          2) cos(x)\n"
           "    3) e^x             4) ln|x|\n"
           "    5) x³-3x²+2x-1    6) sin(x)e^(-x)\n"
           "    0) x²  (default)\n"
           "  Choice: ");
    scanf("%d", &g_func_choice);
}

void numerical_derivative(void) {
    printf(CYAN "\n── Numerical Derivative (Central Difference) ──\n" RESET);
    pick_function();
    double x, h = 1e-5;
    printf("  Evaluate f'(x) at x = "); scanf("%lf", &x);
    double deriv = (eval_func(x+h) - eval_func(x-h)) / (2*h);
    printf(GREEN "  f(x)  = %s\n  f'(%.4f) ≈ %.8f\n" RESET,
           func_name(), x, deriv);
}

/* Simpson's 1/3 rule */
void numerical_integral(void) {
    printf(CYAN "\n── Numerical Integration (Simpson's 1/3 Rule) ──\n" RESET);
    pick_function();
    double a, b;
    int n;
    printf("  Lower limit a = "); scanf("%lf", &a);
    printf("  Upper limit b = "); scanf("%lf", &b);
    printf("  Intervals n (even, e.g. 1000): "); scanf("%d", &n);
    if (n % 2 != 0) n++;

    double h = (b - a) / n;
    double sum = eval_func(a) + eval_func(b);
    for (int i = 1; i < n; i++)
        sum += (i % 2 == 0 ? 2 : 4) * eval_func(a + i*h);
    double result = sum * h / 3.0;
    printf(GREEN "  ∫[%.4f→%.4f] %s dx ≈ %.8f\n" RESET,
           a, b, func_name(), result);
}

/* Newton-Raphson root finding */
void newton_raphson(void) {
    printf(CYAN "\n── Newton-Raphson Root Finder ──\n" RESET);
    pick_function();
    double x;
    printf("  Initial guess x₀ = "); scanf("%lf", &x);
    double h = 1e-6;
    int iter;
    for (iter = 0; iter < 1000; iter++) {
        double fx = eval_func(x);
        double fpx = (eval_func(x+h) - eval_func(x-h)) / (2*h);
        if (fabs(fpx) < 1e-14) { printf(RED "  Zero derivative — method failed.\n" RESET); return; }
        double x1 = x - fx / fpx;
        if (fabs(x1 - x) < 1e-10) { x = x1; break; }
        x = x1;
    }
    printf(GREEN "  Root of %s: x ≈ %.10f  (after %d iterations)\n" RESET,
           func_name(), x, iter);
}

/* ============================================================
   3. MECHANICS
   ============================================================ */

void kinematics_1d(void) {
    printf(CYAN "\n── 1-D Kinematics ──\n" RESET);
    double u, a_acc, t;
    printf("  Initial velocity u (m/s)  = "); scanf("%lf", &u);
    printf("  Acceleration     a (m/s²) = "); scanf("%lf", &a_acc);
    printf("  Time             t (s)    = "); scanf("%lf", &t);

    double v  = u + a_acc * t;
    double s  = u*t + 0.5*a_acc*t*t;
    double v2 = u*u + 2*a_acc*s;

    printf(GREEN
        "  Results:\n"
        "    Final velocity  v  = %.4f m/s\n"
        "    Displacement    s  = %.4f m\n"
        "    v² check        v² = %.4f m²/s²\n"
        RESET, v, s, v2);
}

void projectile_motion(void) {
    printf(CYAN "\n── Projectile Motion ──\n" RESET);
    double v0, theta_deg, h0;
    printf("  Launch speed   v₀  (m/s)  = "); scanf("%lf", &v0);
    printf("  Launch angle   θ   (deg)  = "); scanf("%lf", &theta_deg);
    printf("  Initial height h₀  (m)    = "); scanf("%lf", &h0);

    double theta = theta_deg * PI / 180.0;
    double vx = v0 * cos(theta);
    double vy = v0 * sin(theta);

    /* Time of flight: h0 + vy*t - 0.5*G*t^2 = 0 */
    double disc = vy*vy + 2*G*h0;
    if (disc < 0) { printf(RED "  No real solution.\n" RESET); return; }
    double t_flight = (vy + sqrt(disc)) / G;

    double range    = vx * t_flight;
    double t_peak   = vy / G;
    double h_peak   = (t_peak >= 0) ? h0 + vy*t_peak - 0.5*G*t_peak*t_peak : h0;

    printf(GREEN
        "  Results:\n"
        "    Time of flight  = %.4f s\n"
        "    Range           = %.4f m\n"
        "    Peak height     = %.4f m\n"
        "    vₓ              = %.4f m/s\n"
        "    v_y0            = %.4f m/s\n"
        RESET, t_flight, range, h_peak, vx, vy);
}

void circular_motion(void) {
    printf(CYAN "\n── Uniform Circular Motion ──\n" RESET);
    double r, v;
    printf("  Radius r (m)    = "); scanf("%lf", &r);
    printf("  Speed  v (m/s)  = "); scanf("%lf", &v);

    double T    = 2*PI*r / v;
    double omega = v / r;
    double a_c  = v*v / r;
    printf(GREEN
        "  Results:\n"
        "    Period          T  = %.4f s\n"
        "    Angular vel.    ω  = %.4f rad/s\n"
        "    Centripetal acc a  = %.4f m/s²\n"
        RESET, T, omega, a_c);
}

/* ============================================================
   4. THERMODYNAMICS
   ============================================================ */

void ideal_gas_law(void) {
    printf(CYAN "\n── Ideal Gas Law  PV = nRT ──\n" RESET);
    printf("  Which variable to find? [P/V/n/T]: ");
    char var[4]; scanf("%s", var);

    if (strcasecmp(var,"P")==0) {
        double n,V,T; printf("  n (mol)="); scanf("%lf",&n);
        printf("  V (L)="); scanf("%lf",&V);
        printf("  T (K)="); scanf("%lf",&T);
        printf(GREEN "  P = %.4f atm\n" RESET, n*0.082057*T/V);
    } else if (strcasecmp(var,"V")==0) {
        double n,P,T; printf("  n (mol)="); scanf("%lf",&n);
        printf("  P (atm)="); scanf("%lf",&P);
        printf("  T (K)="); scanf("%lf",&T);
        printf(GREEN "  V = %.4f L\n" RESET, n*0.082057*T/P);
    } else if (strcasecmp(var,"n")==0) {
        double P,V,T; printf("  P (atm)="); scanf("%lf",&P);
        printf("  V (L)="); scanf("%lf",&V);
        printf("  T (K)="); scanf("%lf",&T);
        printf(GREEN "  n = %.4f mol\n" RESET, P*V/(0.082057*T));
    } else {
        double P,V,n; printf("  P (atm)="); scanf("%lf",&P);
        printf("  V (L)="); scanf("%lf",&V);
        printf("  n (mol)="); scanf("%lf",&n);
        printf(GREEN "  T = %.4f K\n" RESET, P*V/(n*0.082057));
    }
}

void heat_transfer(void) {
    printf(CYAN "\n── Sensible Heat Transfer  Q = mcΔT ──\n" RESET);
    double m, c, dT;
    printf("  Mass m (kg)                     = "); scanf("%lf", &m);
    printf("  Specific heat c (J/(kg·K))      = "); scanf("%lf", &c);
    printf("  Temperature change ΔT (K or °C) = "); scanf("%lf", &dT);
    printf(GREEN "  Q = %.4f J  (%.4f kJ)\n" RESET, m*c*dT, m*c*dT/1000);
}

void carnot_efficiency(void) {
    printf(CYAN "\n── Carnot Efficiency  η = 1 - T_cold/T_hot ──\n" RESET);
    double Th, Tc;
    printf("  Hot reservoir  T_H (K) = "); scanf("%lf", &Th);
    printf("  Cold reservoir T_C (K) = "); scanf("%lf", &Tc);
    if (Tc >= Th) { printf(RED "  T_C must be < T_H.\n" RESET); return; }
    double eta = 1.0 - Tc/Th;
    printf(GREEN "  Carnot efficiency η = %.4f (%.2f%%)\n" RESET, eta, eta*100);
}

/* ============================================================
   5. STATISTICS
   ============================================================ */

void statistics_calc(void) {
    printf(CYAN "\n── Descriptive Statistics ──\n" RESET);
    int n;
    printf("  How many data points? "); scanf("%d", &n);
    if (n < 1) return;
    double *data = malloc(n * sizeof(double));
    if (!data) { printf(RED "  Memory error.\n" RESET); return; }
    printf("  Enter values:\n");
    for (int i = 0; i < n; i++) { printf("  [%d] ", i+1); scanf("%lf", &data[i]); }

    /* Mean */
    double sum = 0;
    for (int i = 0; i < n; i++) sum += data[i];
    double mean = sum / n;

    /* Variance & std dev */
    double var_sum = 0;
    for (int i = 0; i < n; i++) var_sum += (data[i]-mean)*(data[i]-mean);
    double variance = var_sum / (n-1);   /* sample variance */
    double stddev   = sqrt(variance);

    /* Min / Max */
    double mn = data[0], mx = data[0];
    for (int i = 1; i < n; i++) {
        if (data[i] < mn) mn = data[i];
        if (data[i] > mx) mx = data[i];
    }

    /* Median — sort copy */
    double *sorted = malloc(n * sizeof(double));
    memcpy(sorted, data, n * sizeof(double));
    for (int i = 0; i < n-1; i++)
        for (int j = i+1; j < n; j++)
            if (sorted[j] < sorted[i]) { double t=sorted[i]; sorted[i]=sorted[j]; sorted[j]=t; }
    double median = (n%2==0) ? (sorted[n/2-1]+sorted[n/2])/2.0 : sorted[n/2];

    printf(GREEN
        "  Results:\n"
        "    Count    n  = %d\n"
        "    Mean     μ  = %.6f\n"
        "    Median      = %.6f\n"
        "    Std Dev  s  = %.6f\n"
        "    Variance s² = %.6f\n"
        "    Min         = %.6f\n"
        "    Max         = %.6f\n"
        "    Range       = %.6f\n"
        RESET, n, mean, median, stddev, variance, mn, mx, mx-mn);

    free(data); free(sorted);
}

void normal_distribution(void) {
    printf(CYAN "\n── Normal Distribution: P(X ≤ x) ──\n" RESET);
    double mu, sigma, x;
    printf("  Mean  μ = "); scanf("%lf", &mu);
    printf("  Std   σ = "); scanf("%lf", &sigma);
    printf("  Value x = "); scanf("%lf", &x);
    double z = (x - mu) / sigma;
    /* erfc approximation */
    double p = 0.5 * erfc(-z / sqrt(2.0));
    printf(GREEN
        "  z-score = %.4f\n"
        "  P(X ≤ %.4f) = %.6f (%.2f%%)\n"
        RESET, z, x, p, p*100);
}

/* ============================================================
   MAIN MENU
   ============================================================ */

void print_banner(void) {
    printf(BOLD CYAN
    "\n╔══════════════════════════════════════════════════════╗\n"
    "║   University Math & Physics Calculator  [C Edition]  ║\n"
    "╚══════════════════════════════════════════════════════╝\n" RESET);
}

void print_menu(void) {
    printf(BOLD "\n  ── ALGEBRA ─────────────────────────────────\n" RESET);
    printf("   1) Quadratic Equation Solver\n");
    printf("   2) Linear System Solver (Gaussian Elim.)\n");
    printf(BOLD "  ── CALCULUS ─────────────────────────────────\n" RESET);
    printf("   3) Numerical Derivative\n");
    printf("   4) Numerical Integral (Simpson)\n");
    printf("   5) Root Finder (Newton-Raphson)\n");
    printf(BOLD "  ── MECHANICS ────────────────────────────────\n" RESET);
    printf("   6) 1-D Kinematics\n");
    printf("   7) Projectile Motion\n");
    printf("   8) Uniform Circular Motion\n");
    printf(BOLD "  ── THERMODYNAMICS ───────────────────────────\n" RESET);
    printf("   9) Ideal Gas Law\n");
    printf("  10) Heat Transfer (Q = mcΔT)\n");
    printf("  11) Carnot Efficiency\n");
    printf(BOLD "  ── STATISTICS ───────────────────────────────\n" RESET);
    printf("  12) Descriptive Statistics\n");
    printf("  13) Normal Distribution CDF\n");
    printf(BOLD "  ─────────────────────────────────────────────\n" RESET);
    printf("   0) Exit\n\n");
    printf("  Your choice: ");
}

int main(void) {
    int choice;
    print_banner();
    do {
        print_menu();
        scanf("%d", &choice);
        switch (choice) {
            case  1: quadratic_solver();    break;
            case  2: linear_system_solver();break;
            case  3: numerical_derivative();break;
            case  4: numerical_integral();  break;
            case  5: newton_raphson();      break;
            case  6: kinematics_1d();       break;
            case  7: projectile_motion();   break;
            case  8: circular_motion();     break;
            case  9: ideal_gas_law();       break;
            case 10: heat_transfer();       break;
            case 11: carnot_efficiency();   break;
            case 12: statistics_calc();     break;
            case 13: normal_distribution(); break;
            case  0: printf(CYAN "\n  Goodbye!\n\n" RESET); break;
            default: printf(RED "  Unknown option.\n" RESET);
        }
    } while (choice != 0);
    return 0;
}
