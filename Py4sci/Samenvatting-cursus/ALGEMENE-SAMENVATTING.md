# ALGEMENE-SAMENVATTING

Deze algemene samenvatting bevat de samenvattingen van Hoofdstuk 1 t/m Hoofdstuk 12, samengevoegd uit de per-hoofdstuk bestanden.

---

## Hoofdstuk 1 — Numerical Limitations (waarom numeriek rekenen soms “vals speelt”)

Dit hoofdstuk is de *grondwet* van de hele cursus: elke numerieke methode is een trade-off tussen **benaderingsfout** (truncation) en **rekenfout** (rounding). Als je op het examen moet uitleggen *waarom methode A beter is dan B*, zit het argument bijna altijd hier: **stabiliteit + conditioning + foutopbouw + kost**.

---

## Hoofdstuk 2 — Systems of Linear Equations

Dit hoofdstuk leert je **lineaire stelsels** oplossen
$$
\mathbf{A}\mathbf{x}=\mathbf{b},
$$
en vooral: **hoe je dat doet op een manier die (i) snel is en (ii) numeriek stabiel** is in floating-point (link met H1).

De rode draad: transformeer het probleem naar een vorm die je goedkoop kan oplossen (triangulair), en doe dat met controle op stabiliteit (pivoting). Belangrijke methodes: Gaussian elimination, LU, Cholesky (voor SPD), en pivoting/partial pivoting voor stabiliteit.

---

## Hoofdstuk 3 — Linear Least Squares

In dit hoofdstuk los je problemen op van het type overdetermined ($m>n$) waarbij $\mathbf{A}\mathbf{x}=\mathbf{b}$ geen exacte oplossing heeft. Je zoekt $\hat{\mathbf{x}}$ die $\min \|\mathbf{b}-\mathbf{A}\mathbf{x}\|_2$.

Belangrijke keuzes: normal equations (goedkoop, maar condition number kwadrateert), QR (Householder, stabiel) en SVD (duurste, maar gouden standaard bij rank/diagnose).

---

## Hoofdstuk 4 — Eigenvalue Problems

Eigenproblemen $\mathbf{A}\mathbf{x}=\lambda\mathbf{x}$: methodes zijn o.a. power iteration (dominante mode), inverse iteration (kleinste of dicht bij een shift), Rayleigh quotient iteration (snelle lokale convergentie), deflation (voor meerdere waarden), en QR-iteratie (algemene solver voor alle eigenwaarden). Kies methodes op basis van doel (1 waarde vs allemaal), kosten en stabiliteit.

---

## Hoofdstuk 5 — Nonlinear equations

Oplossen van $f(x)=0$ of $\mathbf{f}(\mathbf{x})=\mathbf{0}$. Methoden: bisection/Brent (robust), secant, Newton (kwadratisch dicht bij root), inverse interpolation, en voor systemen Newton met Jacobiaan (lineaire solves per stap). Belangrijk: conditioning van de wortel en passende stopcriteria (residual/step/tolerance).

---

## Hoofdstuk 6 — Optimization

Optimalisatiemethodes: golden section / Brent (1D, geen gradient), gradient descent, Newton, quasi-Newton (BFGS/L-BFGS), trust-region, Nelder–Mead (direct search). Voor constrained problems gebruik je Lagrange/KKT-ideeën of specifieke solvers (L-BFGS-B, SLSQP, trust-constr). Nonlinear least squares: Gauss–Newton, Levenberg–Marquardt.

---

## Hoofdstuk 7 — Interpolation

Interpolatie: Lagrange, Newton, Vandermonde (gevaarlijk), barycentric (stabiel), splines (piecewise cubic), PCHIP (shape-preserving). Kies global vs piecewise afhankelijk van aantal punten en ruis.

---

## Hoofdstuk 8 — Numerical Integration and Differentiation

Integratie: composite Newton–Cotes, Gaussian quadrature, adaptive quadrature. Differentiatie: finite differences (forward/central), Richardson extrapolation en trade-off truncation vs rounding — er bestaat een optimale $h$.

---

## Hoofdstuk 9 — Ordinary Differential Equations (ODEs)

IVP methodes: explicit RK (RK4, adaptive RK45), implicit methods (Backward Euler, BDF, Radau) voor stiff problems, adaptive step size, and BVP methods (shooting, collocation). Kies op basis van stiffness en kosten.

---

## Hoofdstuk 10 — Partial Differential Equations (PDEs)

PDE numeriek: discretiseer ruimte (FD/FE/spectral), kies tijdstapper (expliciet/impliciet), controleer stabiliteit (CFL/stiffness), en los per stap lineaire/nonlineaire systemen op (sparse direct or iterative solvers zoals CG/GMRES).

---

## Hoofdstuk 11 — Fast Fourier Transform (FFT)

FFT: snelle DFT ($O(N \log N)$), aandachtspunten: spectral leakage (windowing), zero-padding voor lineaire convolutie, aliasing (Nyquist), DFT vs DCT (reflective extension). Use `rfft` for real signals and `next_fast_len` for performance.

---

## Hoofdstuk 12 — Monte Carlo

Monte Carlo: sampling-based integratie met error $\sim 1/\sqrt{N}$, variance reduction (importance sampling, stratification, control variates, antithetic variates), and MCMC (Metropolis–Hastings, burn-in, mixing).

---

De volledige hoofdstukken en details vind je in de individuele bestanden `Hoofdstuk-01/README.md` t/m `Hoofdstuk-12/README.md`.

