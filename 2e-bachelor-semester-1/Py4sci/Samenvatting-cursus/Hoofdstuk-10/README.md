<!-- file: hoofdstuk_10_partial_differential_equations.md -->
# Hoofdstuk 10 — Partial Differential Equations (PDEs)

PDE’s beschrijven velden $u(\mathbf{x},t)$ met ruimte én tijd. Numeriek komt het altijd neer op:
1) discretiseer ruimte (finite differences / finite elements / spectral),
2) kies een tijdstapper (expliciet of impliciet),
3) controleer stabiliteit (CFL / stiffness),
4) los per stap een lineair of nonlineair stelsel op.

De examenkern in Py4Sci-stijl: **mapping PDE → discretisatie → type stelsel → geschikte solver**, plus **kost vs stabiliteit**.

---

## 10.1 Classificatie (waarom je dit meteen wil weten)

Voor 1D voorbeelden:

### Elliptisch (stationair)
Laplace/Poisson:
$$
u_{xx} = f(x)
\quad\text{of}\quad
\nabla^2 u = f.
$$
- Geen tijd.
- Discretisatie $\Rightarrow$ **lineair stelsel** $\mathbf{A}\mathbf{u}=\mathbf{b}$ (vaak SPD en sparse).

### Parabolisch (diffusie)
Heat equation:
$$
u_t = \kappa u_{xx}.
$$
- Eén tijdafgeleide, tweede orde in ruimte.
- Expliciete methodes hebben **strenge stabiliteitslimiet**: $\Delta t = O((\Delta x)^2)$.

### Hyperbolisch (golven/transport)
Wave equation:
$$
u_{tt} = c^2 u_{xx}
$$
of advection:
$$
u_t + c u_x = 0.
$$
- Propagatie; stabiliteit bepaalt meestal $\Delta t = O(\Delta x)$ (CFL).

**Waarom dit examenrelevant is:** het vertelt je meteen of een probleem “stiff” wordt (parabolisch vaak wel) en dus impliciete methodes verdient.

---

## 10.2 Ruimtediscretisatie (finite differences als default)

Rooster in 1D:
$$
x_i = a + i\Delta x,\qquad i=0,\dots,N.
$$

### 10.2.1 Centrale stencil’s (zoals in integratie/differentiatie hoofdstuk)
Eerste afgeleide:
$$
u_x(x_i)\approx \frac{u_{i+1}-u_{i-1}}{2\Delta x}.
$$
Tweede afgeleide:
$$
u_{xx}(x_i)\approx \frac{u_{i+1}-2u_i+u_{i-1}}{(\Delta x)^2}.
$$

### 10.2.2 Boundary conditions (BCs) als “extra vergelijkingen”
- Dirichlet: $u(a)=\alpha$, $u(b)=\beta$ (fixeer randwaarden).
- Neumann: $u_x(a)=\gamma$ (randstencil of ghost points).
- Robin: combinatie.

BC’s bepalen hoe je de eerste/laatste rijen van je matrix $\mathbf{A}$ opbouwt.

---

## 10.3 Stationaire PDE’s: Poisson/Laplace $\Rightarrow$ sparse linear system

Voor Poisson in 1D:
$$
-u_{xx}=f(x),
$$
met centrale stencil krijg je:
$$
\frac{-u_{i+1}+2u_i-u_{i-1}}{(\Delta x)^2}=f_i.
$$

Dit levert een tridiagonale (sparse) matrix $\mathbf{A}$:
$$
\mathbf{A}\mathbf{u}=\mathbf{b}.
$$

### 10.3.1 Welke solver?
- Kleine $N$ (dense of klein sparse): directe solve.
- Grote $N$ (sparse): `spsolve` of iteratief:
  - SPD: Conjugate Gradient (CG),
  - niet-SPD: GMRES/BiCGSTAB (afhankelijk van structuur).

**Waarom:** sparse opslag + sparse matrix-vector products zijn $O(N)$ i.p.v. $O(N^2)$.

### 10.3.2 Conditioning en grid
Conditioning van discretisatie-operators wordt slechter als $\Delta x$ kleiner wordt (meer roosterpunten). Dit beïnvloedt:
- aantal iteraties bij iteratieve solvers,
- gevoeligheid voor rounding.

---

## 10.4 Time-dependent PDE’s: method of lines

Discretiseer eerst ruimte. Dan krijg je een ODE-systeem:
$$
\mathbf{u}'(t)=\mathbf{F}(t,\mathbf{u}(t)).
$$
Daarna gebruik je ODE-methodes (H9) op $\mathbf{u}(t)$.

**Dit is een sleutelbrug tussen PDE en ODE**: PDE numeriek oplossen = “grote ODE” + (sparse) lineaire algebra.

---

## 10.5 Tijdstepping: expliciet vs impliciet

### 10.5.1 Heat equation (parabolisch): waarom expliciet vaak faalt
Neem FTCS:
$$
u_i^{n+1}=u_i^n + \mu\left(u_{i+1}^n-2u_i^n+u_{i-1}^n\right),
\qquad
\mu=\kappa\frac{\Delta t}{(\Delta x)^2}.
$$

Stabiliteit vereist typisch:
$$
\mu \le \frac{1}{2}\quad \text{(1D klassieke bound)}.
$$

**Interpretatie:** bij kleinere $\Delta x$ moet $\Delta t$ kwadratisch kleiner → veel stappen → duur.

#### Impliciet (Backward Euler)
$$
u_i^{n+1}=u_i^n + \mu\left(u_{i+1}^{n+1}-2u_i^{n+1}+u_{i-1}^{n+1}\right).
$$
In matrixvorm:
$$
(\mathbf{I}-\mu\mathbf{L})\mathbf{u}^{n+1}=\mathbf{u}^n.
$$

- Pro: veel stabieler → grotere $\Delta t$ mogelijk.
- Con: per stap een sparse solve (maar dat is vaak goedkoper dan miljoenen expliciete stappen).

#### Crank–Nicolson (CN)
Gemiddelde van expliciet/impliciet:
$$
(\mathbf{I}-\tfrac{\mu}{2}\mathbf{L})\mathbf{u}^{n+1}
=
(\mathbf{I}+\tfrac{\mu}{2}\mathbf{L})\mathbf{u}^{n}.
$$
- Vaak 2e orde in tijd en stabieler dan expliciet.
- Kan oscillaties geven bij ruwe data als je te grote $\Delta t$ neemt (numeriek “ringing”).

---

### 10.5.2 Hyperbolisch (advection/waves): CFL is de baas
Voor advection $u_t+cu_x=0$ is een typische CFL-conditie:
$$
\frac{|c|\Delta t}{\Delta x}\le C_{\max}.
$$

**Fysische betekenis:** informatie beweegt met snelheid $c$; per tijdstap mag een “feature” niet meer dan ongeveer één cel opschuiven, anders mis je de dynamica en wordt het instabiel.

**Schema-keuze is cruciaal:**
- upwind schema’s zijn stabieler (maar numeriek diffuus),
- centrale schema’s kunnen oscilleren zonder extra stabilisatie.

---

## 10.6 Welke methode wanneer? (examengerichte beslisregels)

### Stationair elliptisch (Poisson/Laplace)
- Bouw sparse $\mathbf{A}$, solve $\mathbf{A}\mathbf{u}=\mathbf{b}$.
- SPD? CG (eventueel met preconditioner) of directe factorisatie als klein genoeg.

### Diffusie (heat)
- Kleine problemen: expliciet kan, maar check $\Delta t \sim (\Delta x)^2$.
- Grotere of fijne grids: impliciet/CN wint bijna altijd.

### Advection/waves
- Expliciet met CFL is vaak oké en efficiënt.
- Let op numerieke dispersie/oscillaties → kies passend schema (upwind/flux limiters in geavanceerde setting).

---

## 10.7 Typische valkuilen (die in “insights” gevraagd worden)

1) **Je kiest $\Delta t$ op nauwkeurigheid, maar stabiliteit is strenger** (diffusie).
2) **Circulaire artefacten** als je periodieke BC per ongeluk “meeneemt”.
3) **Sparse vs dense**: een PDE-matrix dense behandelen is dood door geheugen/kost.
4) **BC-implementatie**: vaak de bron van bugs (eerste/laatste rij fout).

---

## SciPy mapping (praktisch)
- `scipy.sparse` voor operators
- `scipy.sparse.linalg.spsolve`, `cg`, `gmres`, ...
- method of lines: gebruik `solve_ivp` met RHS die sparse operator toepast

---

## Links met andere hoofdstukken
- H8: finite differences + truncation/rounding.
- H9: method of lines = groot ODE-probleem.
- H2: impliciete stappen = lineaire solves (sparse LU/CG).
- H4: eigenwaarden van discretisatie-operator sturen stabiliteit (stiffness).
