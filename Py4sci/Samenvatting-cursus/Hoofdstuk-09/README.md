<!-- file: hoofdstuk_09_ordinary_differential_equations.md -->
# Hoofdstuk 9 — Ordinary Differential Equations (ODEs)

We lossen een initial value problem (IVP):
$$
\mathbf{y}'(t)=\mathbf{f}(t,\mathbf{y}(t)),
\qquad
\mathbf{y}(t_0)=\mathbf{y}_0,
$$
waar $\mathbf{y}(t)$ een vector is.

Examenkern:
- discretisatie + local/global truncation error,
- stiff vs non-stiff,
- expliciet vs impliciet,
- adaptieve step size,
- link met H5/H2 (impliciet → nonlineair + lineaire solves).

---

## 9.1 Introduction and useful concepts

- IVP: startwaarde gegeven, evolutie vooruit.
- BVP: randvoorwaarden op twee punten (zie 9.3).

Numerieke methode maakt een rooster:
$$
t_n=t_0+n h,
$$
en benadert $\mathbf{y}(t_n)\approx \mathbf{y}_n$.

---

## 9.2 Numerically solving ODE’s

### Discrete variable methods
We vervangen de differentiaalvergelijking door update-regels:
$$
\mathbf{y}_{n+1} = \mathbf{y}_n + h\,\Phi(t_n,\mathbf{y}_n,h),
$$
waar $\Phi$ de methode-specifieke “slope-combinatie” is.

---

### 9.2.1 Euler methods

**Forward Euler (expliciet):**
$$
\mathbf{y}_{n+1}=\mathbf{y}_n + h\,\mathbf{f}(t_n,\mathbf{y}_n).
$$
- Pro: goedkoop, simpel
- Con: lage orde, kan instabiel zijn

**Backward Euler (impliciet):**
$$
\mathbf{y}_{n+1}=\mathbf{y}_n + h\,\mathbf{f}(t_{n+1},\mathbf{y}_{n+1}).
$$
- Pro: veel stabieler, goed voor stiff problemen
- Con: vereist per stap het oplossen van een nonlineair probleem (H5); vaak Newton + lineaire solves (H2)

---

### 9.2.2 Runge–Kutta (RK) methods (expliciet, hogere orde)

Algemeen RK gebruikt meerdere “stages” $\mathbf{k}_i$ en combineert die.

**RK4 (klassieker):** hoge nauwkeurigheid per stap, maar vaste stapgrootte (tenzij je error-estimates toevoegt).

---

### 9.2.3 Adaptive step size (embedded RK)
In de praktijk wil je controle op fout:
- compute twee oplossingen van verschillende orde in dezelfde stap,
- schat error,
- pas $h$ aan.

**Waarom dit vaak wint:**
- je investeert rekenwerk waar het nodig is (snelle variatie),
- je spaart werk waar het niet nodig is.

---

### 9.2.4 Stability en stiffness (het echte “wanneer welke solver” criterium)

Een probleem is **stiff** als expliciete methodes een extreem klein $h$ moeten nemen voor stabiliteit, zelfs als de oplossing zelf niet zo snel varieert.

- **Non-stiff:** expliciete RK-methodes zijn vaak ideaal.
- **Stiff:** impliciete methodes (Backward Euler, BDF, Radau) winnen gigantisch.

**Kernboodschap:** stapgrootte wordt bij stiff problemen bepaald door stabiliteit, niet door nauwkeurigheid.

---

### 9.2.5 Impliciete methodes in de praktijk (link met H5/H2)
Een impliciete stap vraagt:
$$
\mathbf{y}_{n+1}-\mathbf{y}_n - h\,\mathbf{f}(t_{n+1},\mathbf{y}_{n+1})=\mathbf{0}.
$$
Dit is een nonlineair systeem in $\mathbf{y}_{n+1}$.
Newton geeft:
$$
\mathbf{J}\,\Delta \mathbf{y}=-\mathbf{g}(\mathbf{y}),
$$
waar $\mathbf{J}$ een Jacobiaan is en je lineaire solves nodig hebt (H2).

---

## 9.3 Boundary Value Problems (BVP’s) for ODE’s

BVP: je kent voorwaarden aan twee uiteinden, bv.
$$
y(a)=\alpha,\qquad y(b)=\beta.
$$

### (A) Shooting method
Maak van BVP een IVP met onbekende start-slope $s$ en zoek $s$ zodat de eindvoorwaarde klopt.
Dat is root finding (H5):
$$
F(s)=y(b;s)-\beta=0.
$$

- Pro: conceptueel simpel
- Con: kan instabiel zijn bij gevoelige problemen (slechte conditioning)

### (B) Finite difference / collocation
Discretiseer en los een (groot) algebraïsch systeem op, vaak sparse.
- Pro: robuuster voor lastige BVP’s
- Con: implementatie zwaarder, maar SciPy heeft tooling (`solve_bvp`)

---

## SciPy (typische mapping)

- `scipy.integrate.solve_ivp`
  - non-stiff: `RK45` / `DOP853`
  - stiff: `Radau` / `BDF`
- `scipy.integrate.solve_bvp` voor BVP

---

## Welke methode wanneer?

- **Snelle, gladde dynamica, niet stiff:** expliciete RK + adaptieve stap.
- **Stiff (veel tijdschalen):** impliciet (`BDF`/`Radau`).
- **BVP met eenvoudige fysica:** shooting + root finding (H5).
- **BVP moeilijk/sensitief:** finite difference/collocation (`solve_bvp`).

---

## Links met andere hoofdstukken
- H5: impliciete stappen en shooting vragen root finding.
- H2: Newton-stappen vragen lineaire solves (LU/sparse).
- H10: PDE time-stepping reduceert vaak tot ODE-systemen (method of lines).
