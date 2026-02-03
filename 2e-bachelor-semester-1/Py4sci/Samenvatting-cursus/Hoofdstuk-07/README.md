<!-- file: hoofdstuk_07_interpolation.md -->
# Hoofdstuk 7 — Interpolation

Interpolatie: gegeven discrete data $(x_i,y_i)$, construeer een functie $p(x)$ zodat
$$
p(x_i)=y_i \quad \text{voor alle } i.
$$

Examenkern:
- **global polynomial** vs **piecewise**,
- numerieke stabiliteit (Vandermonde is gevaarlijk),
- wanneer gebruik je interpolatie vs least squares (H3),
- en waarom splines vaak winnen bij veel punten.

---

## 7.1 Introduction

Interpolatie is nuttig als je:
- waarden wil “tussenin” schatten,
- een gladde curve nodig hebt voor verdere analyse (bv. integratie, differentiatie, root finding).

Let op: interpolatie forceert exact door meetpunten → bij ruis is dat vaak een slecht idee (dan wil je H3: regression/least squares).

---

## 7.2 Polynomial interpolation of discrete data

### 7.2.1 Lagrange-form
Voor $n+1$ punten:
$$
p(x)=\sum_{i=0}^n y_i \ell_i(x),
\qquad
\ell_i(x)=\prod_{\substack{j=0\\ j\neq i}}^n \frac{x-x_j}{x_i-x_j}.
$$
- Pro: conceptueel helder
- Con: duur/instabiel als je het naïef evalueert

### 7.2.2 Newton-form + divided differences
Newton-vorm:
$$
p(x)=a_0+a_1(x-x_0)+a_2(x-x_0)(x-x_1)+\cdots
$$
Coeffs $a_k$ via divided differences.
- Pro: makkelijk uitbreidbaar met nieuwe punten; efficiënte evaluatie (Horner-achtig)
- Con: nog steeds global polynomial (Runge-risico)

### 7.2.3 Vandermonde-systeem (waarom je dit meestal NIET wil)
Je kan ook oplossen uit:
$$
\mathbf{V}\mathbf{c}=\mathbf{y},
\qquad
V_{ij}=x_i^j.
$$
Maar $\mathbf{V}$ is vaak extreem slecht geconditioneerd → numeriek ellende (link H1/H2/H3).

**Cursusboodschap:** doe dit alleen als didactiek; in echte code: nee.

### 7.2.4 Interpolation error + Runge phenomenon
Foutterm (als $f$ glad genoeg is):
$$
f(x)-p(x)=\frac{f^{(n+1)}(\xi)}{(n+1)!}\prod_{i=0}^n (x-x_i)
$$
voor een $\xi$ tussen de nodes.

**Runge:** bij gelijk verdeelde nodes kan het global polynomial wild oscilleren aan de randen.
Fix: **Chebyshev nodes** (clusteren naar de randen) → veel stabielere global interpolatie.

### 7.2.5 Barycentric interpolation (de numeriek stabiele manier)
Barycentric formule evalueert Lagrange-interpolatie stabiel en snel.
Cursuspunt: dit is de praktische manier als je toch global polynomial wil.

---

## 7.3 Piecewise polynomial interpolation

### Waarom piecewise?
Bij veel punten is één polynoom van hoge graad meestal een ramp (oscillaties, conditioning). Piecewise houdt de graad laag en stabiliteit hoog.

### 7.3.1 Piecewise linear
Snel, robuust, maar niet glad in afgeleide.

### 7.3.2 Cubic splines
Cubic spline: op elk interval een cubisch polynoom met voorwaarden op continuïteit:
- $p$, $p'$, $p''$ continu op knooppunten.
Je hebt extra randvoorwaarden (bv. “natural spline”: tweede afgeleide 0 aan de randen).

- Pro: zeer glad; meestal “default” voor gladde interpolatie
- Con: kan overshoots geven bij monotone data

### 7.3.3 Shape-preserving (PCHIP)
Piecewise Cubic Hermite Interpolating Polynomial (PCHIP) behoudt monotoniciteit (minder overshoot).
- Pro: goed voor data die fysisch monotone moeten blijven (bv. dichtheid, cumulatieven)
- Con: iets minder “glad” dan klassieke cubic spline (maar meestal een feature)

---

## 7.4 Software (SciPy)

Typische tools:
- `scipy.interpolate.interp1d` (linear/quadratic/cubic, eenvoudig)
- `scipy.interpolate.CubicSpline`
- `scipy.interpolate.PchipInterpolator`
- `scipy.interpolate.BarycentricInterpolator`

---

## Welke methode wanneer?

- **Kleine set punten, zeer gladde functie, je wil hoge orde:** barycentric / Chebyshev nodes (global).
- **Veel punten, je wil stabiel en glad:** cubic spline.
- **Monotone fysische data (geen overshoot toegestaan):** PCHIP.
- **Ruis aanwezig / je wil niet exact door punten:** geen interpolatie → H3 least squares fit.

---

## Links met andere hoofdstukken
- H3: interpolatie vs regression (interpolatie forceert exact; LS filtert ruis).
- H8: integratie/differentiatie van geïnterpoleerde data (splines zijn handig).
- H11: Fourier-achtige interpolatie voor periodieke signalen (FFT-wereld).
