<!-- file: hoofdstuk_06_optimization.md -->
# Hoofdstuk 6 — Optimization

Optimalisatie: vind $\mathbf{x}$ die een functie minimaliseert/maximaliseert:
$$
\min_{\mathbf{x}} \ \phi(\mathbf{x})
\quad\text{of}\quad
\max_{\mathbf{x}} \ \phi(\mathbf{x}).
$$

Examenkern:
- link met H5: optimaliteit $\Rightarrow \nabla\phi(\mathbf{x})=\mathbf{0}$ is een nonlineair stelsel,
- link met H3: least squares is een speciaal geval,
- vooral: **wanneer welke methode** (kost, afgeleiden, convexiteit, dimension, constraints).

---

## 6.1 Introduction

Optimalisatie komt in fysica vaak als:
- energie-minimum (stabiel evenwicht),
- parameterfitting (model ↔ data),
- control/trajectory problems (project).

---

## 6.2 Optimality conditions

### 1D
- Stationair punt: $\phi'(x^*)=0$.
- Minimum: $\phi''(x^*)>0$ (lokaal).

### Multi-D
- Stationair punt: $\nabla \phi(\mathbf{x}^*)=\mathbf{0}$.
- Hessiaan:
  $$
  \mathbf{H}(\mathbf{x})=\nabla^2 \phi(\mathbf{x})
  $$
  Minimum als $\mathbf{H}$ positief definiet is.

**Link met H5:** “solve $\nabla\phi=0$” is exact een nonlineair systeem.

---

## 6.3 Optimization in one dimension

### (A) Bracketing van een minimum
Je zoekt een interval $(a,b,c)$ met $a<b<c$ en $\phi(b)<\phi(a),\phi(c)$.

### (B) Golden section search (robust, geen afgeleiden)
Gebruikt een vaste ratio zodat je één functie-evaluatie hergebruikt per iteratie.
- Pro: gegarandeerde krimp van interval; enkel $\phi$ nodig
- Con: lineaire convergentie; relatief veel evaluaties

### (C) Parabolische interpolatie (sneller, maar kan misgaan)
Fit een parabool door drie punten en neem het toppunt.
- Pro: vaak snel bij “gladde” minima
- Con: kan instabiel zijn zonder safeguards

### (D) Brent voor minima
Combineert golden section (veilig) met paraboolstappen (snel).
**Dit is vaak de default** als je 1D zonder afgeleiden doet.

---

## 6.4 Multidimensional unconstrained optimization

Hier is de keuze vooral: heb je afgeleiden (grad/Hessiaan) en hoe duur is een functie-evaluatie?

### (A) Direct search: Nelder–Mead (zoals examenvraag Q5)
Werkt met een simplex in $\mathbb{R}^n$ (n+1 punten) en gebruikt:
- reflectie,
- expansie,
- contractie,
- shrink.

**Waarom convergeert het vaak traag?**
- Het gebruikt enkel functiewaarden (geen richtinginfo),
- simplex-deformaties kunnen “slenteren” door vlakke valleien,
- geen echte curvature-informatie.

**Hoe verbeteren?**
- Gebruik gradient-based methodes (BFGS, CG, Newton),
- of combineer: Nelder–Mead om “in de buurt” te komen, daarna BFGS/Newton.

**Wanneer toch nuttig?**
- Als $\phi$ noisy/niet-differentieerbaar is,
- als gradients niet beschikbaar of onbetrouwbaar zijn,
- als dimension klein is en evaluaties goedkoop zijn.

---

### (B) Gradient descent (goedkoop per iteratie, kan veel iteraties vragen)
Stap:
$$
\mathbf{x}_{k+1}=\mathbf{x}_k - \alpha_k \nabla \phi(\mathbf{x}_k).
$$
- Pro: simpel; enkel gradient nodig
- Con: sterk afhankelijk van step size; traag bij slechte conditioning (langgerekte dalen)

**Link met H3/H1:** conditioning van Hessiaan bepaalt hoe “scheef” het dal is → hoe traag GD is.

---

### (C) Newton’s method (snel, maar duur per iteratie)
Newton-stap:
$$
\mathbf{H}(\mathbf{x}_k)\Delta\mathbf{x}=-\nabla\phi(\mathbf{x}_k),
\qquad
\mathbf{x}_{k+1}=\mathbf{x}_k+\Delta\mathbf{x}.
$$
- Pro: zeer snel dichtbij minimum (curvature correct)
- Con: Hessiaan bouwen/inverteren is duur; solve kost typisch $\mathcal{O}(n^3)$ (dense); kan naar saddle gaan als Hessiaan niet PD is

Praktisch vaak met:
- line search of trust region,
- Hessiaan-regularisatie.

---

### (D) Quasi-Newton (BFGS / L-BFGS) — vaak de “sweet spot”
BFGS bouwt een benadering $\mathbf{B}_k\approx \mathbf{H}^{-1}$ uit gradient-info.
- Pro: bijna Newton-snelheid zonder Hessiaan; zeer populair
- L-BFGS: geheugen-slim voor grote $n$

---

### (E) Conjugate Gradient (CG) voor grote problemen
Voor kwadratische doelen (of bijna) en grote $n$:
- Pro: lage memory; goede scaling
- Con: minder robuust buiten convex/kwadratisch gebied

---

## 6.5 Non-linear Least Squares (koppeling met H3)
Doel:
$$
\min_{\mathbf{x}} \ \|\mathbf{r}(\mathbf{x})\|_2^2,
\qquad
\mathbf{r}:\mathbb{R}^n\to\mathbb{R}^m.
$$
Dit is de niet-lineaire versie van H3.

### (A) Gauss–Newton
Lineariseer residual:
$$
\mathbf{r}(\mathbf{x}+\Delta\mathbf{x})\approx \mathbf{r}(\mathbf{x})+\mathbf{J}(\mathbf{x})\Delta\mathbf{x}.
$$
Minimaliseer dan de lineaire LS:
$$
\min_{\Delta\mathbf{x}} \ \|\mathbf{r}(\mathbf{x})+\mathbf{J}\Delta\mathbf{x}\|_2^2.
$$
Leidt tot normal equations:
$$
(\mathbf{J}^\mathsf{T}\mathbf{J})\Delta\mathbf{x}=-\mathbf{J}^\mathsf{T}\mathbf{r}.
$$

**Link met H3:** dit is exact “LS in elke iteratie” → QR/SVD argumenten blijven gelden.

### (B) Levenberg–Marquardt (demping/regularisatie)
Voegt een term toe:
$$
(\mathbf{J}^\mathsf{T}\mathbf{J}+\lambda \mathbf{I})\Delta\mathbf{x}=-\mathbf{J}^\mathsf{T}\mathbf{r}
$$
om stabiliteit te verbeteren (vooral als $\mathbf{J}^\mathsf{T}\mathbf{J}$ slecht geconditioneerd is).

---

## 6.6 Constrained optimization

Constraints kunnen zijn:
- bounds: $l_i\le x_i\le u_i$,
- equality: $g(\mathbf{x})=0$,
- inequality: $h(\mathbf{x})\le 0$.

Conceptueel:
- Lagrange multipliers / KKT-voorwaarden (voor inzicht),
- numeriek: penalty/barrier of gespecialiseerde algoritmes (SLSQP, trust-constr, etc.).

---

## Welke methode wanneer?

- **Geen afgeleiden / noisy functie / klein n:** Nelder–Mead.
- **Glad + gradients beschikbaar + medium n:** (L-)BFGS is vaak top.
- **Zeer klein n + Hessiaan beschikbaar:** Newton (met line search/trust region).
- **Niet-lineaire LS (model fit):** Gauss–Newton / Levenberg–Marquardt.
- **Constraints:** kies een constrained solver (bounds: L-BFGS-B; algemene constraints: SLSQP of trust-constr).

---

## Links met andere hoofdstukken
- H5: optimaliteit = nulpunten van $\nabla\phi$.
- H3: (non-)linear least squares als kernsubroutine.
- H2: Newton/Gauss–Newton stappen lossen lineaire systemen op; factorisatie-keuzes bepalen runtime.
