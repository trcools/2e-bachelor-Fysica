<!-- file: hoofdstuk_08_numerical_integration_and_differentiation.md -->
# Hoofdstuk 8 — Numerical Integration and Differentiation

Doel:
- integraal:
  $$
  I=\int_a^b f(x)\,dx
  $$
- afgeleide:
  $$
  f'(x)
  $$

Examenkern:
- foutorde vs kost,
- (on)stabiliteit van differentiatie,
- adaptiviteit,
- Richardson extrapolation als “order booster”.

---

## 8.1 Integration

Integratie is meestal numeriek **vriendelijker** dan differentiatie (differentiatie versterkt ruis).

---

## 8.2 Existence, Uniqueness, Conditioning

Voor integralen:
- bestaan/uniqueness is zelden het probleem (voor nette $f$),
- conditioning hangt vaak af van hoe “wild” of “singulier” $f$ is, en van cancellation.

---

## 8.3 Numerical Quadrature

We bouwen een benadering:
$$
\int_a^b f(x)\,dx \approx \sum_{i=0}^N w_i f(x_i).
$$

### 8.3.1 Newton–Cotes (equispaced)
- Rectangle/left/right: lage orde
- Midpoint: betere orde
- Trapezoidal rule
- Simpson’s rule (parabool-fit)

**Orders (typisch in de cursus):**
- midpoint en trapezoid: fout $O(h^2)$
- Simpson: fout $O(h^4)$
waar $h$ de stapgrootte is.

**Trade-off:**
- hogere orde = minder stappen voor dezelfde truncation error,
- maar bij te klein $h$ kan rounding meespelen (link H1).

### 8.3.2 Composite rules
Splits $[a,b]$ in subintervallen en pas de regel herhaald toe.
- Werkt goed voor “redelijk gladde” $f$.

### 8.3.3 Gaussian quadrature (hoog rendement voor gladde $f$)
Kies nodes/weights zodat de regel exact is voor polynomen tot een hoge graad.
- Pro: extreem efficiënt als $f$ glad is
- Con: minder plug-and-play als $f$ singular/oscillerend is (maar er bestaan varianten)

### 8.3.4 Adaptive quadrature
Pas het interval dynamisch aan: verfijn waar $f$ moeilijk is.
- Pro: efficiënt bij lokale pieken/singulariteit
- Con: overhead/complexiteit; maar vaak de beste “algemene” keuze

---

## 8.4 Other integration problems

### Improper integrals / singular integrands
Strategieën:
- variabeletransformatie om singulariteit te verzwakken,
- splits interval rond de singulariteit,
- adaptive methods.

### Oscillerende integralen
Standaard Newton–Cotes kan falen.
Mogelijke routes:
- speciale oscillatory quadrature,
- herformuleren (Fourier/FFT-ideeën bij periodieke structuren).

### Multi-dimensional integrals
Deterministische quadrature explodeert in kost met dimensie (“curse of dimensionality”).
→ Monte Carlo (H12) wordt dan aantrekkelijk.

---

## 8.5 Numerical Differentiation

### Finite differences
Forward:
$$
f'(x)\approx \frac{f(x+h)-f(x)}{h}
\quad \text{(truncation } O(h)\text{)}
$$
Central:
$$
f'(x)\approx \frac{f(x+h)-f(x-h)}{2h}
\quad \text{(truncation } O(h^2)\text{)}
$$

### Cruciale trade-off: truncation vs rounding
- truncation daalt als $h\to 0$,
- maar rounding/cancellation stijgt omdat je bijna gelijke getallen aftrekt:
  $f(x+h)-f(x)$.

Gevolg: er bestaat een **optimale** $h$; “zo klein mogelijk” is fout (link H1).

### Differentiatie van ruis
Ruis wordt versterkt door afgeleiden → vaak eerst smoothen of een model fitten (H3/H7) en dáárvan differentiëren.

---

## 8.6 Richardson Extrapolation

Idee: combineer twee benaderingen met stap $h$ en $h/2$ om lagere-orde fout weg te elimineren.

Voor een methode met fout:
$$
A(h)=A + C h^p + O(h^{p+1}),
$$
dan:
$$
A \approx \frac{2^p A(h/2)-A(h)}{2^p-1}.
$$

**Gebruik:**
- upgrade trapezoid → Romberg-achtige integratie,
- upgrade differentieformules.

---

## 8.7 SciPy

- `scipy.integrate.quad` (adaptief, 1D)
- `scipy.integrate.quadrature` / `fixed_quad` (Gauss)
- voor ODE-integratie: zie H9 (andere wereld)

---

## Welke methode wanneer?

### Integratie
- **Algemene 1D integralen, weinig gedoe:** adaptief (`quad`-stijl).
- **Zeer glad, hoge nauwkeurigheid met weinig evaluaties:** Gaussian quadrature.
- **Simpele integralen, je wil controle/educatief:** composite trapezoid/Simpson.
- **Hoge dimensie:** ga naar H12 (Monte Carlo).

### Differentiatie
- **Je hebt een analytisch model of fit:** differentieer het model (best).
- **Je hebt zuivere data (weinig ruis):** central differences + Richardson.
- **Ruis aanwezig:** vermijd directe differentiatie; smooth/fit eerst.

---

## Links met andere hoofdstukken
- H7: splines/fit als “voorbewerking” voor differentiatie.
- H9/H10: discretisaties gebruiken exact deze finite differences; stabiliteit hangt af van stapgroottes.
- H12: multi-D integratie → Monte Carlo.
