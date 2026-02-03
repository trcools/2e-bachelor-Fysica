# Samenvatting Cursus Py4sci

Deze hoofdreadme bevat de samenvattingen van Hoofdstuk 1 t/m Hoofdstuk 12, samengevoegd uit de per-hoofdstuk bestanden.

---

## Hoofdstuk 1 — Numerical Limitations (waarom numeriek rekenen soms “vals speelt”)

Dit hoofdstuk is de *grondwet* van de hele cursus: elke numerieke methode is een trade-off tussen **benaderingsfout** (truncation) en **rekenfout** (rounding). Als je op het examen moet uitleggen *waarom methode A beter is dan B*, zit het argument bijna altijd hier: **stabiliteit + conditioning + foutopbouw + kost**.

---

---

## Hoofdstuk 2 — Systems of Linear Equations

Dit hoofdstuk leert je **lineaire stelsels** oplossen
$$
\mathbf{A}\mathbf{x}=\mathbf{b},
$$
en vooral: **hoe je dat doet op een manier die (i) snel is en (ii) numeriek stabiel** is in floating-point (link met H1).

---

---

## Hoofdstuk 3 — Linear Least Squares

In dit hoofdstuk los je problemen op van het type

- **overdetermined**: $m>n$ (meer vergelijkingen dan onbekenden), typisch bij data/metingen,
- waarbij het stelsel $\mathbf{A}\mathbf{x}=\mathbf{b}$ **geen exacte oplossing** heeft.

---

---

## Hoofdstuk 4 — Eigenvalue Problems

In dit hoofdstuk draait alles rond het eigenprobleem
$$
\mathbf{A}\mathbf{x}=\lambda \mathbf{x},
$$
waar $\lambda$ een **eigenwaarde** is en $\mathbf{x}\neq \mathbf{0}$ een bijhorende **eigenvector**.

---

---

## Hoofdstuk 5 — Nonlinear equations

We willen een vergelijking oplossen van de vorm
$$
f(x)=0
$$
(of in meerdere dimensies: $\mathbf{f}(\mathbf{x})=\mathbf{0}$). Het grote verschil met lineaire stelsels is dat:
- er **meerdere** oplossingen kunnen zijn (of geen),
- het **vinden** van een oplossing vaak iteratief gaat,
- stabiliteit/stopcriteria niet triviaal zijn.

---

---

## Hoofdstuk 6 — Optimization

Optimalisatie: vind $\mathbf{x}$ die een functie minimaliseert/maximaliseert:
$$
\min_{\mathbf{x}} \ \phi(\mathbf{x})
\quad\text{of}\quad
\max_{\mathbf{x}} \ \phi(\mathbf{x}).
$$

---

---

## Hoofdstuk 7 — Interpolation

Interpolatie: gegeven discrete data $(x_i,y_i)$, construeer een functie $p(x)$ zodat
$$
p(x_i)=y_i \quad \text{voor alle } i.
$$

---

---

## Hoofdstuk 8 — Numerical Integration and Differentiation

Doel:
- integraal:
	$$
	I=\int_a^b f(x)\,dx
	$$
- afgeleide:
	$$
	f'(x)
	$$

---

---

## Hoofdstuk 9 — Ordinary Differential Equations (ODEs)

We lossen een initial value problem (IVP):
$$
\mathbf{y}'(t)=\mathbf{f}(t,\mathbf{y}(t)),
\qquad
\mathbf{y}(t_0)=\mathbf{y}_0,
$$
waar $\mathbf{y}(t)$ een vector is.

---

---

## Hoofdstuk 10 — Partial Differential Equations (PDEs)

PDE’s beschrijven velden $u(\mathbf{x},t)$ met ruimte én tijd. Numeriek komt het altijd neer op:
1) discretiseer ruimte (finite differences / finite elements / spectral),
2) kies een tijdstapper (expliciet of impliciet),
3) controleer stabiliteit (CFL / stiffness),
4) los per stap een lineair of nonlineair stelsel op.

---

---

## Hoofdstuk 11 — Fast Fourier Transform (FFT)

FFT is een snelle manier om de DFT te berekenen. Het is niet alleen “spectra plotten”; het is een **reken-truc** om dingen te versnellen die anders $O(N^2)$ kosten (convolutie/correlatie) en een **analyse-tool** (frequenties, filtering, aliasing).

---

---

## Hoofdstuk 12 — Monte Carlo

Monte Carlo (MC) methodes lossen numerieke problemen op via steekproeven. Het “superpower”-argument:
- de fout schaalt typisch als $1/\sqrt{N}$, bijna onafhankelijk van dimensie,
- dus bij hoge dimensie (waar deterministische quadrature explodeert) is MC vaak de enige praktische optie.

---

*De volledige hoofdstukken zijn bewaard in de per-hoofdstuk bestanden onder* `Hoofdstuk-01` t/m `Hoofdstuk-12`.

---
