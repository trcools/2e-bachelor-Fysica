<!-- file: hoofdstuk_05_nonlinear_equations.md -->
# Hoofdstuk 5 — Nonlinear equations

We willen een vergelijking oplossen van de vorm
$$
f(x)=0
$$
(of in meerdere dimensies: $\mathbf{f}(\mathbf{x})=\mathbf{0}$). Het grote verschil met lineaire stelsels is dat:
- er **meerdere** oplossingen kunnen zijn (of geen),
- het **vinden** van een oplossing vaak iteratief gaat,
- stabiliteit/stopcriteria niet triviaal zijn.

De examenkern (zie voorbeeldvragen): **hoe verbeteren methodes elkaar conceptueel** en **wanneer kies je wat**.

---

## 5.1 Introduction

Niet-lineaire nulpunten duiken overal op in fysica:
- evenwichtscondities (krachtbalans, momentbalans),
- energie-minima (zelfde als optimalisatie: $\nabla \phi(\mathbf{x})=\mathbf{0}$),
- impliciete tijdstappen (ODE/PDE impliciet $\Rightarrow$ telkens een nonlineair solve-probleem).

---

## 5.2 Number of solutions

### Bracketing (1D)
Als $f$ continu is en je vindt $a<b$ met
$$
f(a)\,f(b)<0,
$$
dan bestaat er minstens één wortel in $(a,b)$ (Intermediate Value Theorem).

**Waarom dit goud waard is:** bracketing-methodes zijn *robuust* (ze garanderen een wortel zolang de continuïteit aanneemt).

### Meerdere wortels
- Als $f$ niet-monotoon is, kunnen er meerdere nulpunten zijn.
- In de praktijk: scan/plot om intervallen te vinden waar $f$ van teken verandert.

---

## 5.3 Sensitivity (conditioning van een wortel)

Neem een eenvoudige wortel $x^*$ met $f(x^*)=0$ en $f'(x^*)\neq 0$.
Een kleine perturbatie $\delta f$ in de functie kan een wortelshift $\delta x$ geven met (eerste orde):
$$
0 = f(x^*+\delta x)+\delta f \approx f(x^*)+f'(x^*)\delta x + \delta f
\quad\Rightarrow\quad
\delta x \approx -\frac{\delta f}{f'(x^*)}.
$$

**Interpretatie:**
- Als $|f'(x^*)|$ klein is (platte functie), is de wortel **slecht geconditioneerd**: minieme fouten → grote $\delta x$.
- Bij een meervoudige wortel (typisch $f'(x^*)=0$) wordt het nog gevoeliger én convergentie wordt trager.

Link met H1/H2: zelfs een perfecte iteratiemethode kan niet “beter zijn dan de conditioning” van het probleem.

---

## 5.4 Convergence rates and stopping criteria

### Convergentie-orde
Je meet vaak:
$$
|e_{k+1}| \approx C |e_k|^p,
\quad\text{met } e_k=x_k-x^*.
$$
- $p=1$: lineair
- $p=2$: kwadratisch (Newton)
- $p\approx 1.618$: secant

### Stopcriteria (praktisch én examenvriendelijk)
Gebruik zelden maar één criterium; combineer typisch:
- **Interval** klein: $|b-a|\le \text{tol}$
- **Stap** klein: $|x_{k+1}-x_k|\le \text{tol}(1+|x_{k+1}|)$
- **Residual** klein: $|f(x_k)|\le \text{tol}$

Let op: kleine residual garandeert geen kleine fout als $f'(x^*)$ klein is (conditioning).

---

## 5.5 Solving nonlinear equations in one dimension

### (A) Bisection (bracketing, altijd veilig)
Start met $[a,b]$ met $f(a)f(b)<0$. Neem
$$
m=a+\frac{b-a}{2}.
$$
Kies het subinterval waar het teken verandert.

**Waarom die lijn exact zo staat** (zoals in examenvraag Q7):  
$m=a+(b-a)/2$ vermijdt kleine “numerieke” verrassingen en is de standaard mid-point (zelfde als $(a+b)/2$ maar explicieter in floating-point redenering).

- Convergentie: lineair, fout halveert per iteratie
- Pro: gegarandeerde convergentie (bij continuïteit)
- Con: traag

**Wanneer gebruiken?**
- Als je zekerheid wil (fysica-probleem waar “geen gekke dingen” mogen gebeuren).
- Als $f$ evalueren goedkoop is.
- Als je initieel enkel een bracket hebt, geen goede startwaarde.

---

### (B) Secant method (open method, sneller, geen afgeleide nodig)
Gebruik twee startpunten $x_{k-1},x_k$ en benader $f'$ door een secant:
$$
x_{k+1} = x_k - f(x_k)\frac{x_k-x_{k-1}}{f(x_k)-f(x_{k-1})}.
$$

- Convergentie: superlineair ($p\approx 1.618$) bij “mooie” wortels
- Pro: sneller dan bisection; geen $f'$
- Con: geen garantie op behoud van bracket; kan divergeren

**Conceptuele verbetering t.o.v. bisection (zoals Q2.1):**
- bisection gebruikt alleen tekeninfo; secant gebruikt *lokaal lineair model* van $f$ → grotere, “slimmere” stappen.

---

### (C) Newton’s method (afgeleide, heel snel als je dichtbij zit)
Iteratie:
$$
x_{k+1}=x_k-\frac{f(x_k)}{f'(x_k)}.
$$

- Convergentie: kwadratisch bij eenvoudige wortel en goede start
- Pro: extreem efficiënt als $f'$ beschikbaar is en start goed is
- Con: vraagt $f'$; kan divergeren of naar verkeerde wortel schieten; faalt bij $f'(x_k)\approx 0$

**Conceptuele verbetering t.o.v. secant (zoals Q2.2):**
- secant benadert $f'$; Newton gebruikt de echte afgeleide → sneller en betrouwbaarder *als* $f'$ goed berekenbaar is.

---

### (D) Inverse interpolation / inverse quadratic interpolation (IQI)
In plaats van $f(x)$ te interpoleren en nul te zoeken, interpoleer je **$x$ als functie van $f$**:
$$
x \approx q(f).
$$
Dan is de wortel meteen $x^*\approx q(0)$.

- Pro: vaak sneller dan secant; gebruikt geen afgeleiden
- Con: minder robuust; implementatie complexer

**Waarom IQI voor nulpunten, maar “gewone” interpolatie voor minima (link met voorbeeldvraag Q4):**
- Bij een nulpunt is $f=0$ het *doel* → het is logisch om rechtstreeks $x(f=0)$ te schatten.
- Bij een minimum heb je typisch $f'(x)=0$ of je wil $x$ minimaliseren → “parabool-fit” in $x$-ruimte is natuurlijker (je modelleert $f(x)$ als parabool en neemt het toppunt).

---

### (E) Hybride topkeuze: Brent-achtige methodes
In de praktijk wil je:
- de robuustheid van bracketing,
- de snelheid van secant/IQI.

Brent’s methode combineert bracketing + secant + IQI en valt terug op bisection als het gevaarlijk wordt.

**Wanneer kiezen?**
- Default in productie: “zo goed als altijd werkt” zonder veel tuning.

---

## 5.6 Systems of nonlinear equations

We lossen:
$$
\mathbf{f}(\mathbf{x})=\mathbf{0},\qquad \mathbf{f}:\mathbb{R}^n\to\mathbb{R}^n.
$$

### Newton voor systemen
Lineariseer rond $\mathbf{x}_k$:
$$
\mathbf{f}(\mathbf{x}_k+\Delta\mathbf{x})\approx \mathbf{f}(\mathbf{x}_k)+\mathbf{J}(\mathbf{x}_k)\Delta\mathbf{x},
$$
waar $\mathbf{J}$ de Jacobiaan is:
$$
J_{ij}(\mathbf{x})=\frac{\partial f_i}{\partial x_j}.
$$

Stel $\mathbf{f}(\mathbf{x}_{k+1})=\mathbf{0}$:
$$
\mathbf{J}(\mathbf{x}_k)\Delta\mathbf{x}=-\mathbf{f}(\mathbf{x}_k),\qquad
\mathbf{x}_{k+1}=\mathbf{x}_k+\Delta\mathbf{x}.
$$

**Link met H2:** elke Newton-stap is een lineair stelsel solve. Keuze LU/Cholesky/sparse bepaalt je totale kost.

### Praktische stabilisatie
- **Damping/line search:** $\mathbf{x}_{k+1}=\mathbf{x}_k+\alpha\Delta\mathbf{x}$ met $0<\alpha\le 1$ om overshoots te vermijden.
- **Goede start** is essentieel (Newton is lokaal).

---

## Welke methode wanneer?

- **Je hebt alleen een bracket en wil zekerheid:** bisection of Brent.
- **Je hebt geen afgeleide maar wil sneller:** secant of Brent.
- **Je hebt $f'$ (of Jacobiaan) en een goede start:** Newton.
- **Je wil “best of both worlds” zonder te veel gedoe:** Brent (1D) of gedempte Newton (multi-D).

---

## Links met andere hoofdstukken
- H2: Newton voor systemen = herhaald lineaire stelsels oplossen.
- H6: optimalisatie gebruikt vaak solves van $\nabla \phi(\mathbf{x})=\mathbf{0}$ → nonlineair.
- H9/H10: impliciete ODE/PDE-stappen geven nonlineaire systemen.
