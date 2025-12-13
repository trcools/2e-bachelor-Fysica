# Hoofdstuk 2 — Systems of Linear Equations

Dit hoofdstuk leert je **lineaire stelsels** oplossen
$$
\mathbf{A}\mathbf{x}=\mathbf{b},
$$
en vooral: **hoe je dat doet op een manier die (i) snel is en (ii) numeriek stabiel** is in floating-point (link met H1).

De rode draad (en dit is letterlijk de “type 1”-examenvraag zoals Q1):  
**transformeer** het probleem naar een vorm die je **goedkoop** kan oplossen (triangulair), en doe dat met controle op stabiliteit (pivoting).

---

## 2.1 Introduction and Notation

### Lineaire modellen in fysica
In de cursus worden o.a. genoemd:
- Newton: $\mathbf{F}=m\mathbf{a}$
- Ohm: $R=U/I$
- Hooke: $\mathbf{F}_s=k\mathbf{x}$

Algemeen: een systeem van $m$ lineaire vergelijkingen met $n$ onbekenden:
$$
\begin{cases}
a_{11}x_1 + a_{12}x_2 + \dots + a_{1n}x_n = b_1\\
a_{21}x_1 + a_{22}x_2 + \dots + a_{2n}x_n = b_2\\
\vdots\\
a_{m1}x_1 + a_{m2}x_2 + \dots + a_{mn}x_n = b_m
\end{cases}
$$
en in matrixvorm:
$$
\mathbf{A}\mathbf{x}=\mathbf{b},
$$
waar $\mathbf{A}\in\mathbb{R}^{m\times n}$ (of $\mathbb{C}^{m\times n}$), $\mathbf{x}\in\mathbb{R}^n$, $\mathbf{b}\in\mathbb{R}^m$.

---

## 2.2 Solving Linear Systems

### 2.2.1 Triangular systems: forward/back substitution
Als je het stelsel kunt omvormen naar een triangulair stelsel, is oplossen goedkoop.

**Lower triangular** $\mathbf{L}\mathbf{x}=\mathbf{b}$ (forward substitution):
$$
x_1=\frac{b_1}{l_{11}},\qquad
x_i=\frac{b_i-\sum_{j=1}^{i-1}l_{ij}x_j}{l_{ii}}\quad (i=2,\dots,n).
$$

**Upper triangular** $\mathbf{U}\mathbf{x}=\mathbf{b}$ (back substitution):
$$
x_n=\frac{b_n}{u_{nn}},\qquad
x_i=\frac{b_i-\sum_{j=i+1}^{n}u_{ij}x_j}{u_{ii}}\quad (i=n-1,\dots,1).
$$

**Waarom dit belangrijk is:** de hele rest van het hoofdstuk gaat over “hoe krijg ik $\mathbf{A}$ naar (een product van) zulke driehoeksmatrices?”

---

### 2.2.2 Elementary elimination matrices (Gauss transformations)
Gaussian elimination elimineert systematisch elementen onder de diagonaal via rijoperaties.

Bij stap $k$ kies je een pivot (in de cursus: “the pivot”) en elimineer je $a_{ik}$ voor $i>k$ door rij $i$ te vervangen door:
$$
\text{row}_i \leftarrow \text{row}_i + m_i\,\text{row}_k,
\qquad
\text{met } m_i=-\frac{a_{ik}}{a_{kk}}.
$$

Dat kan in matrixvorm met een **elementary elimination matrix** $\mathbf{M}_k$:
- Effect: voeg een veelvoud van rij $k$ toe aan rijen $k+1,\dots,n$ zodat de kolom onder de diagonaal nul wordt.
- In de cursus staat expliciet de vorm:
$$
\mathbf{M}_k=\mathbf{I}-\mathbf{m}_k\mathbf{e}_k^\mathsf{T},
$$
waar $\mathbf{m}_k=[0,\dots,0,m_{k+1},\dots,m_n]^\mathsf{T}$ en $\mathbf{e}_k$ de $k$-de kolom van $\mathbf{I}$ is.

Belangrijke eigenschap (cursus):
$$
\mathbf{M}_k^{-1}=\mathbf{I}+\mathbf{m}_k\mathbf{e}_k^\mathsf{T}.
$$
Dus de inverse heeft **dezelfde structuur** maar met **omgekeerde tekens** (in de cursus wordt $\mathbf{M}_k^{-1}$ vaak aangeduid als $\mathbf{L}_k$).

---

### 2.2.3 Gaussian elimination $\Rightarrow$ LU factorization
Als je alle eliminatiestappen samenneemt:
$$
\mathbf{M}=\mathbf{M}_{n-1}\cdots \mathbf{M}_1,
$$
dan wordt
$$
\mathbf{M}\mathbf{A}=\mathbf{U},
$$
met $\mathbf{U}$ upper triangular.

Definieer dan:
$$
\mathbf{L}=\mathbf{M}^{-1},\qquad \mathbf{U}=\mathbf{M}\mathbf{A}.
$$
Dan krijg je de **LU-factorisatie**:
$$
\mathbf{A}=\mathbf{L}\mathbf{U}.
$$

**Waarom LU zo centraal is (examenvraag-stijl):**
- Factoriseer één keer (duur),
- Los daarna op met twee triangulaire solves (goedkoop):
  1. $\mathbf{L}\mathbf{y}=\mathbf{b}$ (forward)
  2. $\mathbf{U}\mathbf{x}=\mathbf{y}$ (back)

**Extra winst (cursus benadrukt dit):** als je meerdere rechterleden hebt (zelfde $\mathbf{A}$, andere $\mathbf{b}$), hergebruik je $\mathbf{L},\mathbf{U}$ en betaal je per extra $\mathbf{b}$ enkel substituties.

---

### 2.2.4 Partial pivoting (stabiliteit + vermijden van breakdown)
De cursus geeft 2 problemen bij “naïeve” Gaussian elimination:

1) **Breakdown** als pivot $a_{kk}=0$ (je moet delen door 0).  
   Oplossing: wissel rijen zodat je een niet-nul pivot hebt → **pivoting**.

2) **Numerieke instabiliteit** in floating-point: te grote multipliers versterken rounding errors.  
   Oplossing: **partial pivoting**: kies in kolom $k$ de entry met **grootste absolute waarde** op/onder de diagonaal als pivot. Dan blijven multipliers in grootte $\le 1$.

Met pivoting verschijnt een permutatiematrix $\mathbf{P}$ en typisch krijg je (notatie zoals in veel software):
$$
\mathbf{P}\mathbf{A}=\mathbf{L}\mathbf{U}.
$$

De cursus vermeldt ook expliciet een SciPy-conventie: als je $\mathbf{P}$ “in $\mathbf{L}$ absorbeert”, kan je schrijven dat $\mathbf{L}\mathbf{U}=\mathbf{A}$ (met een $\mathbf{L}$ die een permutatie van lower-triangular is). Dit is consistent met `scipy.linalg` met optie `permute_l=True`.

**Praktisch besluit (waarom het “beter” is):** partial pivoting maakt LU **veel robuuster** zonder de orde van de kost te veranderen.

---

### 2.2.5 Gauss–Jordan elimination (waarom meestal niet)
Je kunt ook elimineren tot een **diagonale** matrix (of zelfs $\mathbf{I}$) door *ook boven* de diagonaal weg te werken: dat is Gauss–Jordan.

Cursusboodschap:
- Ja, het kan.
- Maar de extra kost levert meestal niet genoeg voordeel op voor het oplossen van $\mathbf{A}\mathbf{x}=\mathbf{b}$, zeker niet t.o.v. LU + substitutie.

---

## 2.3 Special types of linear systems

### 2.3.1 Symmetric positive definite (SPD) $\Rightarrow$ Cholesky
Als $\mathbf{A}$ **symmetric positive definite** is, bestaat:
$$
\mathbf{A}=\mathbf{L}\mathbf{L}^\mathsf{T}
$$
(Cholesky factorization).

De cursus somt expliciet voordelen op:
- De $n$ wortels zijn van **positieve** getallen → algoritme is goed-gedefinieerd
- **Geen pivoting** nodig
- Je gebruikt enkel de **lower triangle** van $\mathbf{A}$ (minder opslag)
- Kost: ongeveer $n^3/6$ multiplications (en vergelijkbaar aantal additions)

Conclusie (cursus): Cholesky is ongeveer **half zoveel werk en opslag** als algemene LU.

**Wanneer kies je Cholesky?**  
Als je uit de fysica/matrixstructuur kunt argumenteren dat $\mathbf{A}$ SPD is (bv. bepaalde energie/Hessiaan-achtige matrices, normal equations later in H3).

---

### 2.3.2 Computational complexity (kostargumenten die je moet kunnen verwoorden)
De cursus geeft de typische flop-counts:

- LU-factorisatie van $n\times n$: ongeveer $n^3/3$ flops
- Volledige matrix-inversie: ongeveer $n^3$ flops (dus ~3× duurder)
- Oplossen met forward + backward substitution na LU: ongeveer $n^2$ flops (voor grote $n$ verwaarloosbaar t.o.v. factorisatie)
- Cramer’s rule: “astronomisch duur”

**Belangrijkste praktijkregel (cursus zegt dit letterlijk in spirit):**  
Bereken $\mathbf{A}^{-1}$ bijna nooit expliciet; los $\mathbf{A}\mathbf{x}=\mathbf{b}$ op via factorisatie + substitutie (sneller én nauwkeuriger).

---

## 2.4 Sensitivity and Conditioning

Hier leer je het verschil tussen:
- “ik heb een $\mathbf{x}$ uitgerekend”
- “mijn $\mathbf{x}$ is betrouwbaar”

### 2.4.1 Vector norms ($p$-normen zoals in de cursus)
Voor $p>0$:
$$
\|\mathbf{x}\|_p=\left(\sum_{i=1}^n |x_i|^p\right)^{1/p}.
$$

Belangrijke gevallen:
- $p=1$ (Manhattan):
  $$\|\mathbf{x}\|_1=\sum_{i=1}^n |x_i|$$
- $p=2$ (Euclidisch):
  $$\|\mathbf{x}\|_2=\left(\sum_{i=1}^n |x_i|^2\right)^{1/2}$$
- $p=\infty$:
  $$\|\mathbf{x}\|_\infty=\max_i |x_i|$$

### 2.4.2 Matrix norms (induced norms die SciPy ook gebruikt)
In de cursus worden typisch gebruikt:
- 1-norm (max kolomsom):
  $$\|\mathbf{A}\|_1=\max_j \sum_i |a_{ij}|$$
- $\infty$-norm (max rijsom):
  $$\|\mathbf{A}\|_\infty=\max_i \sum_j |a_{ij}|$$

Sleutel-eigenschap voor foutbounds:
$$
\|\mathbf{A}\mathbf{x}\|\le \|\mathbf{A}\|\;\|\mathbf{x}\|.
$$

### 2.4.3 Condition number (cursusnotatie met `cond`)
Voor een gekozen norm:
$$
\mathrm{cond}(\mathbf{A})=\|\mathbf{A}\|\;\|\mathbf{A}^{-1}\|.
$$
Specifiek bv.:
$$
\mathrm{cond}_\infty(\mathbf{A})=\|\mathbf{A}\|_\infty\;\|\mathbf{A}^{-1}\|_\infty.
$$

Interpretatie: hoe hard kan het probleem kleine inputfouten **amplifiëren**?

### 2.4.4 Error estimation (perturbatie in $\mathbf{b}$)
Neem:
$$
\mathbf{A}\mathbf{x}=\mathbf{b},\qquad
\mathbf{A}\mathbf{x}'=\mathbf{b}+\Delta\mathbf{b},
\qquad \Delta\mathbf{x}=\mathbf{x}'-\mathbf{x}.
$$
Dan:
$$
\mathbf{A}\Delta\mathbf{x}=\Delta\mathbf{b}
\quad\Rightarrow\quad
\Delta\mathbf{x}=\mathbf{A}^{-1}\Delta\mathbf{b}.
$$
Normeren geeft de standaard bound (idee zoals in de cursus):
$$
\frac{\|\Delta\mathbf{x}\|}{\|\mathbf{x}\|}
\le
\mathrm{cond}(\mathbf{A})\;
\frac{\|\Delta\mathbf{b}\|}{\|\mathbf{b}\|}.
$$

**Interpretatie (examenvriendelijk):** bij grote $\mathrm{cond}(\mathbf{A})$ kan een klein meetfoutje in $\mathbf{b}$ tot een grote fout in $\mathbf{x}$ leiden.

### 2.4.5 Residual (wat controleer je in de praktijk?)
Als je numeriek $\mathbf{x}_{\text{num}}$ krijgt, definieer:
$$
\mathbf{r}=\mathbf{A}\mathbf{x}_{\text{num}}-\mathbf{b}.
$$
Dan:
$$
\mathbf{x}_{\text{num}}-\mathbf{x}=-\mathbf{A}^{-1}\mathbf{r}
\quad\Rightarrow\quad
\|\mathbf{x}_{\text{num}}-\mathbf{x}\|\le \|\mathbf{A}^{-1}\|\;\|\mathbf{r}\|.
$$

**Belangrijke nuance:** een kleine residual betekent pas “goede oplossing” als $\mathbf{A}$ niet slecht geconditioneerd is.

---

## 2.5 Software (SciPy/NumPy zoals in de cursus)

- Solve (algemeen): `scipy.linalg.solve(A, b)`
- LU (met pivoting): `scipy.linalg.lu(A)` of efficiënter `lu_factor` + `lu_solve`
- Cholesky (SPD): `scipy.linalg.cholesky(A, lower=True)` en dan twee substituties
- Normen/residuals: `scipy.linalg.norm(...)`

---

## “Welke methode wanneer?” (de mapping fysica-probleem → methode)

1) **Algemeen dense $n\times n$ stelsel**
   - Kies: LU met partial pivoting (standaard in `solve`)
   - Waarom: robuust; kost $\sim n^3/3$ flops

2) **Zelfde $\mathbf{A}$, veel verschillende $\mathbf{b}$**
   - Kies: één keer LU-factorisatie, daarna herhaald `lu_solve`
   - Waarom: factorisatie duur, solve per $\mathbf{b}$ goedkoop ($\sim n^2$)

3) **$\mathbf{A}$ is symmetric positive definite**
   - Kies: Cholesky $\mathbf{A}=\mathbf{L}\mathbf{L}^\mathsf{T}$
   - Waarom: geen pivoting, halve opslag, $\sim n^3/6$ multiplications

4) **Iemand wil $\mathbf{A}^{-1}$ expliciet**
   - Meestal: niet doen
   - Waarom (cursus): $\sim n^3$ flops en minder nauwkeurig dan factorisatie + substitutie

5) **Je vertrouwt de oplossing niet**
   - Check: residual $\mathbf{r}$ én een idee van $\mathrm{cond}(\mathbf{A})$
   - Waarom: kleine residual is niet voldoende bij slechte conditioning

---

## Links naar andere notebooks (zoals “type 1”-vragen)

- **Naar H1 Numerical limitations:** pivoting en “stabiele factorisaties” zijn een direct antwoord op rounding-error amplificatie.
- **Naar H3 Linear least squares:** SPD en Cholesky komen terug via normal equations; en de keuze LU vs QR vs SVD is precies “kost vs stabiliteit”.
- **Naar eigenwaarden/SVD later:** conditioning en (bijna) singulariteit worden daar “structureel zichtbaar” via singular values.

---
