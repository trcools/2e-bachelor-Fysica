# Hoofdstuk 3 — Linear Least Squares

In dit hoofdstuk los je problemen op van het type

- **overdetermined**: $m>n$ (meer vergelijkingen dan onbekenden), typisch bij data/metingen,
- waarbij het stelsel $\mathbf{A}\mathbf{x}=\mathbf{b}$ **geen exacte oplossing** heeft.

Je zoekt dan $\mathbf{x}$ die de **residual** klein maakt:
$$
\mathbf{r}(\mathbf{x})=\mathbf{b}-\mathbf{A}\mathbf{x},
\qquad
\min_{\mathbf{x}} \|\mathbf{r}(\mathbf{x})\|_2.
$$

De examenfocus is bijna altijd: **welke methode wanneer en waarom** (kost vs stabiliteit), plus de connectie met H2 (factorisaties) en H4 (SVD/eigen).

---

## 3.1 Wat is een least squares probleem?

### Overdetermined: $m>n$
$\mathbf{A}\in\mathbb{R}^{m\times n}$, $\mathbf{b}\in\mathbb{R}^{m}$, $\mathbf{x}\in\mathbb{R}^{n}$ met $m>n$.

Het LS-probleem is:
$$
\min_{\mathbf{x}\in\mathbb{R}^n} \|\mathbf{b}-\mathbf{A}\mathbf{x}\|_2^2.
$$

Waarom de **kwadraat**-norm? Omdat
- het wiskundig glad is (afleidbaar),
- en fysisch vaak overeenkomt met “energie”/“error energy”.

### Geometrische interpretatie (zeer examenvriendelijk)
De kolommen van $\mathbf{A}$ spannen een subruimte op: $\mathcal{C}(\mathbf{A})$ (column space).

Least squares kiest $\hat{\mathbf{x}}$ zodat $\mathbf{A}\hat{\mathbf{x}}$ de **orthogonale projectie** is van $\mathbf{b}$ op $\mathcal{C}(\mathbf{A})$.

Dus in de optimum geldt:
$$
\mathbf{r}=\mathbf{b}-\mathbf{A}\hat{\mathbf{x}}
\perp
\mathcal{C}(\mathbf{A})
\quad\Longleftrightarrow\quad
\mathbf{A}^\mathsf{T}\mathbf{r}=\mathbf{0}.
$$

---

## 3.2 Normal equations (klassiek, goedkoop, maar minder stabiel)

Uit $\mathbf{A}^\mathsf{T}\mathbf{r}=\mathbf{0}$ volgt:
$$
\mathbf{A}^\mathsf{T}(\mathbf{b}-\mathbf{A}\hat{\mathbf{x}})=\mathbf{0}
\quad\Rightarrow\quad
(\mathbf{A}^\mathsf{T}\mathbf{A})\hat{\mathbf{x}}=\mathbf{A}^\mathsf{T}\mathbf{b}.
$$

Dit heet de **normal equations**.

### Wanneer werkt dit mooi?
Als $\mathbf{A}$ **full column rank** heeft (rank $n$), dan is $\mathbf{A}^\mathsf{T}\mathbf{A}$:
- symmetric,
- **positive definite**,
dus je kunt **Cholesky** gebruiken (link met H2):
$$
\mathbf{A}^\mathsf{T}\mathbf{A}=\mathbf{L}\mathbf{L}^\mathsf{T}.
$$

### Waarom is dit numeriek riskanter?
Belangrijk inzicht (cursusklassieker):

- In 2-norm geldt typisch:
$$
\mathrm{cond}_2(\mathbf{A}^\mathsf{T}\mathbf{A})=\mathrm{cond}_2(\mathbf{A})^2.
$$

Dus je **condition number kwadrateert** → rounding errors worden veel sterker versterkt dan bij methodes die rechtstreeks met $\mathbf{A}$ werken.

### Samengevat
- **Pro**: relatief goedkoop (zeker als $m\gg n$) en eenvoudig; Cholesky is snel.
- **Con**: kan slecht zijn bij ill-conditioned $\mathbf{A}$; gevaarlijk bij bijna-rank-deficient.

---

## 3.3 QR-factorisatie (standaard “goede keuze” voor LS)

Je factoriseert:
$$
\mathbf{A}=\mathbf{Q}\mathbf{R},
$$
waar
- $\mathbf{Q}\in\mathbb{R}^{m\times n}$ kolom-orthonormaal is (in “thin QR”),
- $\mathbf{R}\in\mathbb{R}^{n\times n}$ upper triangular.

Dan wordt:
$$
\|\mathbf{b}-\mathbf{A}\mathbf{x}\|_2
=
\|\mathbf{b}-\mathbf{Q}\mathbf{R}\mathbf{x}\|_2.
$$

Omdat $\mathbf{Q}$ orthonormaal is, behoudt het 2-normen onder transformatie met $\mathbf{Q}^\mathsf{T}$ (isometrie). Je krijgt:
$$
\min_{\mathbf{x}}\|\mathbf{Q}^\mathsf{T}\mathbf{b}-\mathbf{R}\mathbf{x}\|_2.
$$

In het full-rank geval los je gewoon het triangulaire stelsel:
$$
\mathbf{R}\hat{\mathbf{x}}=\mathbf{Q}^\mathsf{T}\mathbf{b}
$$
met back substitution.

### Hoe bouw je QR numeriek stabiel?
De cursus benadrukt typisch **Householder transformations** (reflecties):
- ze zijn orthogonaal,
- ze zijn numeriek zeer stabiel,
- ze elimineren kolom per kolom tot je $\mathbf{R}$ hebt.

(Alternatief: Givens rotations, vooral handig als je sparseness wil behouden, maar Householder is de default in dense LS.)

### Waarom QR beter is dan normal equations?
- QR vermijdt $\mathbf{A}^\mathsf{T}\mathbf{A}$ en dus het kwadrateren van de conditioning.
- Orthogonale transformaties zijn “rounding-friendly”.

---

## 3.4 SVD (duurste, maar “gouden standaard” voor diagnose en rank issues)

Singular value decomposition:
$$
\mathbf{A}=\mathbf{U}\mathbf{\Sigma}\mathbf{V}^\mathsf{T},
$$
met
- $\mathbf{U}\in\mathbb{R}^{m\times m}$ orthogonaal,
- $\mathbf{V}\in\mathbb{R}^{n\times n}$ orthogonaal,
- $\mathbf{\Sigma}$ diagonaal (singular values $\sigma_1\ge \sigma_2\ge\dots\ge 0$).

### Least squares oplossing via pseudoinverse
De (minimum-norm) LS-oplossing kan geschreven worden met de pseudoinverse:
$$
\hat{\mathbf{x}}=\mathbf{A}^+\mathbf{b}
=
\mathbf{V}\mathbf{\Sigma}^+\mathbf{U}^\mathsf{T}\mathbf{b},
$$
waar $\mathbf{\Sigma}^+$ de inverse neemt van niet-nul singular values (en 0 laat staan).

### Waarom SVD zo nuttig?
- Detecteert **rank deficiency**: kleine $\sigma_i$ betekenen “bijna afhankelijk”.
- Geeft de beste numerieke controle over ill-conditioning.
- Maakt regularisatie (zoals truncated SVD) conceptueel makkelijk.

### Nadeel
- Kostbaar (grootste constante factoren). Gebruik SVD wanneer je *moet* (stabiliteit/diagnose), niet standaard voor alles.

---

## 3.5 Rank bepalen (rechtstreeks link met voorbeeldvraag Q6)

Je zoekt een **numerieke rank**: hoeveel richtingen zijn “significant” boven floating-point noise?

### Methode 1: SVD (meest robuust)
- Compute $\sigma_1\ge\dots\ge\sigma_n$.
- Kies tolerance, bv. $\sigma_i > \tau$ met $\tau$ typisch gekoppeld aan machine precision en schaal (in de cursus vaak “relative threshold”).
- Rank = aantal singular values boven de threshold.

**Voordeel**: zeer betrouwbaar.  
**Nadeel**: duur.

### Methode 2: Rank-revealing QR (QR met pivoting)
Je doet kolompivoting:
$$
\mathbf{A}\mathbf{P}=\mathbf{Q}\mathbf{R},
$$
waar $\mathbf{P}$ kolommen herschikt zodat diagonaal van $\mathbf{R}$ afneemt.
Dan lees je rank af uit de grootte van $|r_{ii}|$ (met tolerance).

**Voordeel**: goedkoper dan SVD, vaak “goed genoeg”.  
**Nadeel**: minder robuust dan SVD in lastige gevallen.

---

## 3.6 Kostvergelijking (rechtstreeks link met voorbeeldvraag Q8)

Voor $\mathbf{A}\in\mathbb{R}^{m\times n}$ met $m>n$:

1) **Normal equations**: vorm $\mathbf{A}^\mathsf{T}\mathbf{A}$ en $\mathbf{A}^\mathsf{T}\mathbf{b}$ en los op (Cholesky)
- typisch goedkoopste in flops,
- maar slechtste in stabiliteit (conditioning kwadrateert).

2) **QR via Householder**
- iets duurder dan normal equations,
- veel stabieler.

3) **SVD**
- duurst,
- meest robuust (zeker bij rank deficiency / ill-conditioning).

**Dus voor $\mathbf{A}$ van grootte $2n\times n$ (zoals in de voorbeeldvraag):**
$$
\text{Normal equations (min work)} \;<\; \text{QR (middel)} \;<\; \text{SVD (max work)}.
$$

---

## 3.7 “Welke methode wanneer?” (de beslisboom die je mondeling moet kunnen verdedigen)

### Kies normal equations + Cholesky als…
- je $\mathbf{A}$ redelijk **goed geconditioneerd** is,
- je puur snelheid wil,
- je weet/verwacht dat rank issues geen rol spelen,
- en je liefst $\mathbf{A}^\mathsf{T}\mathbf{A}$ toch al nodig hebt (maar wees eerlijk over stabiliteit).

**Argument**: goedkoop; SPD $\Rightarrow$ Cholesky snel (link H2).

### Kies QR (Householder) als default voor LS als…
- je een “normale” LS-fit doet op meetdata,
- je stabiliteit wil zonder extreme kost,
- je geen expliciete rank deficiency verwacht.

**Argument**: orthogonale transformaties zijn stabiel; geen conditioning-kwadratering.

### Kies SVD als…
- je rank deficiency vermoedt (multicollineariteit, bijna lineair afhankelijke kolommen),
- je condition number groot is,
- je een betrouwbare numerieke rank wil,
- je een diagnose/regularisatie nodig hebt.

**Argument**: singular values geven meteen inzicht + beste numerieke robuustheid.

---

## Links naar andere hoofdstukken / notebooks (type-1 examenvragen)

### Link met Hoofdstuk 2 (LU vs QR, precies zoals voorbeeldvraag Q3)
- **LU** is de standaard voor square $\mathbf{A}\mathbf{x}=\mathbf{b}$ (direct solve).
- **QR** is de standaard voor least squares omdat je de LS-structuur benut en stabiliteit wint.
- Je *kan* een square systeem ook met QR oplossen: $\mathbf{A}=\mathbf{Q}\mathbf{R}\Rightarrow \mathbf{R}\mathbf{x}=\mathbf{Q}^\mathsf{T}\mathbf{b}$.
  - Voordeel: betere numerieke stabiliteit in sommige gevallen (orthogonale transformaties).
  - Nadeel: meestal duurder dan LU voor pure square solve.

### Link met Hoofdstuk 1 (numerical limitations)
Normal equations zijn het klassieke voorbeeld waar je *een wiskundig nette stap* zet die numeriek slecht is omdat $\mathrm{cond}$ kwadrateert.

### Link met eigenwaarden/SVD notebook
SVD is de “rank/conditioning lens” en verklaart waarom sommige LS-problemen inherent gevoelig zijn.

### Link met optimization notebook
Least squares is een speciaal geval van convex optimization:
$$
\min_{\mathbf{x}} \|\mathbf{b}-\mathbf{A}\mathbf{x}\|_2^2.
$$
Later zie je dezelfde ideeën terug (gradients/Hessians, condition numbers, step choices).

---
