# Hoofdstuk 1 — Numerical Limitations (waarom numeriek rekenen soms “vals speelt”)

Dit hoofdstuk is de *grondwet* van de hele cursus: elke numerieke methode is een trade-off tussen **benaderingsfout** (truncation) en **rekenfout** (rounding). Als je op het examen moet uitleggen *waarom methode A beter is dan B*, zit het argument bijna altijd hier: **stabiliteit + conditioning + foutopbouw + kost**.

---

## 1.1 Approximations in scientific computation

### Absolute vs relative error
- **Absolute error**:
  $$e_{\text{abs}} = x_{\text{approx}} - x_{\text{true}}$$
- **Relative error**:
  $$e_{\text{rel}} = \frac{x_{\text{approx}} - x_{\text{true}}}{x_{\text{true}}}$$

Interpretatie (examenvriendelijk): als $e_{\text{rel}} \sim 10^{-p}$, dan heb je ongeveer **$p$ correcte significant digits**.

### Precision vs accuracy
- **Precision** = hoeveel digits je *opschrijft/bewaart*.
- **Accuracy** = hoeveel digits *effectief correct zijn*.
Kernboodschap: meer precision garandeert geen accuracy (bv. een lang getal kan toch totaal naast $\pi$ zitten).

### Truncation error vs rounding error
- **Truncation error**: komt uit de *wiskundige benadering* (afkappen van reeksen, finite differences, discretisatie, …).
- **Rounding error**: komt uit het feit dat computers met **eindige precisie** rekenen en voortdurend afronden.

**Belangrijk patroon voor “wanneer welke methode?”**  
Als je de stapgrootte $h$ kleiner maakt (of $N$ groter, meer iteraties, fijnere grid):
- truncation error daalt meestal,
- rounding error en kost stijgen meestal,
- en soms wordt het resultaat zelfs slechter door foutopstapeling.

Dus: “gewoon fijner nemen” is géén universele win.

---

## 1.2 Computer arithmetic

### 1.2.1 Floating-point number systems (het model achter `float`)
Floating-point lijkt op wetenschappelijke notatie:
$$
\pm (d_0 + d_1\beta^{-1} + \dots + d_{p-1}\beta^{-(p-1)})\,\beta^{E}
$$
met:
- $\beta$ = base/radix (typisch 2),
- $p$ = precision (aantal “mantissa digits”),
- $E \in [L,U]$ = exponent range.

**Normalisatie** (zoals in IEEE 754): de leidende digit wordt “vastgezet” (typisch $d_0=1$ voor niet-nul). Daardoor “win” je effectief één bit aan informatie (de bekende *hidden bit* in base 2).

Praktisch in Python:
- `float` is (typisch) **double precision**.
- `numpy` laat ook **single precision** toe (`np.float32`), met minder geheugen maar merkbaar meer rounding error.

### 1.2.2 Properties (discreet, eindig, en soms gemeen)

**(1) Discreet getallenrooster**  
Tussen twee representabele floats zit een “gat”. Veel decimalen bestaan niet exact in binair.

**(2) Overflow/underflow**
- Te groot $\to$ `inf` (overflow).
- Te klein $\to$ afronding naar `0.0` (underflow), eventueel via subnormals/denormals.

**(3) $0.1$ is niet exact**
Klassiek gevolg: `0.1 + 0.1 + 0.1 != 0.3`.  
Dus: `==` is bijna altijd fout bij floats.

**(4) Afronding maakt optellen niet-associatief**
Door rounding geldt vaak:
$$(a+b)+c \neq a+(b+c)$$
Concreet gevolg: **somvolgorde** van een lange reeks verandert de uitkomst (en soms merkbaar).
- “Natural order” kan slechter zijn dan “reverse order”.
- Groepering beïnvloedt cancellation en rounding-opstapeling.

> Examenzin die altijd scoort: *Een algoritme kan wiskundig correct zijn, maar numeriek instabiel doordat rounding zich opstapelt afhankelijk van de route die je neemt.*

### 1.2.3 Good Practices for computer arithmetic (de survival kit)

#### A. Cancellation (catastrophic cancellation)
- **Vermijd**: subtractie van bijna gelijke getallen.  
  Voorbeeldidee: $\sqrt{x+1} - \sqrt{x}$ voor groot $x$ $\to$ twee bijna gelijke grote getallen $\to$ verschil verliest significant digits.
- **Fix**: herschrijf algebraïsch naar een stabiele vorm:
  $$
  \sqrt{x+1} - \sqrt{x} = \frac{1}{\sqrt{x+1} + \sqrt{x}}
  $$

#### B. Addition (small + large)
- **Vermijd**: een heel klein getal optellen bij een enorm groot getal (het kleine kan volledig “verdwijnen” omdat het onder de float-spacing zit).
- **Praktische regel**: sommeer een lijst **van klein naar groot** (of gebruik algoritmes die dit effect beperken).

#### C. Float-vergelijkingen: altijd met tolerantie
- Absolute tolerantie is belangrijk rond 0.
- Relatieve tolerantie is belangrijk bij grote schalen.
In code: `np.isclose` / `math.isclose` i.p.v. `==`.

---

## “Wanneer welke aanpak?” — beslisregels die je echt gebruikt in fysicaproblemen

### 1) Kies je datatype bewust
- **`float64` (double)**: default voor fysica/numeriek werk $\to$ betere accuracy, meestal nog snel genoeg.
- **`float32` (single)**: nuttig bij grote arrays/GPU/memory pressure, maar verwacht:
  - meer rounding,
  - sneller instabiliteit,
  - grotere gevoeligheid voor somvolgorde/cancellation.

### 2) Schaal je probleem (units/normalization) vóór je gaat rekenen
Als grootheden extreem groot/klein zijn, krijg je sneller overflow/underflow en slechtere conditioning.
- Rescale variabelen zodat typische waarden $O(1)$ zijn.
- Dit is geen cosmetica: het kan de numerieke betrouwbaarheid drastisch verbeteren.

### 3) Prefer “numerically stable” formuleringen boven “directe” formuleringen
Zelfde wiskunde, andere route $\to$ andere rounding error.
- Herformuleer uitdrukkingen met subtractie van bijna gelijken.
- Orden sommen slim.
- (Later in de cursus: kies matrixfactorisaties die stabiliteit geven.)

---

## Links naar andere hoofdstukken (de bruggen die examenvragen graag testen)

- **Hoofdstuk 2 (Linear systems)**  
  Pivoting en keuze van factorisatie (LU vs Gauss-Jordan vs Cholesky) zijn antwoorden op rounding error, stabiliteit en conditioning.

- **Hoofdstuk 3 (Least squares)**  
  “Normal equations” vs “QR” vs “SVD” is een klassieker:
  - normal equations kunnen conditioning verslechteren $\to$ rounding wordt gevaarlijker,
  - QR en SVD zijn duurder maar stabieler.

- **Niet-lineaire oplossingen / optimalisatie / integratie / ODE’s**  
  Overal zie je hetzelfde thema: kleinere $h$ verlaagt truncation, maar kan rounding/foutopstapeling verhogen; stabiliteit bepaalt of “harder werken” ook echt “beter antwoord” geeft.

---

## Mini-checklist (mondeling examen waardig)
Als je moet verdedigen waarom je aanpak numeriek “goed” is, zeg expliciet:
1. Welke fouten bestaan hier? (truncation + rounding)
2. Waar zit cancellation/ill-conditioning?
3. Welke herformulering / datatype / somvolgorde maakt het stabieler?
4. Wat kost dat (tijd/geheugen) en waarom is die kost het waard?
