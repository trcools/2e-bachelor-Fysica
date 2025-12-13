# Hoofdstuk 4 — Eigenvalue Problems

In dit hoofdstuk draait alles rond het eigenprobleem
$$
\mathbf{A}\mathbf{x}=\lambda \mathbf{x},
$$
waar $\lambda$ een **eigenwaarde** is en $\mathbf{x}\neq \mathbf{0}$ een bijhorende **eigenvector**.

De examenkern is (zoals altijd in deze cursus): **welke methode kies je wanneer, en waarom (kost + stabiliteit + doel: 1 eigenwaarde of allemaal)?**  
En: hoe dit terugkoppelt naar H2 (lineaire stelsels oplossen) en H3 (SVD/conditioning/rank).

---

## 4.1 Introductie, concept en nuttige eigenschappen

### 4.1.1 Characteristic polynomial (theorie, maar numeriek meestal NIET doen)
Uit $\mathbf{A}\mathbf{x}=\lambda \mathbf{x}$ volgt
$$
(\mathbf{A}-\lambda \mathbf{I})\mathbf{x}=\mathbf{0}.
$$
Een niet-triviale oplossing bestaat enkel als
$$
\det(\mathbf{A}-\lambda \mathbf{I})=0.
$$
De polynoom
$$
p(\lambda)=\det(\mathbf{A}-\lambda \mathbf{I})
$$
heet de **characteristic polynomial**; zijn wortels zijn de eigenwaarden.

**Belangrijk (cursus zegt dit expliciet):** eigenwaarden vinden via de wortels van $p(\lambda)$ is *geen goede numerieke strategie* voor matrices van niet-triviale grootte, o.a. omdat
- de coëfficiënten van $p(\lambda)$ extreem gevoelig kunnen zijn voor kleine perturbaties in $\mathbf{A}$,
- rounding errors in $p(\lambda)$ de wortels volledig kunnen verpesten,
- wortels van een hoge-graads polynoom vinden zelf al een lastig numeriek probleem is.

**Conclusie:** in praktijk gebruik je iteratieve methodes (power/inverse/Rayleigh/QR) of library-routines.

---

### 4.1.2 Properties and transformations (dit stuurt de algoritmes)

De cursus geeft een lijst “eigenwaarden blijven hetzelfde, of transformeren voorspelbaar”:

- **Symmetric/Hermitian:**  
  Als $\mathbf{A}$ symmetric/Hermitian is, dan zijn alle eigenwaarden **reëel**.

- **Shift:** als $\mathbf{A}\mathbf{x}=\lambda\mathbf{x}$ en $\sigma$ is een scalair:
  $$
  (\mathbf{A}-\sigma\mathbf{I})\mathbf{x}=(\lambda-\sigma)\mathbf{x}.
  $$
  Eigenwaarden schuiven met $\sigma$, eigenvectoren blijven hetzelfde.

- **Inversion:** $\mathbf{A}^{-1}$ heeft dezelfde eigenvectoren, eigenwaarden worden
  $$
  \lambda \mapsto \frac{1}{\lambda}.
  $$

- **Powers:** $\mathbf{A}^k$ heeft dezelfde eigenvectoren, eigenwaarden worden
  $$
  \lambda \mapsto \lambda^k.
  $$

- **Polynomials:** voor een polynoom $p(t)$:
  $$
  p(\mathbf{A})\mathbf{x}=p(\lambda)\mathbf{x}.
  $$
  Dus eigenwaarden transformeren als $\lambda\mapsto p(\lambda)$, eigenvectoren blijven die van $\mathbf{A}$.

- **Similarity:** $\mathbf{B}$ is similar aan $\mathbf{A}$ als er een invertibele $\mathbf{T}$ bestaat zodat
  $$
  \mathbf{B}=\mathbf{T}^{-1}\mathbf{A}\mathbf{T}.
  $$
  Dan hebben $\mathbf{A}$ en $\mathbf{B}$ **dezelfde eigenwaarden**, en eigenvectoren worden systematisch “meegetransformeerd”.

**Waarom dit belangrijk is:** QR-iteratie en veel “matrix-reducties” zijn gebouwd op similarity-transformaties (eigenwaarden blijven identiek, maar de matrix wordt eenvoudiger).

---

## 4.2 Eigenwaarden en eigenvectoren berekenen (methodes)

### 4.2.1 Power iteration (dominante eigenwaarde)
Doel: schat de **dominante** eigenwaarde (grootste modulus) en eigenvector.

Idee:
- Kies een willekeurige $\mathbf{x}_0\neq \mathbf{0}$,
- herhaal
  $$
  \mathbf{x}_{k}=\mathbf{A}\mathbf{x}_{k-1},
  $$
- normaliseer elke stap (bv. met $\|\cdot\|_\infty$ of $\|\cdot\|_2$) om overflow/onderflow te vermijden.

**Waarom dit werkt (cursusproof in woorden):** schrijf
$$
\mathbf{x}_0=\sum_{j=1}^n \alpha_j \mathbf{v}_j,
$$
met $\mathbf{v}_j$ eigenvectoren. Dan
$$
\mathbf{x}_k=\mathbf{A}^k\mathbf{x}_0=\sum_{j=1}^n \alpha_j \lambda_j^k \mathbf{v}_j.
$$
Als er één unieke eigenwaarde $\lambda_1$ is met maximale modulus en $\alpha_1\neq 0$, dan domineert die term en convergeert de richting van $\mathbf{x}_k$ naar $\mathbf{v}_1$.

**Convergentiesnelheid (intuïtie):**
$$
\text{foutfactor} \sim \left|\frac{\lambda_2}{\lambda_1}\right|^k.
$$
Dus traag als $|\lambda_2|\approx|\lambda_1|$.

**Wanneer gebruiken?**
- Je wil alleen de grootste eigenwaarde/eigenvector.
- Matrix is groot en je wil iets heel simpels.
- Je aanvaardt lineaire (soms trage) convergentie.

**Wanneer niet?**
- Je wil meerdere eigenwaarden.
- Dominante eigenwaarde is niet uniek of spectrum ligt “dicht opeen”.

---

### 4.2.2 Inverse iteration (kleinste eigenwaarde, of eigenwaarde dicht bij een shift)
De cursus: inverse iteration convergeert naar de eigenvector van de **grootste eigenwaarde van $\mathbf{A}^{-1}$**, dus naar de eigenvector van de **kleinste** eigenwaarde van $\mathbf{A}$.

Praktisch doe je niet expliciet $\mathbf{A}^{-1}$, je lost per iteratie een lineair stelsel op:
$$
\mathbf{A}\mathbf{y}_k=\mathbf{x}_{k-1},\qquad \mathbf{x}_k=\frac{\mathbf{y}_k}{\|\mathbf{y}_k\|}.
$$

**Shifted inverse iteration (superbelangrijk):** wil je een eigenwaarde nabij $\sigma$? Gebruik
$$
(\mathbf{A}-\sigma\mathbf{I})\mathbf{y}_k=\mathbf{x}_{k-1}.
$$
Dan convergeer je naar de eigenvector van de eigenwaarde van $\mathbf{A}$ die het **dichtst bij $\sigma$** ligt.

**Opmerking uit de cursus:** bij shifted inverse iteration krijg je in feite de inverse van de **geshiftte** eigenwaarde; om terug naar de eigenwaarde van $\mathbf{A}$ te gaan:
- neem het reciproke,
- tel $\sigma$ terug erbij.

**Kost/implementatie-inzicht (link met H2):**
- Als $\sigma$ vast blijft: factoriseer $\mathbf{A}-\sigma\mathbf{I}$ één keer (LU) en hergebruik in elke iteratie.
- Dus “duur + veel goedkope solves” (zoals H2).

**Wanneer gebruiken?**
- Je wil een eigenwaarde/eigenvector **in een specifiek gebied** van het spectrum.
- Je hebt een redelijke shift (bv. uit fysisch inzicht, of uit een ruwe schatting).

---

### 4.2.3 Rayleigh quotient iteration (snelle convergentie met slimme shift)
De cursus linkt dit aan een LS-probleem:
$$
\mathbf{x}\,\lambda \cong \mathbf{A}\mathbf{x}.
$$
De beste LS-schatting van $\lambda$ is de **Rayleigh quotient**:
$$
\rho(\mathbf{x})=\frac{\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}}{\mathbf{x}^\mathsf{T}\mathbf{x}}
\quad
(\text{voor reële } \mathbf{A}, \mathbf{x}).
$$

Rayleigh quotient iteration gebruikt deze als shift:
1. $\sigma_k=\rho(\mathbf{x}_k)$
2. los
   $$
   (\mathbf{A}-\sigma_k\mathbf{I})\mathbf{y}_k=\mathbf{x}_k
   $$
3. normaliseer: $\mathbf{x}_{k+1}=\mathbf{y}_k/\|\mathbf{y}_k\|$

**Waarom dit “episch goed” kan zijn:**
- Als je al een redelijke eigenvector-guess hebt, gaat dit vaak veel sneller dan gewone power/inverse.
- (Klassieke extra kennis: voor symmetric matrices heb je vaak zeer snelle, zelfs kubische convergentie wanneer je dichtbij zit.)

**Wanneer gebruiken?**
- Je wil heel snel “finetunen” naar een eigenpaar.
- Je hebt al een goede startvector (bv. van inverse iteration of fysische mode-vorm).

---

### 4.2.4 Deflation (meerdere eigenwaarden na elkaar, maar met valkuilen)
Deflation probeert na het vinden van $(\lambda_1,\mathbf{x}_1)$ een nieuwe matrix te maken waarin die eigenwaarde “verwijderd” is.

De cursus toont de constructie met een vector $\mathbf{u}_1$ zodat
$$
\mathbf{u}_1^\mathsf{T}\mathbf{x}_1=\lambda_1,
$$
en dan
$$
\mathbf{A}_{\text{deflated}}=\mathbf{A}-\mathbf{x}_1\mathbf{u}_1^\mathsf{T}.
$$

Daarna kan je opnieuw power iteration doen om de “volgende” eigenwaarde te vinden.

**Cursus-waarschuwing:** deflation wordt snel
- omslachtig,
- numeriek minder accuraat,
- en in de praktijk gebruik je betere methodes om veel eigenwaarden te vinden.

**Wanneer toch nuttig?**
- Klein probleem, didactisch.
- Je wil “een paar” eigenwaarden en accepteert dat je later inverse iteration met shifts nodig hebt.

---

### 4.2.5 QR iteration (workhorse voor alle eigenwaarden)
De cursus noemt dit “in practice, the fastest and most used method” voor alle eigenwaarden (dense case).

QR-iteratie definieert een reeks:
1. QR-factorisatie:
   $$
   \mathbf{A}_m=\mathbf{Q}_m\mathbf{R}_m
   $$
2. Update:
   $$
   \mathbf{A}_{m+1}=\mathbf{R}_m\mathbf{Q}_m
   $$

**Cruciale eigenschap (waarom eigenwaarden behouden blijven):**
Omdat $\mathbf{Q}_m$ orthogonaal is,
$$
\mathbf{A}_{m+1}=\mathbf{R}_m\mathbf{Q}_m
=
\mathbf{Q}_m^\mathsf{T}\mathbf{A}_m\mathbf{Q}_m,
$$
dus $\mathbf{A}_{m+1}$ is **similar** aan $\mathbf{A}_m$ → eigenwaarden blijven dezelfde.

De iteratie convergeert naar een (quasi-)triangulaire vorm; de diagonaal convergeert naar de eigenwaarden.

**Shifts in QR iteration (versnellen):**
Neem bv. als shift de rechtsonder entry:
$$
\mu_m=(\mathbf{A}_m)_{nn}.
$$
Doe QR op $\mathbf{A}_m-\mu_m\mathbf{I}$ en zet daarna shift terug:
$$
\mathbf{A}_{m+1}=\mathbf{R}_m\mathbf{Q}_m+\mu_m\mathbf{I}.
$$
Cursus: dit laat off-diagonale elementen sneller naar 0 gaan.

**Wanneer gebruiken?**
- Je wil **alle** eigenwaarden (en eventueel eigenvectoren) van een dense matrix.
- Je wil een beproefde, robuuste methode (library).

---

## 4.3 Singular Value Decomposition berekenen (koppeling met eigenwaarden)
De cursus herhaalt:
$$
\mathbf{A}=\mathbf{U}\mathbf{\Sigma}\mathbf{V}^\mathsf{T},
$$
met $\mathbf{U}$ orthogonaal ($m\times m$), $\mathbf{V}$ orthogonaal ($n\times n$), en $\mathbf{\Sigma}$ diagonaal ($m\times n$) met $\sigma_i\ge 0$.

**Belangrijke link (ook in H3):**
- $\sigma_i^2$ zijn eigenwaarden van $\mathbf{A}^\mathsf{T}\mathbf{A}$,
- right singular vectors $\mathbf{v}_i$ zijn eigenvectoren van $\mathbf{A}^\mathsf{T}\mathbf{A}$,
- left singular vectors $\mathbf{u}_i$ zijn eigenvectoren van $\mathbf{A}\mathbf{A}^\mathsf{T}$,
en
$$
\mathbf{A}\mathbf{v}_i=\sigma_i \mathbf{u}_i.
$$

In het symmetrische voorbeeld uit de cursus geldt zelfs $\mathbf{U}=\mathbf{V}$.

**Wanneer denk je “SVD” i.p.v. “eigen”?**
- Niet-vierkante matrix.
- Rank/conditioning/least-squares diagnose (H3).
- Je wil robuustheid i.p.v. pure snelheid.

---

## 4.4 Software (SciPy) — wat kies je wanneer?

Cursusboodschap (belangrijk): QR-iteratie was lang de standaard; nieuwere methodes (divide-and-conquer, RRR) zijn vaak sneller voor **alle eigenvectoren**, maar QR heeft een lange “reliability track record”. In SciPy is `linalg.eig` de meest algemene methode.

Praktisch:
- `scipy.linalg.eig`  
  Algemeen eigenprobleem (complex mogelijk), gebruikt QR-achtige aanpak.
- Voor matrices met speciale structuur zijn er snellere routines (cursus geeft o.a. tridiagonal symmetric varianten).
- Voor SVD:
  - `scipy.linalg.svd`
  - `scipy.linalg.svdvals`

**Vuistregels**
- Symmetric/Hermitian? Gebruik de symmetric/Hermitian routines (sneller + reële eigenwaarden + stabieler).
- Tridiagonal symmetric? Gebruik tridiagonal routine (nog sneller).
- Algemeen? `eig`.

---

## 4.5 Physics Example — spring-and-mass system (waarom eigenwaarden fysisch tellen)

De cursus zet een gekoppeld veer-massa systeem om naar een matrixvorm.

Voor twee massa’s met veren krijg je (in tweede orde) iets van de vorm:
$$
\ddot{\mathbf{x}}=\mathbf{A}\mathbf{x},
$$
met $\mathbf{x}=[x_1,x_2]^\mathsf{T}$ en een matrix $\mathbf{A}$ die afhangt van $k$ en $M_1,M_2$ (in de cursus staat expliciet een $2\times 2$ matrix met termen zoals $-2k/M_1$, $k/M_1$, etc.).

De eigenwaarden $\lambda$ komen uit
$$
\det(\mathbf{A}-\lambda \mathbf{I})=0.
$$

**Fysische interpretatie (normale modes):**
- Voor een stabiel veer-systeem verwacht je oscillaties met $\omega$.
- Typisch geldt dan $\lambda=-\omega^2$ (teken hangt af van hoe $\mathbf{A}$ exact gedefinieerd is).
- Eigenvectoren geven de **mode-vormen** (hoe massa’s relatief bewegen).
- Eigenwaarden geven de **frequenties** (via $\omega=\sqrt{-\lambda}$ als $\lambda<0$).

Dit is precies waarom eigenproblemen zo vaak opduiken in fysica: *coupled linear dynamics = modes = eigenvectors*.

---

## “Welke methode wanneer?” — examengerichte beslisregels

### Je wil één dominante eigenwaarde/eigenvector
- **Power iteration**
- Waarom: goedkoop per iteratie (matrix-vector product), eenvoudig.
- Valkuil: traag als spectrum dicht is; faalt als dominante eigenwaarde niet uniek is.

### Je wil de kleinste eigenwaarde of eentje dicht bij een schatting $\sigma$
- **(Shifted) inverse iteration**
- Waarom: door $(\mathbf{A}-\sigma\mathbf{I})^{-1}$ wordt “dichtst bij $\sigma$” dominant.
- Kost: per iteratie een solve; met vaste shift kan je LU hergebruiken.

### Je hebt al een goede eigenvector-guess en wil razendsnel convergeren
- **Rayleigh quotient iteration**
- Waarom: shift wordt automatisch “best mogelijke” eigenwaarde-schatting voor de huidige vector.

### Je wil meerdere eigenwaarden, maar het probleem is klein/onderwijs
- **Deflation** (met gezond wantrouwen)
- Waarom: conceptueel oké, maar wordt numeriek snel minder mooi.

### Je wil alle eigenwaarden (dense matrix, standaard library-werk)
- **QR iteration / `scipy.linalg.eig`**
- Waarom: robuust, standaard, convergent naar (quasi-)triangulaire vorm waarvan de diagonaal de eigenwaarden geeft.
- Gebruik shifts voor versnelling (bibliotheken doen dit slim).

---

## Links met andere hoofdstukken (type-1 examenvragen)

- **Link met H2 (lineaire stelsels):** inverse/Rayleigh iteraties bestaan uit “herhaald een lineair stelsel oplossen” → LU/Cholesky keuzes bepalen runtime en stabiliteit.
- **Link met H3 (least squares & SVD):** singular values via eigenwaarden van $\mathbf{A}^\mathsf{T}\mathbf{A}$; rank/conditioning komt rechtstreeks uit SVD.
- **Link met H1 (numerical limitations):** “characteristic polynomial roots” is het poster-child van een wiskundig correcte maar numeriek slechte route.

---
