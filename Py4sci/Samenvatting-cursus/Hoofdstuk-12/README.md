<!-- file: hoofdstuk_12_monte_carlo.md -->
# Hoofdstuk 12 — Monte Carlo

Monte Carlo (MC) methodes lossen numerieke problemen op via steekproeven. Het “superpower”-argument:
- de fout schaalt typisch als $1/\sqrt{N}$, bijna onafhankelijk van dimensie,
- dus bij hoge dimensie (waar deterministische quadrature explodeert) is MC vaak de enige praktische optie.

Tegelijk: $1/\sqrt{N}$ is traag. Het examen draait dus om:
1) **wanneer MC de juiste hamer is**,  
2) **hoe je error en betrouwbaarheid kwantificeert**,  
3) **hoe je variance verlaagt** (variance reduction),  
4) **wat MCMC doet en wanneer je het nodig hebt**.

---

## 12.1 Randomness in computation: PRNG en reproducibility

### Pseudo-random number generators (PRNG)
Computers genereren geen “echte” random getallen, maar pseudo-random:
- deterministisch algoritme,
- lijkt statistisch random (als de generator goed is),
- volledig bepaald door een **seed**.

**Waarom dit exam-relevant is:** reproduceerbaarheid.
- Zelfde seed $\Rightarrow$ zelfde resultaten (essentieel in wetenschappelijke code).
- Andere seed $\Rightarrow$ statistisch equivalente resultaten (als de generator goed is).

### Praktische good practice (zoals in moderne NumPy/SciPy)
- Maak één generator-object.
- Vermijd “global state” (zoals `np.random.seed(...)` overal).

---

## 12.2 Monte Carlo integratie: de basis-estimator

We willen een integraal over een domein $D$:
$$
I=\int_D f(\mathbf{x})\,d\mathbf{x}.
$$

### 12.2.1 Uniform sampling op een domein
Als je uniform $\mathbf{x}_i$ in $D$ trekt en $V=\mathrm{vol}(D)$:
$$
\hat{I}_N = V \,\frac{1}{N}\sum_{i=1}^N f(\mathbf{x}_i).
$$

**Unbiasedness (intuïtief):**
$$
\mathbb{E}[\hat{I}_N]=I
$$
als je echt uniform samplet.

### 12.2.2 Foutschatting (standard error)
Definieer de steekproefgemiddelde:
$$
\bar{f}=\frac{1}{N}\sum_{i=1}^N f(\mathbf{x}_i),
$$
en steekproefvariantie:
$$
s_f^2 = \frac{1}{N-1}\sum_{i=1}^N (f(\mathbf{x}_i)-\bar{f})^2.
$$
Dan is de standaardfout op $I$:
$$
\mathrm{SE}(\hat{I}_N) \approx V\,\frac{s_f}{\sqrt{N}}.
$$

**Cruciale interpretatie:**
- Als je error 10× kleiner wil, heb je ~100× meer samples nodig.
- Daardoor is variance reduction vaak belangrijker dan “gewoon meer samples”.

### 12.2.3 Confidence intervals (CI)
Voor groot $N$ (CLT) is:
$$
\hat{I}_N \approx \mathcal{N}(I,\mathrm{SE}^2).
$$
Een (ongeveer) 95% interval:
$$
\hat{I}_N \pm 1.96\,\mathrm{SE}.
$$

**Wat je mondeling moet kunnen zeggen:**
- Dit is asymptotisch (werkt beter voor grote $N$),
- Bij heavy tails of kleine $N$ kan dit misleidend zijn.

---

## 12.3 Curse of dimensionality: waarom MC soms wint

Deterministische quadrature in $d$ dimensies:
- nodes groeien vaak exponentieel met $d$.

MC:
- error $\sim 1/\sqrt{N}$, niet exponentieel afhankelijk van $d$.

**Beslisregel (exam-waardig):**
- $d\le 3$ en $f$ glad: deterministisch (H8) wint vaak.
- $d\gg 1$: MC is vaak de enige haalbare.

---

## 12.4 Variance reduction: sneller dan $1/\sqrt{N}$ in de praktijk

Het doel is niet “andere asymptotiek” (meestal blijft $1/\sqrt{N}$), maar een veel kleinere constante: $\mathrm{Var}$ omlaag.

### 12.4.1 Importance sampling
Kies een pdf $p(\mathbf{x})$ die lijkt op waar $|f|$ groot is.
Dan:
$$
I=\int_D f(\mathbf{x})\,d\mathbf{x}
=\int_D \frac{f(\mathbf{x})}{p(\mathbf{x})} p(\mathbf{x})\,d\mathbf{x}
\approx \frac{1}{N}\sum_{i=1}^N \frac{f(\mathbf{x}_i)}{p(\mathbf{x}_i)},
\quad \mathbf{x}_i\sim p.
$$

**Wanneer gebruiken?**
- als $f$ “spiky” is of lokaal groot (anders mis je de piek met uniform sampling),
- bij zeldzame events.

**Valkuilen:**
- als $p$ te klein is waar $f$ groot is, krijg je enorme gewichten $\frac{f}{p}$ → variance explodeert.
- je moet $p$ kunnen samplen én evalueren.

---

### 12.4.2 Stratified sampling
Splits $D$ in strata $D_j$ met volumes $V_j$.
Sample in elk stratum apart en combineer:
$$
\hat{I}=\sum_j V_j \bar{f}_j.
$$

**Wanneer gebruiken?**
- als $f$ sterk varieert per regio,
- als je “garandeert” dat elke regio vertegenwoordigd is.

---

### 12.4.3 Control variates
Vind $g$ met bekende integraal $G=\int g$ en sterke correlatie met $f$.
Neem estimator:
$$
\hat{I} = \hat{I}_f - c(\hat{I}_g - G).
$$
Kies $c$ optimaal via covariantie/variantie (in praktijk vaak geschat uit samples).

**Wanneer gebruiken?**
- als je een benaderende “goed integrerende” functie kent die op $f$ lijkt.

---

### 12.4.4 Antithetic variates (symmetrie benutten)
Sample in paren $(\mathbf{x},\tilde{\mathbf{x}})$ die negatief correleren (bv. $u$ en $1-u$ bij uniforme variabelen).
- vermindert variance als $f$ monotone/structureel is.

---

## 12.5 Monte Carlo voor andere taken

### 12.5.1 Monte Carlo voor $\pi$ / geometrische probabiliteit
Klassiek: area/volume ratio’s.
Dit is vooral didactisch om unbiasedness en $1/\sqrt{N}$ te tonen.

### 12.5.2 Random walks en diffusion
Random walks linken MC aan fysische processen (diffusie, Brownse beweging).

---

## 12.6 Markov Chain Monte Carlo (MCMC): sampelen uit moeilijke verdelingen

Soms wil je geen integraal rechtstreeks, maar samples uit een verdeling $\pi(\mathbf{x})$ (bv. posterior in Bayes).

Als direct sampelen niet lukt, gebruik je een Markov chain met stationaire verdeling $\pi$.

### 12.6.1 Metropolis–Hastings
- huidige state: $\mathbf{x}$
- voorstel (proposal): $\mathbf{x}'\sim q(\mathbf{x}'|\mathbf{x})$
- acceptatie:
$$
\alpha(\mathbf{x}\to \mathbf{x}')=
\min\left(
1,
\frac{\pi(\mathbf{x}')\,q(\mathbf{x}|\mathbf{x}')}{\pi(\mathbf{x})\,q(\mathbf{x}'|\mathbf{x})}
\right).
$$

Voor symmetrische proposal (random walk):
$$
\alpha=
\min\left(1,\frac{\pi(\mathbf{x}')}{\pi(\mathbf{x})}\right).
$$

### 12.6.2 Burn-in, mixing, autocorrelation
MCMC-samples zijn **gecorreleerd**.
- burn-in: eerste stuk weggooien (chain nog niet in stationaire regime),
- autocorrelatie betekent dat “effectieve N” kleiner is dan het aantal stappen,
- tuning van proposal-stepsize beïnvloedt acceptance rate en mixing.

**Examengerichte insight:**
- Te kleine stappen: hoge acceptatie maar langzaam exploreren.
- Te grote stappen: lage acceptatie, chain “stuck”.
Je wil een “sweet spot” (afhankelijk van dimensie/proposal).

---

## 12.7 Welke methode wanneer?

### Gebruik deterministische quadrature (H8) als…
- dimensie laag is,
- $f$ glad is,
- je hoge precisie wil met weinig evaluaties.

### Gebruik plain Monte Carlo als…
- dimensie hoog is,
- domein complex is,
- je snel een ruwe schatting + foutbalk wil.

### Gebruik variance reduction als…
- plain MC te traag convergeert,
- $f$ spiky / rare-event / sterk variabel is,
- je structuur kunt exploiteren (symmetrie, bekende vergelijkbare $g$).

### Gebruik MCMC als…
- je samples uit $\pi(\mathbf{x})$ nodig hebt en direct sampling moeilijk is,
- je Bayesiaanse inferentie doet,
- je integralen wil als expectations onder $\pi$:
  $$
  \mathbb{E}_\pi[h(\mathbf{x})]\approx \frac{1}{N}\sum_{i=1}^N h(\mathbf{x}_i),
  \quad \mathbf{x}_i\sim \pi.
  $$

---

## Links met andere hoofdstukken

- H8: multi-D integratie → MC is het natuurlijke alternatief voor curse of dimensionality.
- H1: errorbars, rounding, reproducibility (seed) zijn essentieel.
- H3/H6: stochastic sampling vs fitting/optimisation; control variates lijken op “modelfout compenseren”.
- Statistische interpretatie: CLT, variantie, confidence intervals (wordt vaak in “insights” gevraagd na coding).
