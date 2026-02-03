<!-- file: hoofdstuk_11_fast_fourier_transform.md -->
# Hoofdstuk 11 — Fast Fourier Transform (FFT)

FFT is een snelle manier om de DFT te berekenen. Het is niet alleen “spectra plotten”; het is een **reken-truc** om dingen te versnellen die anders $O(N^2)$ kosten (convolutie/correlatie) en een **analyse-tool** (frequenties, filtering, aliasing).

Examenkern:
- DFT/IDFT + interpretatie
- FFT-kost: $O(N\log N)$
- spectral leakage + windows
- convolution theorem + zero-padding (circulair vs lineair!)
- DFT vs DCT: wanneer welke

---

## 11.1 DFT: definitie en interpretatie

Voor samples $x_n$, $n=0,\dots,N-1$:

DFT:
$$
X_k=\sum_{n=0}^{N-1} x_n\,e^{-2\pi i kn/N}.
$$

Inverse:
$$
x_n=\frac{1}{N}\sum_{k=0}^{N-1} X_k\,e^{2\pi i kn/N}.
$$

### Frequentieschaling
Als sample rate $f_s$ (samples per seconde), dan correspondeert bin $k$ met frequentie:
$$
f_k=\frac{k}{N}f_s
$$
(en bij real signals is er symmetrie zodat je vaak `rfft` gebruikt).

### Periodiciteitsassumptie
DFT behandelt je data als één periode van een periodiek signaal met periode $N$.
- Als $x_0\neq x_{N-1}$ in “trend”, creëer je een rand-discontinuïteit.
- Dat veroorzaakt leakage (zie 11.4).

---

## 11.2 FFT: waarom het snel is

Direct DFT: $O(N^2)$.

FFT (Cooley–Tukey): splits in even/odd indices en hergebruik kleinere DFT’s:
- $O(N\log N)$.

**Praktische impact:**
- $N=10^6$ is haalbaar met FFT, onmogelijk met directe DFT.

---

## 11.3 Real FFT en efficiëntie

Voor reële $x_n$ geldt conjugate symmetry:
$$
X_{N-k}=\overline{X_k}.
$$
Dus je kan de helft opslaan/berekenen:
- `rfft` en `irfft` zijn sneller en zuiniger.

---

## 11.4 Spectral leakage, resolution, windows

### 11.4.1 Wat is leakage?
Als je een zuivere sinus hebt met frequentie die niet exact op een bin valt, dan “lekt” de energie naar andere bins.

Intuïtie: je neemt een eindig venster van een oneindig signaal. Dat is een vermenigvuldiging met een window in tijd → in frequentie is dat een convolutie met de Fourier transform van de window. Dus je spectrum wordt uitgesmeerd.

### 11.4.2 Windowing (Hann/Hamming/Blackman/...)
Je vermenigvuldigt $x_n$ met een window $w_n$ die naar 0 gaat aan de randen.

- Pro: veel minder leakage (side lobes omlaag)
- Con:
  - main lobe breder → slechtere frequentieresolutie
  - amplitude-bias (je moet vaak “window correction” doen als amplitudes exact moeten)

**Wanneer windowen?**
- Bijna altijd bij real-world data die niet exact periodiek is in het sample-interval.
- Als je vooral frequenties wil detecteren i.p.v. perfecte amplitudes.

**Wanneer niet?**
- Als je signaal exact periodic in het interval is (zeldzaam), dan is boxcar oké en geeft hoogste resolutie.

---

## 11.5 Aliasing (Nyquist)

Als je sample rate $f_s$ is, dan is de Nyquist frequentie:
$$
f_N=\frac{f_s}{2}.
$$
Frequenties boven $f_N$ folden terug (aliasing). FFT kan dat niet “fixen”; dat is een sampling probleem.

**Examenvriendelijk inzicht:** aliasing is geen numerieke fout van FFT; het is informatieverlies door sampling.

---

## 11.6 Convolution en correlation met FFT (de grote speed-up)

Discrete lineaire convolutie:
$$
(y*g)_n=\sum_m y_m g_{n-m}.
$$

DFT-convolution theorem:
$$
\mathrm{DFT}(y*g)=\mathrm{DFT}(y)\,\mathrm{DFT}(g).
$$

### Circulair vs lineair (de valkuil)
FFT op lengte $N$ geeft **circulaire convolutie**:
- indices wrappen modulo $N$.

Voor lineaire convolutie moet je **zero-padding** doen tot minstens:
$$
N_{\text{pad}} \ge N_y + N_g - 1.
$$

**Typische “insights”-vraag na coding:**
- Waarom zag ik wrap-around/artefacten?  
  → omdat je circulair deed zonder padding.

### Correlatie
Correlatie lijkt op convolutie maar met een omkering/conjugatie:
$$
(r_{xy})_k = \sum_n x_n \overline{y_{n-k}}.
$$
FFT-truc blijft: FFT’s + puntgewijs product + inverse FFT.

---

## 11.7 DFT vs DCT (wanneer welke, zoals jullie vaak vragen)

### DFT
- basisfuncties: complexe exponenten
- impliciet periodiek: “wrap-around”

### DCT (Discrete Cosine Transform)
DCT komt overeen met een **even (reflective) extension** van je data:
- minder rand-discontinuïteit als je signaal niet-periodiek is,
- vaak compacter spectrum voor gladde niet-periodieke functies.

**Wanneer is DCT beter?**
- data op $[0,L]$ met “reflectie”-achtig randgedrag,
- compressie (energie in lage modes),
- PDE’s met Neumann-achtige BC’s (cosinus-basis natuurlijk).

**Wanneer DFT beter?**
- echt periodieke signals,
- complexe fase-informatie,
- convolution/correlation in periodieke settings.

---

## 11.8 Welke methode wanneer? (beslisregels)

- **Spectrale analyse (frequenties vinden):** FFT (`rfft` voor real).
- **Filtering / smoothing / band-pass:** FFT + filtermasker (maar let op windowing & edge effects).
- **Snelle convolutie met lange kernel:** FFT + zero-padding.
- **Niet-periodieke, gladde data op een interval:** DCT vaak beter dan DFT.
- **Als randvoorwaarden periodic zijn (PDE/spectraal):** DFT/FFT is natuurlijk.

---

## 11.9 Typische valkuilen (die punten kosten)
1) Geen zero-padding → circulaire artefacten.
2) Geen window → leakage en “spookfrequenties”.
3) Verkeerde frequentie-as (fftfreq) → verkeerde interpretatie.
4) Aliasing verwarren met leakage.
5) Vergeten dat FFT output complex is: magnitude $|X_k|$, phase $\arg(X_k)$.

---

## SciPy mapping
- `scipy.fft.fft`, `ifft`
- `scipy.fft.rfft`, `irfft`
- `scipy.fft.fftfreq`, `rfftfreq`
- `scipy.fft.next_fast_len` (snelle padding lengte)
- `scipy.fft.dct`, `idct` (voor DCT)

---

## Links met andere hoofdstukken
- H7/H8: interpolatie/differentiatie kunnen spectraal (FFT) bij periodieke problemen.
- H10: spectrale methodes en snelle Poisson-solvers bij periodieke BC.
- H1: interpretatieproblemen (leakage/aliasing) domineren; numerieke rounding is zelden de bottleneck.
