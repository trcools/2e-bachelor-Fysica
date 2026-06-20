# Samenvatting: Lie-groepen en Lie-algebra's van rotaties

*Gebaseerd op Groups and Representations, hoofdstuk 4*

---

## 0. Grote lijn — wat hoort bij wat?

```
            isomorf              Lie algebra            complexificatie
  SO(2)  ────────────  U(1)  ───────────────  (triviaal, abels: geen generatoren nodig)

  SO(3)  ──/≄/──  SU(2)         so(3)  ≃  su(2)  ───────────────  sl(2)
   (NIET isomorf als groep, WEL zelfde Lie algebra lokaal)
```

**De kernboodschap van dit hoofdstuk:** een Lie-groep bestaat uit twee lagen informatie:

|                                                       | Wat het beschrijft                                                    | Bepaalt                            |
| ----------------------------------------------------- | --------------------------------------------------------------------- | ---------------------------------- |
| **Lie-algebra** (so(3), su(2), sl(2))           | infinitesimale (kleine) transformaties dicht bij het neutraal element | *lokale* structuur               |
| **Lie-groep** (SO(2), SO(3), U(1), U(2), SU(2)) | alle (ook "grote") transformaties                                     | *globale/topologische* structuur |

Twee groepen kunnen dezelfde Lie-algebra hebben (so(3) ≃ su(2)) zonder zelf isomorf te zijn (SO(3) ≇ SU(2)). Dat is precies wat SO(3) en SU(2) interessant maakt.

---

## 1. SO(2) en U(1) — het abelse opwarmertje

### Definities

$$
SO(2) = \left\{ R(\theta) = \begin{pmatrix}\cos\theta & -\sin\theta \\ \sin\theta & \cos\theta\end{pmatrix} \;\middle|\; \theta \in [0,2\pi) \right\}, \qquad U(1) = \{e^{i\theta} \mid \theta \in [0,2\pi)\}
$$

- $SO(2)$: reële $2\times2$ orthogonale matrices met determinant $+1$ (rotaties in het vlak).
- $U(1)$: complexe getallen met norm 1 (fasefactoren).

### Isomorfisme

$$
SO(2) \cong U(1), \qquad R(\theta) \;\longleftrightarrow\; e^{i\theta}
$$

Beide groepen zijn **abels** (commutatief) en **1-dimensionaal** (parametrisatie door één hoek $\theta$).

### Irreps

Schur's lemma $\Rightarrow$ irreps van een abelse groep zijn altijd **1-dimensionaal**. Elk irrep is een homomorfisme $\varphi: U(1) \to U(1)$:

$$
\boxed{\varphi_n(e^{i\theta}) = e^{in\theta}, \qquad n \in \mathbb{Z}}
$$

- **Oneindig veel** irreps, gelabeld door een **geheel getal** $n$ (in tegenstelling tot eindige groepen: daar is het aantal irreps = aantal conjugatieklassen).
- Fysische naam: $n$ = (elektrische) lading; hangt samen met $U(1)$-ijktransformaties en Noether's theorema.
- De rotatiematrix $R(\theta)$ zelf is een **reducibele** 2D-representatie: $R(\theta) \sim \varphi_{+1} \oplus \varphi_{-1} = \mathrm{diag}(e^{i\theta}, e^{-i\theta})$.

Er bestaat **geen** Lie-algebra-discussie nodig voor $SO(2)/U(1)$: de groep is al abels, dus de "infinitesimale generator" is gewoon het getal $i$ (1-dimensionaal, triviaal commuterend).

---

## 2. SO(3) en SU(2) — de hoofdmoot

### 2.1 Definities

$$
SO(3) = \{M \in \mathbb{R}^{3\times3} \mid MM^T = M^TM = \mathbb{1},\ \det M = 1\}
$$

$$
SU(2) = \{U \in \mathbb{C}^{2\times2} \mid UU^\dagger = U^\dagger U = \mathbb{1},\ \det U = 1\}
$$

|                           | $SO(3)$                                         | $SU(2)$                                             |
| ------------------------- | ------------------------------------------------- | ----------------------------------------------------- |
| Matrices                  | reëel$3\times3$                                | complex$2\times2$                                   |
| Conditie                  | orthogonaal,$\det=1$                            | unitair,$\det=1$                                    |
| Bewaart                   | inproduct op$\mathbb{R}^3$                      | inproduct op$\mathbb{C}^2$                          |
| Dimensie (als variëteit) | 3                                                 | 3                                                     |
| Topologie                 | $\mathbb{RP}^3$ (niet enkelvoudig samenhangend) | $S^3$ (3-sfeer, **enkelvoudig samenhangend**) |

### 2.2 Parametrisaties

**$SO(3)$ — drie equivalente manieren:**

1. **Opeenvolgende rotaties rond assen:** $R(\theta_x,\theta_y,\theta_z) = R_z(\theta_z)R_y(\theta_y)R_x(\theta_x)$ — let op: volgorde is belangrijk, rotaties commuteren niet!
2. **As-hoek parametrisatie:** elke rotatie = rotatie met hoek $\theta \in [0,\pi]$ rond eenheidsvector $\vec n$, genoteerd $R_{\vec n}(\theta)$. Rechterhandregel.
   - $R_{\vec n}(\pi) = R_{-\vec n}(\pi)$: dit randgeval is cruciaal voor de topologie (§2.3)!
3. **Euler-hoeken** $(\psi,\theta,\varphi)$: rotatie rond $z$, dan (nieuwe) $x$, dan (nieuwe) $z$.

**$SU(2)$ — twee equivalente manieren:**

1. **Direct via complexe parameters:**

$$
U = \begin{pmatrix}\alpha & \beta \\ -\bar\beta & \bar\alpha\end{pmatrix}, \qquad \alpha,\beta\in\mathbb{C},\ |\alpha|^2+|\beta|^2=1
$$


$$
= \alpha_1\mathbb{1} + i\beta_2\sigma_X + i\beta_1\sigma_Y + i\alpha_2\sigma_Z
$$


2. **Via matrixexponentiaal (Pauli-matrices):**

$$
R_{\vec n}(\theta) = e^{i\theta\sigma_{\vec n}} = \cos\!\left(\frac{\theta}{2}\right)\mathbb{1} + 2i\sin\!\left(\frac{\theta}{2}\right)\sigma_{\vec n}, \qquad \theta \in [0,2\pi)
$$

   met $\sigma_{\vec n} = n_x t_x + n_y t_y + n_z t_z$, $t_i := \sigma_i/2$.

**Let op de gelijkenis met $SO(3)$'s as-hoek vorm** — dit is *geen* toeval (zie §2.3), maar de hoekranges verschillen: $SO(3)$ heeft $\theta\in[0,\pi]$, $SU(2)$ heeft $\theta\in[0,2\pi)$!

### 2.3 Topologie: waarom $SO(3) \ncong SU(2)$

- $SU(2)$ via $(\alpha_1,\alpha_2,\beta_1,\beta_2)$ met $\alpha_1^2+\alpha_2^2+\beta_1^2+\beta_2^2=1$ ligt op de **3-sfeer $S^3$** $\subset \mathbb{R}^4$: compact én **enkelvoudig samenhangend**.
- Equivalent: $SU(2)$ = bal $B(2\pi)$ van straal $2\pi$, waarbij elk punt op de rand ($r=2\pi$) wordt geïdentificeerd met **hetzelfde** element $-\mathbb{1}$.
- $SO(3)$ is ook een bal, maar van straal $\pi$, met **antipodale** punten op de rand geïdentificeerd ($R_{\vec n}(\pi) = R_{-\vec n}(\pi)$) $\Rightarrow$ niet enkelvoudig samenhangend.
- Gevolg: $SU(2) \to SO(3)$ is een **2-op-1 dekkende afbeelding** ($\pm U \mapsto$ dezelfde rotatie). $SU(2)$ is de *universele dekgroep* van $SO(3)$.


![](GROEPEN-EN-REPRESENTATIES/Figuren/topologie_SO3_vs_SU2.png)


**Vuistregel:** lokaal (dicht bij identiteit) zijn $SO(3)$ en $SU(2)$ niet te onderscheiden; vandaar dat hun Lie-algebra's identiek zijn. Pas op globale schaal (grote rotaties) zie je het verschil.

---

## 3. De Lie-algebra's: so(3), su(2) en sl(2)

### 3.1 Generatoren via infinitesimale rotaties

**so(3):** Taylor-ontwikkel $R_x(\delta\theta)$ rond $\delta\theta=0$:

$$
R_x(\delta\theta) \approx \mathbb{1}_3 + i\delta\theta\, L_x, \qquad L_x = \begin{pmatrix}0&0&0\\0&0&-i\\0&i&0\end{pmatrix}, \text{ en analoog } L_y, L_z
$$

**su(2):** Taylor-ontwikkel $R_{\vec n}(\delta\theta) = e^{i\delta\theta\sigma_{\vec n}}$:

$$
R_{\vec n}(\delta\theta) \approx \mathbb{1}_2 + i\delta\theta\left(n_x\frac{\sigma_X}{2} + n_y\frac{\sigma_Y}{2} + n_z\frac{\sigma_Z}{2}\right), \qquad L_i = \frac{\sigma_i}{2}
$$

### 3.2 De gedeelde commutatierelaties

Beide voldoen aan **dezelfde** commutatierelaties (met $i,j,k \equiv x,y,z \equiv 1,2,3$):

$$
\boxed{[L_i, L_j] = i\sum_k \varepsilon_{ijk} L_k}
$$

Dit is precies waarom men zegt $\mathfrak{so}(3) \cong \mathfrak{su}(2)$ — als abstracte Lie-algebra's zijn ze identiek; enkel de matrixrepresentatie verschilt ($3\times3$ Hermitisch vs. $2\times2$ Hermitisch).

|                    | $\mathfrak{so}(3)$                                    | $\mathfrak{su}(2)$                |
| ------------------ | ------------------------------------------------------- | ----------------------------------- |
| Generatoren        | $L_x,L_y,L_z$ ($3\times3$)                          | $\sigma_i/2$ ($2\times2$ Pauli) |
| Type matrices      | reële vectorruimte, complexe (Hermitische) generatoren | idem                                |
| Commutatierelaties | $[L_i,L_j]=i\varepsilon_{ijk}L_k$                     | identiek                            |
| Reële dimensie    | 3                                                       | 3                                   |

> **Subtiliteit:** so(3) en su(2) zijn *reële* vectorruimten (je neemt enkel reële lineaire combinaties van de generatoren), ook al zijn de generatoren zelf complex/Hermitisch.

### 3.3 Complexificatie: sl(2)

Om de irreps te vinden, definieer je **ladderoperatoren**:

$$
L_\pm := L_1 \pm iL_2
$$

Deze zijn **complexe** combinaties — ze verlaten de reële Lie-algebra su(2) en leven in de **complexificatie** sl(2):

$$
\mathfrak{sl}(2) \cong \mathfrak{su}(2) \oplus i\,\mathfrak{su}(2) \cong \mathfrak{su}(2) \otimes_{\mathbb{R}} \mathbb{C}
$$

**Definitie (sl(2)):** complexe 3-dimensionale Lie-algebra met generatoren $L_0, L_+, L_-$ en

$$
[L_0, L_\pm] = \pm 2L_\pm, \qquad [L_+, L_-] = L_0
$$

|             | so(3)                                  | su(2)                                  | sl(2)                           |
| ----------- | -------------------------------------- | -------------------------------------- | ------------------------------- |
| Generatoren | $L_1,L_2,L_3$                        | $L_1,L_2,L_3$ (Pauli/2)              | $L_0, L_+, L_-$               |
| Veld        | reëel                                 | reëel                                 | **complex**               |
| Relatie     | $\cong$ su(2)                        | $\cong$ so(3)                        | complexificatie van su(2)       |
| Notatie     | $\mathfrak{so}(3,\mathbb{R})$ achtig | $\mathfrak{sl}(2,\mathbb{R})$-achtig | $\mathfrak{sl}(2,\mathbb{C})$ |

**Waarom dit ertoe doet:** je rekent in de praktijk altijd met sl(2) (de ladderoperatoren), maar er bestaat een **1-op-1 correspondentie** tussen irreps van sl(2) en irreps van su(2) — omdat su(2) de *compacte reële vorm* van sl(2) is. Je vindt dus alle fysische (su(2)-)irreps door eerst sl(2) op te lossen.

---

## 4. Irreps van sl(2) / su(2) — highest-weight constructie

### 4.1 De Casimir-operator

$$
\vec L^2 := L_xL_x + L_yL_y + L_zL_z
$$

commuteert met **alle** generatoren: $[\vec L^2, L_i] = 0$. Volgens Schur's lemma is $\vec L^2$ dus een schalair op elk irrep — dit levert het label $l$.

### 4.2 Constructie in 4 stappen

Gemeenschappelijke eigenbasis $|l,m\rangle$ van $\vec L^2$ (eigenwaarde $l(l+1)$) en $L_3$ (eigenwaarde $m$):

1. **Start** bij het hoogste gewicht $|l,M\rangle$, geannihileerd door $L_+$: $L_+|l,M\rangle = 0$.
2. **Werk af** met $L_-$: genereert de keten $|l,M\rangle, |l,M-1\rangle, \dots$
3. $l$ is **constant** binnen één irrep (Casimir-eigenwaarde verandert niet onder $L_\pm$).
4. **Eindigheid** vereist een laagste gewicht $M'$ met $L_-|l,M'\rangle=0$.

Uitwerken van $L_-L_+|l,M\rangle=0$ en $L_+L_-|l,M'\rangle=0$ geeft:

$$
M = l, \qquad M' = -l \qquad \Rightarrow \qquad \boxed{\dim V_l = M - M' + 1 = 2l+1}
$$

### 4.3 Resultaat

$$
\boxed{\text{Irreps van } \mathfrak{sl}(2)/\mathfrak{su}(2): \quad l = 0, \tfrac12, 1, \tfrac32, 2, \dots, \qquad \dim = 2l+1}
$$

gespannen door $\{|l,m\rangle : m = l, l-1, \dots, -l\}$. Matrixelementen van de ladderoperatoren:

$$
L_\pm|l,m\rangle = N_\pm |l,m\pm1\rangle, \qquad N_\pm = \sqrt{l(l+1)-m(m\pm1)}
$$

**Fysische naamgeving:**

| Dimensie | $l$        | $l(l+1)$   | Naam                           |
| -------- | ------------ | ------------ | ------------------------------ |
| 1        | 0            | 0            | triviaal, scalair (singlet)    |
| 2        | $\tfrac12$ | $\tfrac34$ | fundamenteel, spinor (doublet) |
| 3        | 1            | 2            | adjoint, vector (triplet)      |

### 4.4 Het cruciale verschil tussen SU(2) en SO(3): integer vs. half-integer $l$

Dit is de **belangrijkste consequentie** van de topologische discussie in §2.3:

| $l$                                                | Irrep van su(2)/sl(2)? | Irrep van$SU(2)$ (groep)? | Irrep van $SO(3)$ (groep)?                  |
| ---------------------------------------------------- | ---------------------- | --------------------------- | ---------------------------------------------- |
| $0,1,2,\dots$ (**integer**)                  | ✓                     | ✓                          | ✓ (gewone representatie)                      |
| $\tfrac12,\tfrac32,\dots$ (**half-integer**) | ✓                     | ✓                          | ✗ — enkel**projectieve** representatie |

- **Alle** irreps van su(2) tillen op tot echte irreps van $SU(2)$ (de groep is enkelvoudig samenhangend $\Rightarrow$ exponentiële afbeelding is bijectief op irrep-niveau).
- Voor $SO(3)$ geldt dat enkel **integer** $l$ een echte representatie geeft. Half-integer $l$ (spinoren!) zijn alleen welgedefinieerd op de dekgroep $SU(2)$, en geven op $SO(3)$ slechts representaties op tekens na ($R_{\vec n}(2\pi) \mapsto -\mathbb{1}$ i.p.v. $+\mathbb{1}$) — vandaar "projectief".

Dit is precies de wiskundige reden waarom elektronen (spin $-\frac12$) een rotatie van $2\pi$ nodig hebben om met een extra factor $-1$ terug te keren, en pas na $4\pi$ écht naar de identiteit terugkeren.

---

## 5. Een overzichtstabel om alles samen te vatten

|                          | $SO(2)$                                                         | $U(1)$                     | $SO(3)$                                               | $SU(2)$                                                   |
| ------------------------ | ----------------------------------------------------------------- | ---------------------------- | ------------------------------------------------------- | ----------------------------------------------------------- |
| **Type matrices**  | reëel$2\times2$,<br />orthogonaal,<br /> $\det=1$            | $|z|=1 \subset \mathbb{C}$ | reëel$3\times3$,<br /> orthogonaal, <br />$\det=1$ | complex$2\times2$, <br />unitair, <br />$\det=1$        |
| **Abels?**         | ja                                                                | ja                           | nee                                                     | nee                                                         |
| **Parametrisatie** | $\theta\in[0,2\pi)$                                             | $\theta\in[0,2\pi)$        | as-hoek, Euler-hoeken,<br />$R_z R_y R_x$             | $(\alpha,\beta)\in S^3 $of $e^{i\theta\sigma_{\vec n}}$ |
| **Topologie**      | cirkel$S^1$                                                     | cirkel$S^1$                | $\mathbb{RP}^3$ (niet 1-samenhangend)                 | 3-sfeer$S^3$ (1-samenhangend)                             |
| **Isomorf met**    | $U(1)$                                                          | $SO(2)$                    | —<br />(dekgroep is$SU(2)$)                          | —<br /> ($SU(2) \to SO(3)$ is 2-op-1)                    |
| **Lie-algebra**    | triviaal (1D, abels)                                              | triviaal                     | $\mathfrak{so}(3)$                                    | $\mathfrak{su}(2) \cong \mathfrak{so}(3)$                 |
| **Irreps**         | $\varphi_n=e^{in\theta}$,<br /> $br$ (1-dim, $\infty$ veel) | idem                         | $l=0,1,2,\dots$ (enkel integer!)                      | $l=0,\tfrac12,1,\tfrac32,\dots$ (alle)                    |
| **Irrep-dimensie** | 1                                                                 | 1                            | $2l+1$                                                | $2l+1$                                                    |

---

## 6. Snelle "hoe leg ik dit uit"-checklist

1. **Wat is het verschil tussen SO(2)/U(1) en SO(3)/SU(2)?**
   → SO(2)/U(1) zijn abels, dus alle irreps zijn 1-dimensionaal en gelabeld door $n\in\mathbb{Z}$. SO(3)/SU(2) zijn niet-abels, dus irreps kunnen hoger-dimensionaal zijn, gelabeld door spin $l$.
2. **Wat is het verschil tussen SO(3) en SU(2) als groep?**
   → Zelfde Lie-algebra (lokaal identiek), maar verschillende topologie globaal: $SU(2)=S^3$ is enkelvoudig samenhangend, $SO(3)$ niet. $SU(2)$ is de 2-op-1 dekgroep van $SO(3)$.
3. **Waarom reken je met sl(2) i.p.v. su(2)?**
   → De ladderoperatoren $L_\pm=L_1\pm iL_2$ zijn complexe combinaties, dus ze leven niet meer in de reële algebra su(2) maar in haar complexificatie sl(2). Er is een 1-op-1 correspondentie tussen hun irreps.
4. **Waarom heeft SO(3) geen half-integer irreps, maar SU(2) wel?**
   → Omdat half-integer irreps van su(2)/sl(2) niet "sluiten" als representatie van de groep $SO(3)$ zelf (enkel als projectieve representatie); een rotatie over $2\pi$ geeft een extra factor $-1$. Dit komt rechtstreeks uit de niet-triviale topologie van $SO(3)$ (§2.3–4.4).
5. **Hoe parametriseer je elke groep?**
   → SO(2)/U(1): één hoek $\theta$. SO(3): drie hoeken (Euler) óf as + hoek $(\vec n,\theta\in[0,\pi])$. SU(2): vier reële getallen op $S^3$, of equivalent as + hoek $(\vec n,\theta\in[0,2\pi))$ via $e^{i\theta\sigma_{\vec n}/1}$ ; met dubbele overdekking t.o.v. SO(3).
