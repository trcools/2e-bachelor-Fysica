# H1 — The origins of quantum theory

## Niet te kennen (zoals in je notities)
- 1.1 *(geen afleidingen)*
- “Finite nuclear mass” in 1.4
- Experimenten van o.a. G.P. Thompson / Jansson / … *(stukje is handgeschreven en niet volledig leesbaar)*

---

## 1.2 Photo-elektrisch effect (Einstein)
- Licht bestaat uit “energiepakketjes” (fotonen) met
  $$
  E = h\nu = \frac{hc}{\lambda}.
  $$

---

## 1.3 Compton-effect (relativistisch nodig)

### Relativistische energie en impuls
Voor een deeltje met rustmassa $m$ en snelheid $v$:
$$
E=\frac{mc^2}{\sqrt{1-\frac{v^2}{c^2}}}, \qquad
p=\frac{mv}{\sqrt{1-\frac{v^2}{c^2}}}.
$$

Hieruit volgt de energie-impulsrelatie:
$$
E^2 = p^2c^2 + m^2c^4
\quad\Rightarrow\quad
E=\sqrt{m^2c^4+p^2c^2}.
$$

Voor een foton geldt $m=0$, dus:
$$
E=pc, \qquad p=\frac{E}{c}=\frac{h}{\lambda}.
$$

### Botsing: behoud van impuls en energie
Situatie: invallend foton $(E_0,p_0)$ botst met elektron in rust; uitgaand foton $(E_1,p_1)$ onder hoek $\theta$; elektron recoilt met impuls $p_2$ onder hoek $\varphi$.

**1) Impulsbehoud (vectorieel):**
$$
\vec p_0=\vec p_1+\vec p_2.
$$
Componenten (zoals in je notities):
$$
p_0 = p_1\cos\theta + p_2\cos\varphi, \qquad
0 = p_1\sin\theta - p_2\sin\varphi.
$$
Hieruit volgt (cosinusregel):
$$
p_2^2=p_0^2+p_1^2-2p_0p_1\cos\theta.
$$

**2) Energiebehoud:**
$$
E_0 + mc^2 = E_1 + \sqrt{m^2c^4+p_2^2c^2}.
$$

Definieer de kinetische energie van het elektron na de botsing:
$$
T_2=\sqrt{m^2c^4+p_2^2c^2}-mc^2.
$$
Dan:
$$
T_2 = E_0 - E_1 = c(p_0-p_1).
$$

Uitwerken geeft (zoals in je afleiding):
$$
p_2^2 = (p_0-p_1)^2 + 2mc(p_0-p_1).
$$

**Combineer** met
$p_2^2=p_0^2+p_1^2-2p_0p_1\cos\theta$ en gebruik $p=h/\lambda$:

$$
\lambda_1-\lambda_0
= \frac{h}{mc}(1-\cos\theta)
= 2\frac{h}{mc}\sin^2\!\left(\frac{\theta}{2}\right).
$$

Met de Compton-golflengte van het elektron
$$
\lambda_C=\frac{h}{mc}
$$
wordt dit:
$$
\Delta\lambda = \lambda_1-\lambda_0 = 2\lambda_C\sin^2\!\left(\frac{\theta}{2}\right).
$$

---

## 1.4 Atomaire spectra en Bohr-model

### Coulombkracht = middelpuntzoekende kracht
Voor 1 elektron rond een kern met ladingsgetal $Z$:
$$
\frac{Ze^2}{4\pi\varepsilon_0\,r^2}=\frac{mv^2}{r}.
$$

### Bohr-postulaat (kwantisatie van impulsmoment)
$$
L = |\vec r\times\vec p| = mvr = n\hbar,\qquad n=1,2,3,\dots
$$

### Straal en snelheid
Hieruit volgt:
$$
v_n=\frac{Ze^2}{4\pi\varepsilon_0\hbar}\frac{1}{n}
= \frac{Z\alpha c}{n},
\qquad
\alpha=\frac{e^2}{4\pi\varepsilon_0\hbar c}.
$$

En:
$$
r_n=\frac{4\pi\varepsilon_0\hbar^2}{Zm e^2}\,n^2
= \frac{a_0}{Z}\,n^2,
\qquad
a_0=\frac{4\pi\varepsilon_0\hbar^2}{m e^2}.
$$

### Energieën
Met $T=\tfrac12 mv^2$ en $V=-\dfrac{Ze^2}{4\pi\varepsilon_0 r}$ krijg je:
$$
E_n=T+V=-\frac{m}{2\hbar^2}\left(\frac{Ze^2}{4\pi\varepsilon_0}\right)^2\frac{1}{n^2}.
$$

(Je notities vermelden ook: als $Z$ niet te groot is, dan is $|E_n| \ll mc^2$ → niet-relativistische aanpak oké.)

### Spectrum / overgangsfrequentie
Voor een overgang $n_i\to n_f$:
$$
h\nu = E_{n_i}-E_{n_f}.
$$
Equivalent in “spectroscopische vorm” (Rydberg-achtig):
$$
\tilde\nu = \frac{1}{\lambda} = RZ^2\left(\frac{1}{n_f^2}-\frac{1}{n_i^2}\right).
$$

### Franck–Hertz (zoals je samenvat)
- Elektronen worden versneld door een spanning $V$.
- Bij bepaalde spanningen treedt **inelastische** botsing op: elektron verliest een vaste energie $E_{\text{exc}}$ om het atoom te exciteren.
- Dit toont: atomen hebben **discrete** excitatie-energieën (kwantisatie).

---

## 1.5 Stern–Gerlach experiment

### Magnetisch moment van een baanstroom (orbitaal)
Modelleer het elektron als een stroomlus:
$$
\mu = IA.
$$
Met $I=\dfrac{ev}{2\pi r}$ en $A=\pi r^2$:
$$
\mu = \frac{evr}{2}.
$$
Omdat $L=mvr$:
$$
\vec\mu_\ell = -\frac{e}{2m}\vec L = -\mu_B\frac{\vec L}{\hbar},
\qquad
\mu_B=\frac{e\hbar}{2m}.
$$

### Energie in een magneetveld en kracht
Potentiële energie:
$$
W=-\vec\mu\cdot\vec B.
$$
Kracht:
$$
\vec F = -\nabla W = \nabla(\vec\mu\cdot\vec B).
$$

In Stern–Gerlach: $B_x=B_y=0$ en $B_z$ is **inhomogeen** (afhankelijk van $z$), dus:
$$
F_x=F_y=0,\qquad
F_z = \mu_z\,\frac{\partial B_z}{\partial z}.
$$

### Klassieke verwachting vs observatie
- Klassiek zou $\mu_z$ continu kunnen variëren → continue spreiding.
- Experiment: **discrete** bundels → kwantisatie van hoekimpulscomponenten.

Quantumresultaten (netjes geformuleerd):
$$
|\vec L|=\sqrt{\ell(\ell+1)}\,\hbar,\qquad
L_z = m_\ell\hbar,\quad m_\ell=-\ell,\dots,+\ell.
$$

### Spin (intrinsiek impulsmoment)
Voor elektronen:
$$
s=\frac12,\qquad
|\vec S|=\sqrt{s(s+1)}\,\hbar,\qquad
S_z=m_s\hbar,\; m_s=\pm\frac12.
$$

Spinmagnetisch moment (met $g_s\approx 2$):
$$
\vec\mu_s = -g_s\,\mu_B\frac{\vec S}{\hbar}.
$$

Totaal magnetisch moment:
$$
\vec\mu = \vec\mu_\ell + \vec\mu_s.
$$









# H2 — De golffunctie en het onzekerheidsbeginsel


## Niet te kennen (zoals in je notities)
- (deze 2 bullets zijn in de scan niet 100% leesbaar)
  - oplossingen in 3D met … + golfpakket … (2.4)
  - twee staten … (experiment/…) … breking, resonantie-absorptie (2.5)

---

## Superpositie en interferentie

Neem een superpositie:
$$
\psi = c_1\psi_1 + c_2\psi_2,
$$
met $c_1,c_2$ complexe constanten en schrijf
$$
\psi_1 = |\psi_1|e^{i\alpha_1}, \qquad \psi_2 = |\psi_2|e^{i\alpha_2}.
$$

De probabiliteitsdichtheid is
$$
|\psi|^2=\psi^\*\psi.
$$

Uitwerken geeft:
$$
|\psi|^2
=|c_1|^2|\psi_1|^2+|c_2|^2|\psi_2|^2
+c_1^\*c_2|\psi_1||\psi_2|e^{i(\alpha_2-\alpha_1)}
+c_2^\*c_1|\psi_1||\psi_2|e^{-i(\alpha_2-\alpha_1)}.
$$

Dus:
$$
|\psi|^2
=|c_1|^2|\psi_1|^2+|c_2|^2|\psi_2|^2
+2\,\mathrm{Re}\!\left\{c_1^\*c_2|\psi_1||\psi_2|e^{i(\alpha_2-\alpha_1)}\right\}.
$$

➡️ Die laatste term is de **interferentieterm** (komt door het complexe faseverschil).

---

## 2.3 Golfpakketten (1D) en groepssnelheid

Een golfpakket kan je schrijven als een superpositie van vlakke golven:
$$
\psi(x,t)=\frac{1}{\sqrt{2\pi\hbar}}
\int dp_x\;
\exp\!\left[\frac{i}{\hbar}\big(p_x x - E(p_x)t\big)\right]\,
\phi(p_x).
$$

Hier is $\phi(p_x)$ een “gewicht”-functie die typisch **gepiekt** is rond $p_x=p_0$.

### Groepssnelheid
$$
v_g=\left.\frac{dE(p_x)}{dp_x}\right|_{p_x=p_0}.
$$

### Taylor rond $p_0$
$$
E(p_x)=E(p_0)+(p_x-p_0)\left.\frac{dE}{dp_x}\right|_{p_0}
+\frac{(p_x-p_0)^2}{2!}\left.\frac{d^2E}{dp_x^2}\right|_{p_0}+\cdots
$$

In eerste orde:
$$
E(p_x)\approx E(p_0)+v_g(p_x-p_0).
$$

---

## Fouriertransformatie (paar)

(Zoals in je notities:)
$$
\psi(x,t)=(2\pi\hbar)^{-1/2}\int_{-\infty}^{+\infty} e^{i p_x x/\hbar}\,\Phi(p_x,t)\,dp_x,
$$
$$
\Phi(p_x,t)=(2\pi\hbar)^{-1/2}\int_{-\infty}^{+\infty} e^{-i p_x x/\hbar}\,\psi(x,t)\,dx.
$$

---

## Hulprelatie: Gauss-integralen

Je gebruikt o.a. de klassieke identiteit
$$
\int_{-\infty}^{+\infty} e^{-x^2}\,dx=\sqrt{\pi}.
$$

In je notities wordt dit afgeleid via:
$$
\left(\int_{-\infty}^{+\infty} e^{-x^2}\,dx\right)
\left(\int_{-\infty}^{+\infty} e^{-y^2}\,dy\right)
=\int\!\!\int_{\mathbb{R}^2} e^{-(x^2+y^2)}\,dx\,dy,
$$
en dan naar poolcoördinaten:
$$
\int_0^{2\pi} d\theta \int_0^\infty r e^{-r^2}\,dr
=2\pi\cdot\frac12=\pi,
$$
dus $I^2=\pi \Rightarrow I=\sqrt{\pi}$.

(Er staat ook een uitwerking met “completing the square” voor integralen van de vorm
$\int e^{-\alpha u^2-\beta u}\,du$.)

---

## Gaussisch golfpakket

Neem een Gauss in impulsmruimte:
$$
\Phi(p_x)=\frac{1}{\sqrt{\Delta p_x\sqrt{\pi}}}\;
\exp\!\left[-\frac{(p_x-p_0)^2}{2(\Delta p_x)^2}\right],
$$
met normalisatieconstante
$$
c=\pi^{-1/4}(\Delta p_x)^{-1/2}
=\frac{1}{\sqrt{\sqrt{\pi}\,\Delta p_x}}.
$$

De tijdsevolutie:
$$
\psi(x,t)=\frac{1}{\sqrt{2\pi\hbar}}
\int_{-\infty}^{+\infty} dp_x\;
e^{i p_x x/\hbar}\,\Phi(p_x)\,e^{-iE(p_x)t/\hbar},
$$
met voor een vrij deeltje:
$$
E(p_x)=\frac{p_x^2}{2m}.
$$

Door de Gauss-integralen krijg je opnieuw een Gauss in $x$ (maar met tijdafhankelijke breedte);
in je notities wordt dit uitgewerkt door de exponent te herschrijven als
$$
-a p_x^2+\beta p_x+\text{const},
$$
en dan de standaard Gauss-integraal toe te passen.

### Onzekerheidsrelatie
Onderaan noteer je:

- voor het **Gaussisch** golfpakket (op $t=0$) is het product minimaal:
  $$
  \Delta x(t=0)\,\Delta p_x = \frac{\hbar}{2},
  $$
- algemeen (Fourier-eigenschap):
  $$
  \Delta x\,\Delta p_x \ge \frac{\hbar}{2}.
  $$




> Examenstof: 3.1 t/m 3.9  
> Niet: **oneindige potentiële energie** (p113) en **Schmidt procedure** (p117–p118).

---

## 3.1 De tijdsafhankelijke Schrödingervergelijking (TDSE)

**Context**
- Geldig voor **niet-relativistische** en **gesloten** systemen.
- Wordt postulaat-achtig ingevoerd (niet “formeel afgeleid” zoals Newton).

**Hamiltoniaan**
- Klassiek (voorbeeld):  
  $$
  E(\mathbf r,\mathbf p,t)=\frac{\mathbf p^2}{2m}+V(\mathbf r,t)
  $$
- Kwantum: vervang door operatoren  
  $$
  \hat{\mathbf r}=\mathbf r,\qquad \hat{\mathbf p}=-i\hbar\nabla,\qquad \hat E=i\hbar\frac{\partial}{\partial t}
  $$
  $$
  \hat H = \frac{\hat{\mathbf p}^2}{2m}+V(\mathbf r,t)
  $$

**TDSE**
- Algemeen:
  $$
  i\hbar\frac{\partial}{\partial t}\Psi(\mathbf r,t)=\hat H\,\Psi(\mathbf r,t)
  $$
- Vrij deeltje (3D):
  $$
  i\hbar\frac{\partial}{\partial t}\Psi(\mathbf r,t)= -\frac{\hbar^2}{2m}\nabla^2\Psi(\mathbf r,t)
  $$

---

## 3.2 Behoud van waarschijnlijkheid

**Waarschijnlijkheidsdichtheid**
$$
P(\mathbf r,t)=|\Psi(\mathbf r,t)|^2=\Psi^*(\mathbf r,t)\Psi(\mathbf r,t)
$$

**Stroomdichtheid (probability current)**
$$
\mathbf j(\mathbf r,t)=\frac{\hbar}{2mi}\Big(\Psi^*\nabla\Psi-(\nabla\Psi^*)\Psi\Big)
=\mathrm{Re}\Big\{\Psi^*\frac{\hbar}{im}\nabla\Psi\Big\}
$$
- Als $\mathrm{Im}\,\Psi=0$ (zuiver reëel), dan volgt vaak $\mathbf j=\mathbf 0$.

**Continuïteitsvergelijking** (voor **reële** potentiaal)
$$
\frac{\partial P}{\partial t}+\nabla\cdot\mathbf j=0
$$

**Normbehoud**
$$
\frac{d}{dt}\int_{\mathbb R^3}|\Psi(\mathbf r,t)|^2\,d^3r = 0
$$
Dus: *eens genormeerd, altijd genormeerd* (mits hermitische $\hat H$).

**Hermiticiteit van $\hat H$**
$$
\int d^3r\,\Psi^*(\hat H\Psi)=\int d^3r\,(\hat H\Psi)^*\Psi
$$
- Reële $V(\mathbf r,t)$  $\Rightarrow$  $\hat H$ hermitisch.
- Imaginair deel in $V$  $\Rightarrow$  typisch niet-hermitisch.

---

## 3.3 Verwachtingswaarden en operatoren

**Basisdefinitie**
$$
\langle \mathbf r\rangle=\int d^3r\,\Psi^*(\mathbf r,t)\,\mathbf r\,\Psi(\mathbf r,t)
$$
Algemeen voor een functie/observeerbare $A$:
$$
\langle A\rangle=\int d^3r\,\Psi^*(\mathbf r,t)\,\hat A\,\Psi(\mathbf r,t)
$$

**“Recept” om $\hat A$ te bouwen**
1. Start met klassieke dynamische variabele $A(\mathbf r,\mathbf p,t)$  
2. Vervang $\mathbf r,\mathbf p$ door operatoren $(\hat{\mathbf r},\hat{\mathbf p})$  
3. “Sandwich”: $\langle A\rangle=\int \Psi^*\hat A\Psi\,d^3r$  
4. Fysisch meetbare $\Rightarrow$ **hermitische** operator (verwachtingswaarde reëel)

**Impulsoperator**
- In coördinatenruimte:
  $$
  \hat{\mathbf p}=-i\hbar\nabla
  $$
en dus
$$
\langle\mathbf p\rangle=\int d^3r\,\Psi^*(-i\hbar\nabla)\Psi
$$

**Energie**
$$
\langle E\rangle=\int d^3r\,\Psi^*\Big(i\hbar\frac{\partial}{\partial t}\Big)\Psi
=\int d^3r\,\Psi^*\,\hat H\,\Psi
$$

**Commutatoren**
$$
[A,B]=AB-BA
$$
Voorbeeld (1D):
$$
[x,p_x]=i\hbar,\qquad [p_x,x]=-i\hbar
$$
Moraal: operatorvolgorde doet ertoe (klassiek vaak niet).

---

## 3.4 Ehrenfest: overgang naar klassieke mechanica

**Klassiek**
$$
\frac{d\mathbf r}{dt}=\frac{\mathbf p}{m},\qquad
\frac{d\mathbf p}{dt}=-\nabla V
$$

**Kwantum (Ehrenfest)**
$$
\frac{d\langle\mathbf r\rangle}{dt}=\frac{\langle\mathbf p\rangle}{m},\qquad
\frac{d\langle\mathbf p\rangle}{dt}=-\langle\nabla V\rangle
$$

**Algemene formule**
$$
\frac{d}{dt}\langle A\rangle=\frac{1}{i\hbar}\langle[\hat A,\hat H]\rangle
+\Big\langle\frac{\partial \hat A}{\partial t}\Big\rangle
$$

Interpretatie (niet-mystiek maar wel mooi):
- Klassiek: één traject.
- Kwantum: verwachtingswaarden volgen “klassiek-achtige” bewegingswetten, maar de toestand blijft statistisch/gespreid.

---

## 3.5 Tijdonafhankelijke Schrödingervergelijking (TISE) — stationaire toestanden

Neem een **tijdonafhankelijke potentiaal** $V(\mathbf r)$.

**Scheiding der variabelen**
$$
\Psi(\mathbf r,t)=\psi_E(\mathbf r)\,T(t)
$$
Leidt tot
$$
T(t)=e^{-iEt/\hbar}
$$
en de **TISE**
$$
\hat H\psi_E(\mathbf r)=E\psi_E(\mathbf r),
\qquad
\hat H=-\frac{\hbar^2}{2m}\nabla^2+V(\mathbf r)
$$

**Stationaire toestand**
$$
\Psi_E(\mathbf r,t)=\psi_E(\mathbf r)\,e^{-iEt/\hbar}
$$
Dan is
$$
|\Psi_E(\mathbf r,t)|^2 = |\psi_E(\mathbf r)|^2
$$
(tijdonafhankelijke dichtheid).

---

## 3.6 Energiekwantisatie

**Waarom kwantisatie?**
- Fysische randvoorwaarden (normaliseerbaarheid, continuïteit, grensvoorwaarden) laten vaak enkel **discrete** $E$ toe.

**Typen spectra**
- **Gebonden toestanden**: discrete energieën $E_1<E_2<\dots$
- **Ongebonden / verstrooiing**: typisch continu spectrum (energie “loopt door”).

---

## 3.7 Eigenschappen van energie-eigenfuncties

**Hermiticiteit ⇒ reële eigenwaarden**
- Voor hermitische $\hat H$ zijn energieën $E$ reëel.

**Orthogonaliteit / orthonormaliteit (schets)**
- Voor discrete niveaus:
  $$
  \int d^3r\,\psi_n^*(\mathbf r)\psi_m(\mathbf r)=\delta_{nm}
  $$
- (Voor continue spectra verschijnt een Dirac-delta i.p.v. Kronecker-delta.)

**Oscillatietheorema (1D, discrete niveaus)**
- Voor $E_1<E_2<E_3<\dots$ en bijhorende $\psi_1(x),\psi_2(x),\dots$:
  - $\psi_n(x)$ heeft **$(n-1)$ nuldoorgangen** (“zeros”) voor eindige $x$.

> Niet in de examenstof: Schmidt-procedure (Gram–Schmidt orthonormalisatie van ontaarde eigenfuncties).

---

## 3.8 Algemene oplossing van de TDSE voor tijdonafhankelijke potentiaal

**Superpositie in energie-eigenfuncties**
$$
\Psi(\mathbf r,t)=\sum_E C_E(t)\,\psi_E(\mathbf r)
$$
Voor tijdonafhankelijke $\hat H$:
$$
C_E(t)=C_E(t_0)\exp\Big[-\frac{i}{\hbar}E(t-t_0)\Big]
$$
Dus
$$
\Psi(\mathbf r,t)=\sum_E c_E\,e^{-iEt/\hbar}\,\psi_E(\mathbf r)
$$

**Coefficients uit beginvoorwaarde**
$$
C_E(t_0)=\int d^3r\,\psi_E^*(\mathbf r)\Psi(\mathbf r,t_0)
$$

**Propagator**
$$
\Psi(\mathbf r,t)=\int d^3r'\,K(\mathbf r,t;\mathbf r',t_0)\,\Psi(\mathbf r',t_0)
$$
met
$$
K(\mathbf r,t;\mathbf r',t_0)=\sum_E \psi_E^*(\mathbf r')\psi_E(\mathbf r)\,
\exp\Big[-\frac{i}{\hbar}E(t-t_0)\Big]
$$

---

## 3.9 Schrödingervergelijking in de impulsruimte

**Fouriertransformatie tussen representaties**
$$
\Psi(\mathbf r,t)=\frac{1}{(2\pi\hbar)^{3/2}}\int d^3p\,\Phi(\mathbf p,t)\,
e^{+i\mathbf p\cdot\mathbf r/\hbar}
$$
$$
\Phi(\mathbf p,t)=\frac{1}{(2\pi\hbar)^{3/2}}\int d^3r\,\Psi(\mathbf r,t)\,
e^{-i\mathbf p\cdot\mathbf r/\hbar}
$$

**TDSE in impulsruimte (integro-differentiaal)**
$$
i\hbar\frac{\partial}{\partial t}\Phi(\mathbf p,t)=\frac{\mathbf p^2}{2m}\Phi(\mathbf p,t)
+\frac{1}{(2\pi\hbar)^{3/2}}\int d^3p'\,V(\mathbf p-\mathbf p',t)\,\Phi(\mathbf p',t)
$$
- Komt uit convolutietheorema; vaak moeilijker dan in $\mathbf r$-ruimte.

**Verwachtingswaarden in impulsruimte**
$$
\langle\mathbf p\rangle=\int d^3p\,\Phi^*(\mathbf p,t)\,\mathbf p\,\Phi(\mathbf p,t)
$$
$$
\langle\mathbf r\rangle=\int d^3p\,\Phi^*(\mathbf p,t)\,(i\hbar\nabla_{\mathbf p})\,\Phi(\mathbf p,t)
$$
Algemeen:
$$
\langle A\rangle=\int d^3p\,\Phi^*(\mathbf p,t)\,A(i\hbar\nabla_{\mathbf p},\mathbf p,t)\,\Phi(\mathbf p,t)
$$

# H3 — Schrödingervergelijking 

---

Tijdsafhankelijke Schrödingervergelijking:
$$
i\hbar\frac{\partial}{\partial t}\Psi(\mathbf r,t)=\hat H\,\Psi(\mathbf r,t)
$$
Hamiltoniaan:
$$
\hat H=-\frac{\hbar^2}{2m}\nabla^2+V(\mathbf r,t)
$$

Opmerking:
- als $\Psi$ reëel is $\Rightarrow \mathbf j=0$ `[zoals genoteerd]`

Vlakke golf (free particle):
$$
\Psi(\mathbf r,t)=A\,e^{i(\mathbf k\cdot\mathbf r-\omega t)}
$$

## — Operatoren + continuïteitsvergelijking

Operatoren:
$$
\hat E=i\hbar\frac{\partial}{\partial t},\qquad
\hat{\mathbf L}=\mathbf r\times \hat{\mathbf p}
$$
Kinetische energie:
$$
\hat T=\frac{\hat p^2}{2m}=-\frac{\hbar^2}{2m}\nabla^2
$$

Continuïteitsvergelijking (uit TDSE):
$$
P(\mathbf r,t)=|\Psi(\mathbf r,t)|^2
$$
$$
\mathbf j(\mathbf r,t)=\frac{\hbar}{2mi}\left(\Psi^*\nabla\Psi-(\nabla\Psi^*)\Psi\right)
$$
$$
\boxed{\;\frac{\partial P}{\partial t}+\nabla\cdot\mathbf j=0\;}
$$

---

##  Expectation values and operators

Configuratieruimte:
$$
\langle \mathbf r\rangle=\int \Psi^*(\mathbf r,t)\,\mathbf r\,\Psi(\mathbf r,t)\,d^3r
$$
Voor $f(\mathbf r,t)$:
$$
\langle f(\mathbf r,t)\rangle=\int \Psi^*\,f(\mathbf r,t)\,\Psi\,d^3r
$$
Voor potentiaal:
$$
\langle V(\mathbf r)\rangle=\int \Psi^*(\mathbf r,t)\,V(\mathbf r)\,\Psi(\mathbf r,t)\,d^3r
$$

Voor dynamische variabele $A(\mathbf p,\mathbf r,t)$:
$$
\langle A\rangle=\int \Psi^*(\mathbf r,t)\,\hat A\,\Psi(\mathbf r,t)\,d^3r
$$

Belangrijke operatoren:
- configuratieruimte:
$$
\hat{\mathbf p}=-i\hbar\nabla
$$
- impulsmruimte:
$$
\hat{\mathbf p}=\mathbf p,\qquad \hat{\mathbf r}=i\hbar\nabla_{\mathbf p}
$$

---

**Stap 1: vanuit golffunctie in p-ruimte**
$$
\psi(\mathbf r) = (2\pi\hbar)^{-3/2}\int \Phi(\mathbf p)\,e^{\frac{i}{\hbar}\mathbf p\cdot\mathbf r}\,d^3p
$$
$$
\Phi(\mathbf p) = (2\pi\hbar)^{-3/2}\int \psi(\mathbf r)\,e^{-\frac{i}{\hbar}\mathbf p\cdot\mathbf r}\,d^3r
$$

**Stap 2: van p-ruimte naar r-ruimte**
$$
\langle \mathbf p\rangle=\int \Phi^*(\mathbf p)\,\mathbf p\,\Phi(\mathbf p)\,d^3p
$$
Substitutie van de Fourierrelaties geeft een drievoudige integraal in $\mathbf r,\mathbf r',\mathbf p$.

**Stap 3: gebruik**
$$
\mathbf p\,e^{\frac{i}{\hbar}\mathbf p\cdot(\mathbf r-\mathbf r')}
= -\,i\hbar\,\nabla_{\mathbf r}\,e^{\frac{i}{\hbar}\mathbf p\cdot(\mathbf r-\mathbf r')}
$$

**Stap 4: p-integraal**
$$
\int e^{\frac{i}{\hbar}\mathbf p\cdot(\mathbf r-\mathbf r')}\,d^3p=(2\pi\hbar)^3\delta(\mathbf r-\mathbf r')
$$

**Stap 5: resultaat**
$$
\boxed{\;\langle \mathbf p\rangle=\int \psi^*(\mathbf r)\,(-i\hbar\nabla)\,\psi(\mathbf r)\,d^3r\;}
$$

---

## Verwachtingswaarde energie + operatoren

Verwachtingswaarde:
$$
\langle \hat H\rangle=\left\langle \frac{\hat p^2}{2m}+V\right\rangle
=\langle \hat T\rangle+\langle V\rangle
$$

Hermitische operatoren (kernpunten):
- eigenwaarden reëel
- eigenfuncties orthogonaal (bij verschillende eigenwaarden)
- $\hat A^\dagger=\hat A$  `[onleesbaar detail]`

---



## Ehrenfest: bewijs voor $\frac{d}{dt}\langle x\rangle$

Start:
$$
\frac{d}{dt}\langle x\rangle=\frac{d}{dt}\int \Psi^*\,x\,\Psi\,d^3r
$$
Na uitwerken:
$$
\boxed{\;\frac{d}{dt}\langle x\rangle=\frac{1}{m}\langle p_x\rangle\;}
$$

---


## — Ehrenfest: bewijs voor $\frac{d}{dt}\langle p_x\rangle$

Start:
$$
\frac{d}{dt}\langle p_x\rangle=\frac{d}{dt}\int \Psi^*(-i\hbar\partial_x)\Psi\,d^3r
$$
Gebruik TDSE + integraalmanipulaties (partiële integratie, randtermen $=0$):
$$
\boxed{\;\frac{d}{dt}\langle p_x\rangle=-\left\langle\frac{\partial V}{\partial x}\right\rangle\;}
$$

---


## 3.5 The time-independent Schrödinger equation (stationary state)

Ansatz:
$$
\Psi(\mathbf r,t)=\psi(\mathbf r)\,f(t)
$$
In TDSE:
$$
i\hbar\frac{\partial}{\partial t}(\psi f)=\hat H(\psi f)
\Rightarrow i\hbar \psi \frac{df}{dt}=f\,\hat H\psi
$$
Deel door $\psi f$:
$$
i\hbar\frac{1}{f}\frac{df}{dt}=\frac{1}{\psi}\hat H\psi = E
$$
Dus:
$$
f(t)=e^{-\frac{i}{\hbar}Et},\qquad \hat H\psi=E\psi
$$
Superpositie geeft de algemene oplossing:
$$
\Psi(\mathbf r,t)=\sum_E c_E\,\psi_E(\mathbf r)\,e^{-\frac{i}{\hbar}Et}
$$

---

## Stationaire toestanden (eigenschappen)

Stationaire oplossing:
$$
\Psi(\mathbf r,t)=\psi(\mathbf r)\,e^{-\frac{i}{\hbar}Et}
$$
met $\hat H\psi=E\psi$ en (voor hermitische $\hat H$) $E\in\mathbb R$.

Verwachtingswaarde energie in eigenstaat:
$$
\langle E\rangle=\int \Psi^*(i\hbar\partial_t)\Psi\,d^3r
=\int \Psi^*\hat H\Psi\,d^3r
=E\int |\psi|^2\,d^3r=E
$$

Eigenschappen stationaire toestand:
- $P(\mathbf r,t)=|\Psi|^2=|\psi(\mathbf r)|^2$ (tijdonafhankelijk)
- $\mathbf j(\mathbf r,t)$ tijdonafhankelijk `[onleesbaar detail]`

Als $\hat A$ niet expliciet tijdsafhankelijk is:
$$
\langle A\rangle=\int \psi^*(\mathbf r)\,\hat A\,\psi(\mathbf r)\,d^3r
$$

---

## 3.6 Energie-kwantisatie

1D-TISE:
$$
-\frac{\hbar^2}{2m}\frac{d^2\psi(x)}{dx^2}+V(x)\psi(x)=E\psi(x)
$$

Rand/regulariteitsvoorwaarden (typisch):
- $\psi(x)$ eindig / normaliseerbaar
- $\psi$ continu; $d\psi/dx$ continu als $V$ eindig

Cases (zoals genoteerd):
- Case 1: $E < V_{\min}$
- Case 2: $V_{\min} < E < V_-$  `[onleesbaar]`
- Case 3: $V_- < E < V_+$
- Case 4: $E > V_+$

Interpretatie (klassiek):
- $V(x)\le E$: klassieke beweging mogelijk → oplossingen oscillatoir
- $V(x)>E$: klassieke beweging onmogelijk → oplossingen exponentieel

---


## Coëfficiënten, propagator, algemene oplossing

Projectiecoëfficiënten:
$$
C_E(t_0)=\int \psi_E^*(\mathbf r)\,\Psi(\mathbf r,t_0)\,d^3r
$$

Propagator:
$$
\Psi(\mathbf r,t)=\int K(\mathbf r,t;\mathbf r',t_0)\,\Psi(\mathbf r',t_0)\,d^3r'
$$
$$
K(\mathbf r,t;\mathbf r',t_0)=\sum_E \psi_E(\mathbf r)\psi_E^*(\mathbf r')\,
\exp\!\left[-\frac{i}{\hbar}E(t-t_0)\right]
$$

Algemene oplossing (tijdonafhankelijke $\hat H$):
$$
\Psi(\mathbf r,t)=\sum_E C_E(t_0)\exp\!\left[-\frac{i}{\hbar}E(t-t_0)\right]\psi_E(\mathbf r)
$$

Normalisatie:
$$
\int |\Psi(\mathbf r,t_0)|^2\,d^3r=1
\quad\Rightarrow\quad
\sum_E |C_E|^2=1
$$

Probabiliteitsdichtheid:
$$
P(\mathbf r,t)=\Psi^*(\mathbf r,t)\Psi(\mathbf r,t)
=\sum_E |C_E|^2|\psi_E(\mathbf r)|^2
+\sum_{E\neq E'} C_E^*C_{E'}\,e^{\frac{i}{\hbar}(E-E')t}\,\psi_E^*\psi_{E'}
$$
(laatste term = interferentie)

---



## 3.3 The Schrödinger equation in momentum space (bewijs ufora)

Definieer de Fouriertransformatie van het potentiaal:
$$
\tilde V(\mathbf p-\mathbf p',t)=(2\pi\hbar)^{-3/2}\int
e^{-\frac{i}{\hbar}(\mathbf p-\mathbf p')\cdot\mathbf r}\,V(\mathbf r,t)\,d^3r
$$

**TDSE in impulsmruimte**
$$
i\hbar\frac{\partial}{\partial t}\Phi(\mathbf p,t)
=\frac{p^2}{2m}\Phi(\mathbf p,t)
+(2\pi\hbar)^{-3/2}\int \tilde V(\mathbf p-\mathbf p',t)\,\Phi(\mathbf p',t)\,d^3p'
$$

**Verwachtingswaarden in impulsmruimte**
$$
\langle \mathbf p\rangle=\int \Phi^*(\mathbf p,t)\,\mathbf p\,\Phi(\mathbf p,t)\,d^3p
$$
$$
\langle \mathbf r\rangle=\int \Phi^*(\mathbf p,t)\,(i\hbar\nabla_{\mathbf p})\,\Phi(\mathbf p,t)\,d^3p
$$
Algemeen:
$$
\langle A\rangle=\int \Phi^*(\mathbf p,t)\,A(i\hbar\nabla_{\mathbf p},\mathbf p,t)\,\Phi(\mathbf p,t)\,d^3p
$$

---


# 4 One-dimensional examples


## 4.1 General formulae

TDSE:
\[
i\hbar \frac{\partial}{\partial t}\psi(x,t)=
\left[-\frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2}+V(x)\right]\psi(x,t)
\]

TISE:
\[
-\frac{\hbar^2}{2m}\frac{d^2}{dx^2}\psi(x)+V(x)\psi(x)=E\psi(x)
\]

\[
\vec{j}(x,t)=\frac{\hbar}{2mi}\left[\psi^*(x,t)\frac{\partial\psi(x,t)}{\partial x}-\psi(x,t)\frac{\partial\psi^*(x,t)}{\partial x}\right]
\]


## 4.2 The free particle

\[
E=\frac{\hbar^2k^2}{2m}
\]
\[
(V=0)\Rightarrow -\frac{\hbar^2}{2m}\frac{d^2\psi(x)}{dx^2}=E\psi(x)
\Rightarrow \frac{d^2\psi(x)}{dx^2}+\frac{2mE}{\hbar^2}\psi(x)=0
\Rightarrow \frac{d^2\psi(x)}{dx^2}+k^2\psi(x)=0
\]

dimensions: \[onleesbaar\] (iets met \(\left(\frac{2mE}{\hbar^2}\right)^{1/2}\))

general solution:
\[
\psi(x)=Ae^{ikx}+Be^{-ikx}
\]

interpretatie:

1) \(Ae^{ikx}\): vrij deeltje met impuls \(+\hbar k\) in de positieve \(x\)-richting

2) \(Be^{-ikx}\): vrij deeltje met impuls \(-\hbar k\) in de negatieve \(x\)-richting

\[
\psi(x,t)=Ae^{i(kx-\omega t)}+Be^{-i(kx+\omega t)}
\]

stationaire toestand: meting van de energie van het vrij deeltje zal steeds de waarde \(E=\hbar\omega\) opleveren.


Geval 1: \(B=0\) (vrij deeltje in de positieve \(x\)-richting)

- \(P(x,t)=P(x)=|A|^2\)

\[
\Delta x\to\infty \quad \text{want}\quad \Delta p_x=0
\]

- waarschijnlijkheidsdichtheid:
\[
j(x,t)=j(x)=\frac{\hbar k}{m}|A|^2=v|A|^2=vP(x)
\]

\[
j=v\cdot P \qquad \text{(flux = snelheid}\cdot\text{dichtheid)}
\]


Geval 2: \(A=0\) analoge resultaten

\[
P(x,t)=P(x)=|B|^2
\]

\[
j(x,t)=j(x)=-v|B|^2=-vP(x)
\]


Geval 3: \(A=B\) (twee vlakke golven die met gelijke amplitude in tegengestelde richting bewegen)

- \[
\psi(x,t)=2(A\cos(kx))\,e^{-i\omega t}
\]

\[
x_m\ \text{nodes}=\pm\frac{\left(\frac{\pi}{2}+m\pi\right)}{k}\ ,\quad m=0,1,2,\ldots
\]

- \[
P(x,t)=P(x)=|2A|^2\cos^2(kx)
\]

- \[
j(x,t)=j(x)=0
\]
(geen netto stroom, want \(\mathrm{im}(\psi^*(x,t)\frac{\partial\psi(x,t)}{\partial x})=0\))


Geval 4: \(A=-B\)

\[
P(x,t)=P(x)=|2A|^2\sin^2(kx)
\]

\[
j(x,t)=j_0=0
\]


Geval 5: \(A\) en \(B\) arbitrair (superpositie van 2 vlakke golven die in tegengestelde richting bewegen resp. met amplitude A en B)

\[
\psi(x,t)=Ae^{i(kx-\omega t)}+Be^{-i(kx+\omega t)}
\]

\[
P(x,t)=P(x)=|A|^2+|B|^2+\left(AB^*e^{2ikx}+A^*Be^{-2ikx}\right)
\]

met interferentie effecten tussen de 2 golven

\[
j(x,t)=j(x)=\frac{\hbar}{m}\ \mathrm{Im}\left(\psi^*(x,t)\frac{\partial\psi(x,t)}{\partial x}\right)
=v\left(|A|^2-|B|^2\right)
\]

netto flux van 2 golven: 1 in de positieve \(x\)-richting (amplitude A) en 1 in de negatieve \(x\)-richting (amplitude B)


* De golffunctie voor een vrij deeltje is niet kwadraat integreerbaar!!

In de praktijk wordt de beweging van het vrij deeltje vaak beperkt tot 1D met behulp van box normalisatie.

vrij deeltje in eindige doos:
met periodieke randvoorwaarden opgelegd: \(\psi_n(x=0)=\psi_n(x=L)\)

mogelijke impuls eigenwaarden en energie eigenwaarden:
\[
k_n=\frac{2\pi n}{L}\ ,\qquad
E_n=\frac{\hbar^2k_n^2}{2m}=\frac{2\pi^2\hbar^2}{mL^2}n^2
\quad (n=0,\pm1,\pm2,\ldots)
\]

de golffunctie kan genormeerd worden in de basis-doos \((0\le x\le L)\):
\[
\psi_n(x)=\frac{1}{\sqrt{L}}e^{ik_nx}
\qquad
\text{met}\quad
\int_0^L \psi_n^*(x)\psi_n(x)\,dx=1
\]


## 4.4 The potential barrier

```
(I)            (II)                 (III)
─────| 0        | a ───────────────────────> x
    |----------|
    V0
      V0
```
*(Schematische barrière: regio I links, regio II hoogt​e \(V_0\) tussen 0 en \(a\), regio III rechts.)*


\[
V(x)=
\begin{cases}
0, & x<0 \\
V_0, & 0<x<a \\
0, & x>a
\end{cases}
\]

**Region I:** \(-\infty<x\le 0\)

\[
\psi(x)=Ae^{i\ell x}+Be^{-i\ell x}
\]
\[
E=\frac{\hbar^2\ell^2}{2m}
\]
\[
j_0=v(|A|^2-|B|^2)=j(\text{incoming})-j(\text{reflected})
\]

**Region II:** \(0\le x\le a\)

2 gevallen:
- **geval 1:** \(E<V_0\)
\[
\psi(x)=Fe^{Kx}+Ge^{-Kx}
\qquad \text{met}\quad
K=\left[\frac{2m(V_0-E)}{\hbar^2}\right]^{1/2}
\]

- **geval 2:** \(E>V_0\)
\[
\psi(x)=Fe^{i\ell' x}+Ge^{-i\ell' x}
\qquad \text{met}\quad
\ell'=\left[\frac{2m(E-V_0)}{\hbar^2}\right]^{1/2}
\]

**Region III:** \(a<x<+\infty\)

\[
\psi(x)=Ce^{i\ell x}
\quad
(\text{want deeltje alleen in de positieve }x\text{-richting beweegt})
\]
\[
E=\frac{\hbar^2\ell^2}{2m}
\]
\[
j_0=v|C|^2=j(\text{doorgelaten})
\]

* Reflectie-coëfficiënt \(R\):
\[
R=\frac{j(\text{reflected})}{j(\text{incoming})}
\]

* Transmissie-coëfficiënt \(T\):
\[
T=\frac{j(\text{doorgelaten})}{j(\text{incoming})}
\]

R en T worden berekend via de waarschijnlijkheidsstromen \(j(x,t)\)
(niet via de positie waarschijnlijkheid \(P(x,t)\)).


\[
R=\frac{|B|^2}{|A|^2}
\qquad ; \qquad
T=\frac{|C|^2}{|A|^2}
\]
\[
R+T=\frac{|B|^2+|C|^2}{|A|^2}=1
\]

Terug naar geval 1: \(E<V_0\)

De golffunctie op \(x=0\) en \(x=a\) is continu, dus:

met \(Ae^{i\ell x}+Be^{-i\ell x}=Fe^{Kx}+Ge^{-Kx}\) op \(x=0\) en \(x=a\) (lijnen), vinden we:
\[
\begin{cases}
A+B=F+G \\
Fe^{Ka}+Ge^{-Ka}=Ce^{i\ell a}
\end{cases}
\]

De \[onleesbaar\] afgeleiden op \(x=0,\ x=a\) moeten ook gelijk zijn om in dezelfde fase te blijven:
\[
\begin{cases}
i\ell(A-B)=K(F-G) \\
K(Fe^{Ka}-Ge^{-Ka})=i\ell Ce^{i\ell a}
\end{cases}
\]

We elimineren \(F\) en \(G\), en bekomen \(\frac{B}{A}\) en \(\frac{C}{A}\):

\[
\frac{B}{A}=
\frac{(\ell^2+K^2)(e^{2Ka}-1)}
{e^{2Ka}(\ell+iK)^2-(\ell-iK)^2}
\]

\[
\frac{C}{A}=
\frac{4i\ell K\,e^{-i\ell a}\,e^{Ka}}
{e^{2Ka}(\ell+iK)^2-(\ell-iK)^2}
\]

\[
T=\frac{|C|^2}{|A|^2}=
\left[1+\frac{(\ell^2+K^2)^2\sinh^2(Ka)}{4\ell^2K^2}\right]^{-1}=
\left[1+\frac{V_0^2\sinh^2(Ka)}{4E(V_0-E)}\right]^{-1}
\]

\[
R=\frac{|B|^2}{|A|^2}=
\left[1+\frac{4\ell^2K^2}{(\ell^2+K^2)^2\sinh^2(Ka)}\right]^{-1}=
\left[1+\frac{4E(V_0-E)}{V_0^2\sinh^2(Ka)}\right]^{-1}
\]

Nu kunnen we nog eens verifiëren dat \(R+T=1\) wat intuitief ook logisch is.


\[
\frac{mV_0a^2}{\hbar^2}
:\ \text{opaciteit van de potentiaal barriere}
\quad
(\text{maat van kwantum mechanische ondoordringbaarheid})
\]

\[
T \approx \frac{16E(V_0-E)}{V_0^2}\,e^{-2Ka}
\quad
(\text{toepassing in scanning tunneling microscope})
\]

Terug naar Case 2: \(E>V_0\)

\[
\psi(x)=Fe^{i\ell' x}+Ge^{-i\ell' x}
\qquad \text{met}\quad
\ell'=\left[\frac{2m(E-V_0)}{\hbar^2}\right]^{1/2}
\]

Analoge oplossing probleem vinden we:

\[
R=\frac{|B|^2}{|A|^2}=
\left[1+\frac{4\ell^2\ell'^2}{(\ell^2-\ell'^2)^2\sin^2(\ell'a)}\right]^{-1}=
\left[1+\frac{4E(E-V_0)}{V_0^2\sin^2(\ell'a)}\right]^{-1}
\]

\[
T=\frac{|C|^2}{|A|^2}=
\left[1+\frac{V_0^2\sin^2(\ell'a)}{4E(E-V_0)}\right]^{-1}
\]

gemakkelijk zien we dat \(R+T=1\)

klassiek: \(R=0\) en \(T=1\)

kwantum: \(R\ne 0\) en \(T<1\)

\[
(\ell a)^2=\left(\frac{2mV_0a^2}{\hbar^2}\right)\left(\frac{E}{V_0}\right)
\]
\[
(Ka)^2=\left(\frac{2mV_0a^2}{\hbar^2}\right)\left(1-\frac{E}{V_0}\right)
\]

T hangt enkel af van 2 dimensieloze grootheden:
de opaciteit \(\frac{mV_0a^2}{\hbar^2}\) en de ratio \(\frac{E}{V_0}\).


## 4.5 The infinite square well

\[
V(x)=
\begin{cases}
0, & -a<x<a \\
\infty, & |x|>a
\end{cases}
\]

\[
\psi(x=\pm a)=0
\]

Schrödinger for \(|x|<a\):
\[
-\frac{\hbar^2}{2m}\frac{d^2\psi(x)}{dx^2}=E\psi(x)
\Rightarrow \frac{d^2\psi(x)}{dx^2}+k^2\psi(x)=0
\qquad
\left(k=\left(\frac{2mE}{\hbar^2}\right)^{1/2}\right)
\]

general solution:
\[
\psi(x)=A\cos(kx)+B\sin(kx)
\]

if we apply the boundary conditions:
\[
A\cos(ka)=0 \quad \text{en} \quad B\sin(ka)=0
\]

1) \(B=0\Rightarrow \cos(ka)=0\)

The only allowed values of \(k\) are:
\[
k_m=\frac{m\pi}{2a}=\frac{m\pi}{L}
\qquad (m=1,3,5,\ldots)
\]

The corresponding eigenfunction
\[
\psi_m(x)=A_m\cos(k_m x)
\]
can be normalized so that:
\[
\int_{-a}^{a}\psi_m^*(x)\psi_m(x)\,dx=1
\]

The normalized eigenfunctions are:
\[
\psi_m(x)=\frac{1}{\sqrt{a}}\cos\!\left(\frac{m\pi}{2a}x\right)
\qquad (m=1,3,5,\ldots)
\]

2) \(A=0\Rightarrow \sin(ka)=0\)

The normalized eigenfunctions are:
\[
\psi_m(x)=\frac{1}{\sqrt{a}}\sin\!\left(\frac{m\pi}{2a}x\right)
\qquad (m=2,4,\ldots)
\]

\[
k_m=\frac{m\pi}{L}
\qquad
\lambda_m=\frac{2\pi}{k_m}=\frac{2L}{m}
\]

The energy is quantized, the energy eigenvalues being:
\[
E_m=\frac{\hbar^2 k_m^2}{2m}
=\frac{\hbar^2\pi^2 m^2}{8ma^2}
=\frac{\hbar^2\pi^2 m^2}{2mL^2}
\qquad (m=1,2,3,\ldots)
\]


## (p.8) vervolg: spectrum + orthogonaliteit + zero-point energy

The energy spectrum obtained consists of an infinite number of discrete energy levels corresponding to bound states.

eigenfunctions \(\psi_n(x)\) and \(\psi_m(x)\) corresponding to different eigenvalues \(E_n\) and \(E_m\) are orthogonal:
\[
\int_{-a}^{a}\psi_n^*(x)\psi_m(x)\,dx=0
\qquad n\ne m
\]

Note that the lowest energy (zero-point energy) is:
\[
E_1=\frac{\hbar^2\pi^2}{8ma^2}
\]
so there is no state of zero-energy!!

uncertainty principle:
\[
\Delta x\approx a \ \Rightarrow\ \Delta p_x \approx \frac{\hbar}{a}
\]
indicates an minimal kinetic energy:
\[
\frac{(\Delta p_x)^2}{2m}=\frac{\hbar^2}{2ma^2}=\frac{4E_1}{\pi^2}
\]

* Parity

The eigenfunctions of the cosine function are such that:
\[
\psi_n(-x)=\psi_n(x)
\]
and are therefore even function of \(x\),

while the eigenfunctions with a sine function in it are such that:
\[
\psi_n(-x)=-\psi_n(x)
\]
and hence are odd.

Reflection through the origin \(x\to -x\) is called the parity operation.

If the potential is symmetric the Hamiltonian
\[
H=-\frac{\hbar^2}{2m}\frac{d^2}{dx^2}+V(x)
\]
does not change when \(x\) is replaced by \(-x\) (invariant under parity operation).


## (p.9) Case 1–2 (degeneracy)  [gedeeltelijk onleesbaar]

Case 1: The eigenvalue \(E\) is non-degenerate \([onleesbaar]\)

The eigenfunctions \(\psi_1(x)\) and \(\psi_2(x)\) corresponding to the same energy eigenvalue \(E\) can thus only differ by a multiplicative constant:
\[
\psi_1(x)=c\,\psi_2(x)
\]

Combining the above with the parity equation yields: \([onleesbaar]\)

and hence we may choose the eigenfunction \(\psi(x)\) have a definite parity being either even or odd.

Case 2: The eigenvalue \(E\) is degenerate \([onleesbaar]\)

In this case more than one linearly independent eigenfunction corresponds to the same energy eigenvalue. \([onleesbaar]\)

However, it is always possible to construct linear combinations of these eigenfunctions which do have definite parity. \([onleesbaar]\)

The Schrödinger equation becomes: \([onleesbaar]\)

\[
\left[-\frac{\hbar^2}{2m}\left(\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}+\frac{\partial^2}{\partial z^2}\right)+V(x,y,z)\right]\psi(x,y,z)=E\psi(x,y,z)
\]

Using parity operation \([onleesbaar]\)

We can have \(P\psi=\pm\psi\) \([onleesbaar]\)

The commutation equation becomes: \([onleesbaar]\)


* wave function representation

general solution TDSE:
\[
\psi(x,t)=\sum_n c_n\,\psi_n(x)\,e^{-iE_nt/\hbar}
\]

\[
c_n=\int \psi_n^*(x)\,\psi(x,t=0)\,dx
\]

\[
\psi(x,t=T)=\psi(x,t=0)
\]

Since \(E_n=nE_1\) we have
\[
\psi(x,t+T)=\sum_n c_n\psi_n(x)\exp\!\left(-i\frac{nE_1(t+T)}{\hbar}\right)
=\sum_n c_n\psi_n(x)e^{-inE_1t/\hbar}\exp\!\left(-i\frac{nE_1T}{\hbar}\right)
\]

Choose
\[
T=\frac{2\pi\hbar}{E_1}
\Rightarrow
\exp\!\left(-i\frac{nE_1T}{\hbar}\right)=\exp(-i2\pi n)=1
\]

and thus
\[
\psi(x,t+T)=\psi(x,t)
\]

----------------------------------------------------------------

4.7 The harmonic oscillator

Energy:
\[
E=\frac{1}{2}m\dot{x}^{\,2}+\frac{1}{2}kx^2
=\frac{1}{2}m\dot{x}^{\,2}+V(x)
\]

Potential:
\[
V(x)=\frac{1}{2}kx^2=\frac{1}{2}m\omega^2x^2
\qquad\text{met}\qquad
\omega=\sqrt{\frac{k}{m}}
\]

* LHO in de klassieke mechanica

\[
\frac{d^2x}{dt^2}=-\omega^2x
\Rightarrow x(t)=x_0\cos(\omega t)
\qquad
(\text{waar }x_0\text{ bepaald wordt door de startcondities})
\]

\[
v_{\max}=\omega x_0
\qquad
T=\frac{2\pi}{\omega}\ \ \text{(interval }(-x_0,x_0)\text{)}
\]

de “turning points” (waar \(v=0\)):
\[
\frac{1}{2}kx_0^2=E
\Rightarrow
x_0=\sqrt{\frac{2E}{k}}=\sqrt{\frac{2E}{m\omega^2}}
\]

* LHO in de quantum mechanica

TISE:
\[
-\frac{\hbar^2}{2m}\frac{d^2\psi(x)}{dx^2}+\frac{1}{2}m\omega^2x^2\,\psi(x)=E\psi(x)
\]

dimensionless eigenvalues:
\[
\lambda=\frac{2E}{\hbar\omega}
\qquad\text{where}\qquad
\omega=\left(\frac{k}{m}\right)^{1/2}
\]

we solicit the dimensionless variable \(\xi=\alpha x\) where
\[
\alpha=\left(\frac{m^2\omega^2}{\hbar^2}\right)^{1/4}
=\left(\frac{m\omega}{\hbar}\right)^{1/2}
\]

The Schrödinger equation then becomes
\[
\frac{d^2\psi(\xi)}{d\xi^2}+(\lambda-\xi^2)\psi(\xi)=0
\]

----------------------------------------------------------------

For large \(|\xi|\):
\[
\psi(\xi)=e^{-\xi^2/2}
\qquad
(\text{satisfy the equation as an asymptotic solution, } \xi\to\infty)
\]

For all values of \(\xi\):
\[
\psi(\xi)=e^{-\xi^2/2}H(\xi)
\]
(we obtain for \(H(\xi)\) the differential equation)

\[
\frac{d^2H}{d\xi^2}-2\xi\frac{dH}{d\xi}+(\lambda-1)H=0
\]
which is called the Hermite equation

* Energy levels

Energy spectrum of the LHO:
\[
E_n=\left(n+\frac{1}{2}\right)\hbar\omega
\qquad (n=0,1,2,\ldots)
\]

LHO with lowest state \((n=0)\):
\[
E_0=\frac{1}{2}\hbar\omega \qquad (\text{zero-point energy})
\]

* Hermite polynomials

\[
\psi_n(\xi)=\text{(const)}\ e^{-\xi^2/2}H_n(\xi)
\]

where the function \(H_n(\xi)\) is a polynomial of order \(n\).

Both \(\psi_n(\xi)\) and \(H_n(\xi)\) have parity \((-1)^n\).

The polynomials \(H_n(\xi)\) satisfy the Hermite equation with \(\lambda=2n+1\):
\[
\frac{d^2H_n}{d\xi^2}-2\xi\frac{dH_n}{d\xi}+2nH_n=0
\]

The Hermite polynomials:
\[
H_n(\xi)=(-1)^n e^{\xi^2}\frac{d^n}{d\xi^n}\left(e^{-\xi^2}\right)
\]

The first few Hermite polynomials \(H_n(\xi)\) are:
\[
H_0(\xi)=1
\]
\[
H_1(\xi)=2\xi
\]
\[
H_2(\xi)=4\xi^2-2
\]
\[
H_3(\xi)=8\xi^3-12\xi
\]
\[
H_4(\xi)=16\xi^4-48\xi^2+12
\]
\[
H_5(\xi)=32\xi^5-160\xi^3+120\xi
\]

Orthogonaliteit:
\[
\int_{-\infty}^{\infty}e^{-\xi^2}H_m(\xi)H_n(\xi)\,d\xi=0
\qquad (m\ne n)
\]

Normering:
\[
\int_{-\infty}^{\infty}e^{-\xi^2}H_n^2(\xi)\,d\xi=\sqrt{\pi}\,2^n n!
\]

Energy-eigenfunctions in eigenwaarden van de lineaire harmonische oscillator:

Energy eigenfunctions:
\[
\psi_n(x)=\left(\frac{\alpha}{\sqrt{\pi}\,2^n n!}\right)^{1/2}
e^{-\alpha^2x^2/2}\,H_n(\alpha x)
\qquad
(\alpha=\sqrt{m\omega/\hbar})
\]
\[onleesbaar\ \text{alternatieve vorm met }x_0=\sqrt{\hbar/(m\omega)}\]

- omdat de Hamiltoniaan voor de LHO Hermitisch is geldt:
\[
\langle m|n\rangle=\int_{-\infty}^{\infty}\psi_m^*(x)\psi_n(x)\,dx=\delta_{mn}
\]

- verwachtingswaarde van \(x\):
\[
\langle x\rangle=\int_{-\infty}^{\infty}\psi_n^*(x)\,x\,\psi_n(x)\,dx=0
\]
(zoals \(\psi_n\) even of oneven is)

----------------------------------------------------------------

Bewijs van (4.178)

\[
I=\int_{-\infty}^{\infty}e^{-\xi^2}\,e^{(-s^2+2s\xi)}\,e^{(-t^2+2t\xi)}\,d\xi
\]

\[
=\int_{-\infty}^{\infty}e^{-\xi^2}\,e^{2\xi(s+t)}\,e^{-(s^2+t^2)}\,d\xi
\]

\[
=e^{-(s^2+t^2)}\int_{-\infty}^{\infty}e^{-\left(\xi^2-2\xi(s+t)\right)}\,d\xi
\]

\[
=e^{-(s^2+t^2)}\int_{-\infty}^{\infty}e^{-(\xi-(s+t))^2}\,e^{(s+t)^2}\,d\xi
\]

\[
=e^{-(s^2+t^2)}\,e^{(s+t)^2}\int_{-\infty}^{\infty}e^{-(\xi-(s+t))^2}\,d\xi
\]

zet \(x=\xi-(s+t)\):
\[
I=e^{2st}\int_{-\infty}^{\infty}e^{-x^2}\,dx
\]

Omdat
\[
\int_{-\infty}^{\infty}e^{-x^2}\,dx=\sqrt{\pi}
\]

zeggen we:
\[
I=\sqrt{\pi}\,e^{2st}
\]
\(\square\)
