# REM Examen 2024 — Volledige Oplossingen

---

## Vraag 1 (/5) — Poynting-vector afleiden

**Gevraagd:** Bewijs dat $\mathbf{S} = \frac{1}{\mu_0}\mathbf{E}\times\mathbf{B}$.

**Startpunt:** De arbeid die het elektromagnetisch veld per tijdseenheid per volume verricht op vrije ladingen is $\mathbf{J}\cdot\mathbf{E}$.

Gebruik de **Ampère-Maxwell-wet** om $\mathbf{J}$ te elimineren:

$$
\nabla\times\mathbf{B} = \mu_0\mathbf{J} + \mu_0\varepsilon_0\frac{\partial\mathbf{E}}{\partial t}
\quad\Rightarrow\quad
\mathbf{J} = \frac{1}{\mu_0}(\nabla\times\mathbf{B}) - \varepsilon_0\frac{\partial\mathbf{E}}{\partial t}
$$

Invullen in $\mathbf{J}\cdot\mathbf{E}$:

$$
\mathbf{J}\cdot\mathbf{E} = \frac{1}{\mu_0}\mathbf{E}\cdot(\nabla\times\mathbf{B}) - \varepsilon_0\mathbf{E}\cdot\frac{\partial\mathbf{E}}{\partial t}
$$

Gebruik de  **vectoridentiteit (productformule)** :

$$
\nabla\cdot(\mathbf{E}\times\mathbf{B}) = \mathbf{B}\cdot(\nabla\times\mathbf{E}) - \mathbf{E}\cdot(\nabla\times\mathbf{B})
$$

zodat $\mathbf{E}\cdot(\nabla\times\mathbf{B}) = \mathbf{B}\cdot(\nabla\times\mathbf{E}) - \nabla\cdot(\mathbf{E}\times\mathbf{B})$.

Gebruik  **Faraday** : $\nabla\times\mathbf{E} = -\partial\mathbf{B}/\partial t$, dus:

$$
\mathbf{E}\cdot(\nabla\times\mathbf{B}) = -\mathbf{B}\cdot\frac{\partial\mathbf{B}}{\partial t} - \nabla\cdot(\mathbf{E}\times\mathbf{B})
$$

Alles samengevoegd:

$$
\mathbf{J}\cdot\mathbf{E} = \frac{1}{\mu_0}\left[-\mathbf{B}\cdot\frac{\partial\mathbf{B}}{\partial t} - \nabla\cdot(\mathbf{E}\times\mathbf{B})\right] - \varepsilon_0\mathbf{E}\cdot\frac{\partial\mathbf{E}}{\partial t}
$$

Herken $\mathbf{B}\cdot\partial_t\mathbf{B} = \partial_t(B^2/2)$ en $\mathbf{E}\cdot\partial_t\mathbf{E} = \partial_t(E^2/2)$:

$$
\boxed{-\frac{\partial u}{\partial t} - \nabla\cdot\mathbf{S} = \mathbf{J}\cdot\mathbf{E}}
$$

met de  **energiedichtheid** :

$$
u = \frac{\varepsilon_0 E^2}{2} + \frac{B^2}{2\mu_0}
$$

en de  **Poynting-vector** :

$$
\mathbf{S} = \frac{1}{\mu_0}\mathbf{E}\times\mathbf{B} \qquad \blacksquare
$$

Dit is de  **stelling van Poynting** : de afname van veldenergie plus de energiestroom $\mathbf{S}$ (per eenheidsoppervlak) is gelijk aan de gedissipeerde vermogen.

---

## Vraag 2 (/5) — Condensator met oppervlakstroom

Stel: twee parallelle platen in het $xy$-vlak, gescheiden door afstand $d$ in $\hat{z}$-richting. Bovenplaat: $-\sigma$, $+\mathbf{K}$; onderplaat: $+\sigma$, $-\mathbf{K}$ (met $\mathbf{K} = K\hat{x}$).

### a) Elektrisch veld tussen de platen

Via de wet van Gauss voor twee oneindig uitgestrekte geladen platen ($+\sigma$ en $-\sigma$) geldt:

$$
\mathbf{E} = \frac{\sigma}{\varepsilon_0}\hat{z}
$$

(de velden buiten heffen op, tussen de platen tellen ze op).

### b) Magnetisch veld tussen de platen

Elke plaat met oppervlaktestroom $K\hat{x}$ creëert een veld $\frac{\mu_0 K}{2}$ in de $\hat{y}$-richting. De twee platen (met tegengestelde stromen) zorgen voor:

* **Tussen de platen:** velden tellen op $\Rightarrow \mathbf{B} = \mu_0 K,\hat{y}$
* **Buiten de platen:** velden heffen op $\Rightarrow \mathbf{B} = 0$

### c) Vermogendissipatie via Poynting-flux

De Poynting-vector tussen de platen:

$$
\mathbf{S} = \frac{1}{\mu_0}\mathbf{E}\times\mathbf{B} = \frac{1}{\mu_0}\frac{\sigma}{\varepsilon_0}\hat{z}\times\mu_0 K\hat{y} = -\frac{\sigma K}{\varepsilon_0},\hat{x}
$$

$\mathbf{S}$ wijst naar binnen (in de $-\hat{x}$-richting): energie stroomt de condensator in via de zijkanten.

Beschouw een rechthoekige dwarsdoorsnede van het volume rond de weerstand (breedte $\ell$ langs $y$, hoogte $d$ langs $z$). De vloeistroom door de zijvlakken:

$$
\Phi_S = |\mathbf{S}|\cdot A_{\text{zijkant}} = \frac{\sigma K}{\varepsilon_0}\cdot d,\ell
$$

Het gedissipeerde vermogen in de weerstand:

$$
P = \mathbf{J}\cdot\mathbf{E}\cdot V_{\text{weerstand}} = E \cdot J \cdot \ell \cdot d \cdot w
$$

met $J = K/w$ (stroom per breedte $w$), zodat $P = \frac{\sigma}{\varepsilon_0}\cdot K \cdot d\cdot\ell = \Phi_S$ ✓

De stelling van Poynting is bevestigd: de Poynting-flux = gedissipeerd vermogen.

---

## Vraag 3 (/5) — Ijktransformatie

### a) Bewijs dat $A', V'$ dezelfde $\mathbf{E}$ en $\mathbf{B}$ geven

Gegeven: $\mathbf{A}' = \mathbf{A} + \nabla\lambda$, $V' = V - \partial\lambda/\partial t$.

**Magnetisch veld:**

$$
\mathbf{B}' = \nabla\times\mathbf{A}' = \nabla\times(\mathbf{A}+\nabla\lambda) = \nabla\times\mathbf{A} + \underbrace{\nabla\times(\nabla\lambda)}_{=,0} = \mathbf{B} \quad \checkmark
$$

**Elektrisch veld:**

$$
\mathbf{E}' = -\nabla V' - \frac{\partial\mathbf{A}'}{\partial t}
= -\nabla!\left(V-\frac{\partial\lambda}{\partial t}\right) - \frac{\partial}{\partial t}(\mathbf{A}+\nabla\lambda)
$$

$$
= -\nabla V + \nabla\frac{\partial\lambda}{\partial t} - \frac{\partial\mathbf{A}}{\partial t} - \frac{\partial(\nabla\lambda)}{\partial t}
$$

Omdat $\nabla$ en $\partial/\partial t$ verwisselbaar zijn: $\nabla(\partial\lambda/\partial t) = \partial(\nabla\lambda)/\partial t$, en deze twee termen heffen op:

$$
\mathbf{E}' = -\nabla V - \frac{\partial\mathbf{A}}{\partial t} = \mathbf{E} \quad \checkmark \qquad \blacksquare
$$

### b) Vind $\lambda$ om van $(V=0,;\mathbf{A}\neq 0)$ naar $(V'\neq 0,;\mathbf{A}'=0)$ te gaan

Voor een stationaire lading in de oorsprong zijn de "merkwaardige" potentialen $V=0$ en $\mathbf{A}\neq 0$. Maar $\mathbf{B}=0$ (stationaire lading), dus $\nabla\times\mathbf{A}=0$: $\mathbf{A}$ is een zuivere gradiënt. Het elektrisch veld $\mathbf{E} = -\partial\mathbf{A}/\partial t$, zodat:

$$
\mathbf{A} = -\frac{q}{4\pi\varepsilon_0}\frac{t}{r^2}\hat{r} = \frac{qt}{4\pi\varepsilon_0}\nabla!\left(\frac{1}{r}\right)
$$

We willen $\mathbf{A}' = \mathbf{A} + \nabla\lambda = 0$, dus $\nabla\lambda = -\mathbf{A}$:

$$
\lambda = -\frac{qt}{4\pi\varepsilon_0 r}
$$

Controle op $V'$:

$$
V' = V - \frac{\partial\lambda}{\partial t} = 0 - \left(-\frac{q}{4\pi\varepsilon_0 r}\right) = \frac{q}{4\pi\varepsilon_0 r} \quad \checkmark
$$

Dit is de vertrouwde Coulomb-potentiaal. ✓

---

## Vraag 4 (/5) — Supergeleider en perfecte geleider

### a) Veldlijnen tekenen

**Supergeleider** ($\mathbf{B} = B\hat{z}$ extern):

Meissner-effect: $\mathbf{B} = 0$ binnenin. De veldlijnen worden **volledig omgebogen** rondom de supergeleider (ze kunnen de supergeleider niet binnendringen). De veldlijnendichtheid neemt toe aan de zijkanten.

```
   ↑↑↑  ↑↑↑↑↑  ↑↑↑
   ↑↑↑  ↑↑↑↑↑  ↑↑↑
   ↑↑  ↑       ↑  ↑
   ↑ ↑ [     ] ↑ ↑    ← supergeleider (B = 0 binnenin)
   ↑  ↑        ↑  ↑
   ↑↑↑  ↑↑↑↑↑  ↑↑↑
```

**Perfecte geleider** ($\mathbf{E} = E\hat{z}$ extern):

$\mathbf{E} = 0$ binnenin. Veldlijnen staan **loodrecht op het oppervlak** en eindigen op oppervlakteladingen. Geen veld in het inwendige.

```
   ↑↑↑↑↑↑↑↑↑↑↑↑↑
   ↑  ↑↑↑↑↑↑↑  ↑
   ↑  [       ] ↑    ← perfecte geleider (E = 0 binnenin, +σ boven, −σ onder)
   ↑  ↑↑↑↑↑↑↑  ↑
   ↑↑↑↑↑↑↑↑↑↑↑↑↑
```

### b) Stroom- en ladingsverdelingen

**Supergeleider** — oppervlaktestroom $\mathbf{K}$:

Beschouw een plaat-supergeleider met normale richting $\hat{x}$ en extern $\mathbf{B} = B\hat{z}$ (parallel aan het oppervlak).

Randvoorwaarde voor $\mathbf{B}$: $(\mathbf{B} *{\text{buiten}} - \mathbf{B}* {\text{binnen}})\times\hat{n} = \mu_0\mathbf{K}$.

$$
B\hat{z}\times\hat{x} = \mu_0\mathbf{K}
\quad\Rightarrow\quad
\mathbf{K} = \frac{B}{\mu_0}\hat{y}
$$

op het rechtervlak, en $\mathbf{K} = -\frac{B}{\mu_0}\hat{y}$ op het linkervlak. Deze Meissner-stroom compenseert het uitwendige veld volledig.

---

**Perfecte geleider** — oppervlaktelading $\sigma$:

Randvoorwaarde voor $\mathbf{E}$: $(\mathbf{E} *{\text{buiten}} - \mathbf{E}* {\text{binnen}})\cdot\hat{n} = \sigma/\varepsilon_0$.

Bovenvlak ($\hat{n} = +\hat{z}$):

$$
\sigma_{\text{boven}} = \varepsilon_0(E - 0) = \varepsilon_0 E > 0
$$

Ondervlak ($\hat{n} = -\hat{z}$):

$$
\sigma_{\text{onder}} = -\varepsilon_0 E < 0
$$

De positieve ladingen zitten op het  **bovenvlak** , de negatieve op het  **ondervlak** . Randeffecten worden genegeerd.

---

## Vraag 5 (/5) — Lorentztransformatie 4-impulsvector & CM-frame

### Lorentztransformatie van de 4-impulsvector

De 4-impuls is $p^\mu = (E/c,, p_x,, p_y,, p_z)$. Bij een Lorentz-boost in de $x$-richting met snelheid $v$ (en $\gamma = 1/\sqrt{1-\beta^2}$, $\beta = v/c$):

$$
\boxed{
\begin{aligned}
E' &= \gamma(E - v,p_x)\
p_x' &= \gamma!\left(p_x - \frac{v,E}{c^2}\right)\
p_y' &= p_y\
p_z' &= p_z
\end{aligned}
}
$$

### Snelheid van het CM-frame

Het middelpunt-van-massa (CM)-frame is het frame waarvoor de **totale impuls nul** is.

Eis: $P'_x = 0$:

$$
\gamma!\left(\sum_i p_i - \frac{v_{\text{CM}}}{c^2}\sum_i E_i\right) = 0
$$

$$
\boxed{v_{\text{CM}} = \frac{c^2,\sum_i p_i}{\sum_i E_i}}
$$

Dit is een lorentzinvariant resultaat: de snelheid van het CM-frame is de verhouding van de totale (ruimtelijke) impuls tot de totale energie (gedeeld door $c^2$). Merk op dat dit de relativistische veralgemening is van $v_{\text{CM}} = \sum m_i v_i / \sum m_i$.

---

## Vraag 6 (/5) — Relativistische botsing

**Gegeven:** Bewegend deeltje (massa $m$) met energie $E_1$, stationair deeltje (massa $m$) met energie $E_2 = mc^2$, en $E_1 = 2E_2 = 2mc^2$.

**Impuls van het bewegende deeltje:**

$$
E_1^2 = (p_1 c)^2 + (mc^2)^2
\quad\Rightarrow\quad
(2mc^2)^2 - (mc^2)^2 = (p_1 c)^2
\quad\Rightarrow\quad
p_1 = \sqrt{3},mc
$$

**Behoud van 4-impuls** (perfecte inelastische botsing):

| Grootheid | Vóór                            | Na      |
| --------- | --------------------------------- | ------- |
| Energie   | $2mc^2 + mc^2 = 3mc^2$          | $E_f$ |
| Impuls    | $\sqrt{3},mc + 0 = \sqrt{3},mc$ | $P_f$ |

dus $E_f = 3mc^2$ en $P_f = \sqrt{3},mc$.

**Invariante massa van het conglomeraat:**

$$
M^2c^4 = E_f^2 - (P_f c)^2 = (3mc^2)^2 - (\sqrt{3},mc\cdot c)^2 = 9m^2c^4 - 3m^2c^4
$$

$$
\boxed{M = m\sqrt{6} \approx 2{,}45,m}
$$

**Snelheid:**

$$
v_f = \frac{P_f c^2}{E_f} = \frac{\sqrt{3},mc\cdot c^2}{3mc^2} = \frac{c}{\sqrt{3}}
$$

$$
\boxed{v_f = \frac{c}{\sqrt{3}} \approx 0{,}577,c}
$$

*Controle:* $\gamma = 1/\sqrt{1-1/3} = \sqrt{3/2}$, en $\gamma M c^2 = \sqrt{3/2}\cdot m\sqrt{6}\cdot c^2 = 3mc^2$ ✓

---

## Vraag 7 (/5) — Methode van beeldladingen (60°-wig)

### a) Aantal beeldladingen en hun posities

Twee geaarde geleidende halfvlakken op $\theta=0$ en $\theta=60°$. Puntlading $+q$ op $(r,;\theta_0=30°)$.

Voor een hoek $\alpha = \pi/n$ (met $n$ een geheel getal) zijn er precies $2n-1$ beeldladingen nodig. Hier: $\alpha = 60° = \pi/3$, dus $n=3$, en er zijn $\boxed{5}$ **beeldladingen** nodig.

De beelden worden systematisch gegenereerd door herhaalde spiegeling in beide vlakken (de symmetriegroep is de diedraalgroep $D_3$ van orde 6):

| Positie                         | Lading | Verkregen door                                       |
| ------------------------------- | ------ | ---------------------------------------------------- |
| $(r,;-30°)$ = $(r,;330°)$ | $-q$ | spiegeling van$+q$ in $\theta=0$                 |
| $(r,;90°)$                   | $-q$ | spiegeling van$+q$ in $\theta=60°$              |
| $(r,;150°)$                  | $+q$ | spiegeling van$-q$ op $330°$ in $\theta=60°$ |
| $(r,;270°)$                  | $+q$ | spiegeling van$-q$ op $90°$ in $\theta=0$     |
| $(r,;210°)$                  | $-q$ | spiegeling van$+q$ op $150°$ in $\theta=0$    |

*Controle:* Op $\theta=0$ zijn telkens paren $(\pm q)$ spiegelbeelden t.o.v. dat vlak → $V=0$ ✓

Op $\theta=60°$: idem ✓

### b) Werkt de methode voor willekeurige $0°<\theta<60°$?

**Ja.** De methode van beeldladingen werkt voor een wig met hoek $\alpha = \pi/n$ ($n$ geheel getal) voor **elke** positie van de lading binnen de wig. De beeldladingen worden geplaatst op $\pm\theta_0 + 2k\cdot60°$ voor $k=0,1,2$ (zonder de oorspronkelijke lading), en liggen altijd buiten de wig (voor $0°<\theta_0<60°$). De randvoorwaarden $V=0$ op beide vlakken zijn voor elke $\theta_0$ bevredigd door dezelfde spiegelingssymmetrie.

---

## Vraag 8 (/5) — Vreemde krachtwet in tensornotatie

De "vreemde" krachtwet is waarschijnlijk de **covariante Lorentzkracht** in tensornotatie:

$$
f^\mu = q,F^{\mu\nu},u_\nu
$$

waarbij $F^{\mu\nu}$ de veldsterktetensor is en $u^\mu = \gamma(c, v_x, v_y, v_z)$ de 4-snelheid.

Met de $(+,-,-,-)$ metriek en $u_\nu = \gamma(c, -v_x, -v_y, -v_z)$:

$$
F^{\mu\nu} = \begin{pmatrix} 0 & -E_x/c & -E_y/c & -E_z/c \ E_x/c & 0 & -B_z & B_y \ E_y/c & B_z & 0 & -B_x \ E_z/c & -B_y & B_x & 0 \end{pmatrix}
$$

### a) $\mu = 1$ (x-component)

$$
f^1 = q,(F^{10}u_0 + F^{11}u_1 + F^{12}u_2 + F^{13}u_3)
$$

$$
= q!\left(\frac{E_x}{c}\cdot\gamma c + 0 + (-B_z)(-\gamma v_y) + B_y(-\gamma v_z)\right)
$$

$$
= q\gamma\left(E_x + v_y B_z - v_z B_y\right) = q\gamma,(\mathbf{E}+\mathbf{v}\times\mathbf{B})_x
$$

### b) $\mu = 2$ en $\mu = 3$

Door hetzelfde te doen:

$$
f^2 = q\gamma,(E_y + v_z B_x - v_x B_z) = q\gamma,(\mathbf{E}+\mathbf{v}\times\mathbf{B})_y
$$

$$
f^3 = q\gamma,(E_z + v_x B_y - v_y B_x) = q\gamma,(\mathbf{E}+\mathbf{v}\times\mathbf{B})_z
$$

**In vectornotatie** (ruimtelijk deel):

$$
\mathbf{f} = q\gamma(\mathbf{E}+\mathbf{v}\times\mathbf{B}) = \frac{d\mathbf{p}}{dt}
$$

Dit is de relativistische Lorentzkracht ✓.

### c) Interpretatie van $\tilde{q}$

Als de krachtwet de **duale veldsterktetensor** $\tilde{F}^{\mu\nu} = \tfrac{1}{2}\varepsilon^{\mu\nu\alpha\beta}F_{\alpha\beta}$ bevat:

$$
f^\mu = \tilde{q},\tilde{F}^{\mu\nu},u_\nu
$$

dan wisselen $\mathbf{E}\leftrightarrow c\mathbf{B}$ (en $c\mathbf{B}\leftrightarrow -\mathbf{E}$) in de componenten. De ruimtelijke kracht wordt dan $\tilde{q},\gamma(c^2\mathbf{B} - \mathbf{v}\times\mathbf{E}/c^2)$... wat overeenkomt met de **kracht op een magnetische monopole** (magnetische lading). $\tilde{q}$ is dus de **magnetische lading** (magnetisch monopolemoment). Dit is de duale symmetrie van de Maxwell-vergelijkingen onder $\mathbf{E}\to c\mathbf{B}$, $c\mathbf{B}\to -\mathbf{E}$.

---

## Vraag 9 (/5) — Coaxiale kabel met magnetisch materiaal

**Gegeven:** Coaxiale kabel met binnenstraal $a$, buitenstraal $b$, stroom $I$ in de binnenste geleider. Het materiaal tussen de geleiders heeft magnetische susceptibiliteit $\chi_m$ (en dus relatieve permeabiliteit $\mu_r = 1+\chi_m$).

**Methode:** Ampère's wet in termen van $\mathbf{H}$ (met alleen vrije stromen):

$$
\oint \mathbf{H}\cdot d\boldsymbol{\ell} = I_{\text{vrij}}
$$

Voor een cirkelvormig Ampèriaans pad op straal $r$ ($a < r < b$):

$$
H\cdot 2\pi r = I
\quad\Rightarrow\quad
\mathbf{H} = \frac{I}{2\pi r},\hat{\varphi}
$$

Het magnetisch veld $\mathbf{B} = \mu\mathbf{H} = \mu_0\mu_r\mathbf{H}$:

$$
\boxed{\mathbf{B} = \frac{\mu_0(1+\chi_m),I}{2\pi r},\hat{\varphi}
\qquad (a < r < b)}
$$

Buiten ($r > b$) en binnenin de geleiders ($r < a$): $\mathbf{B}$ volgt via dezelfde wet (zonder het materiaal).

---

## Vraag 10 (/5) — Lorentztransformatie & CM-frame (variant)

### a) Lorentztransformatie 4-energie-impulsvector (boost langs $x$)

Zie ook Vraag 5. Bij een boost $v$ in de $x$-richting:

# $$\boxed{

\begin{pmatrix} E'/c \ p_x' \ p_y' \ p_z' \end{pmatrix}

\begin{pmatrix} \gamma & -\gamma\beta & 0 & 0\ -\gamma\beta & \gamma & 0 & 0\ 0 & 0 & 1 & 0\ 0 & 0 & 0 & 1 \end{pmatrix}
\begin{pmatrix} E/c \ p_x \ p_y \ p_z \end{pmatrix}
}$$

### b) Snelheid van het nul-impuls-frame

$n$ deeltjes bewegen alle naar rechts met snelheden $v_1, v_2, \ldots$, massa's $m_i$, waarvoor:

$$
E_i = \gamma_i m_i c^2, \qquad p_i = \gamma_i m_i v_i
$$

In het CM-frame geldt $\sum_i p_i' = 0$. Uit de Lorentztransformatie:

$$
\sum_i p_i' = \gamma_{\text{CM}}!\left(\sum_i p_i - \frac{v_{\text{CM}}}{c^2}\sum_i E_i\right) = 0
$$

$$
\boxed{v_{\text{CM}} = \frac{c^2\sum_i \gamma_i m_i v_i}{\sum_i \gamma_i m_i c^2} = \frac{\sum_i \gamma_i m_i v_i}{\sum_i \gamma_i m_i}}
$$

Dit is de relativistische veralgemening van de klassieke CM-snelheid $v_{\text{CM}} = \sum m_i v_i / \sum m_i$.

---

## Vraag 11 (/5) — Gemiddelde intensiteit van een stralende bron

**Gegeven:** $\mathbf{E}$ en $\mathbf{B}$ bevatten termen $\propto r^{-1},,r^{-2},,r^{-3}$ met een $\theta$-afhankelijkheid.

**Sleutelprincipe:** In het verre veld ($r\to\infty$) dragen alleen de $r^{-1}$-termen bij tot de stralingsintensiteit — hogere machten vallen te snel af.

Stel dat de $r^{-1}$-termen de vorm hebben (elektrische dipool als voorbeeld):

$$
\mathbf{E} \approx \frac{A\sin\theta}{r}\cos(\omega t - kr),\hat{\theta}, \qquad
\mathbf{B} = \frac{\mathbf{E}}{c}
$$

**Tijdsgemiddelde Poynting-vector** (gemiddeld over één periode):

$$
\langle\mathbf{S}\rangle = \frac{1}{2\mu_0 c}|E_0(\theta)|^2\frac{1}{r^2},\hat{r}
$$

**Totaal uitgestraald vermogen** (integreer over een bol van straal $r$):

$$
P = \oint \langle\mathbf{S}\rangle\cdot d\mathbf{A} = \frac{|A|^2}{2\mu_0 c}\int_0^\pi\int_0^{2\pi}\sin^2\theta,r^2\sin\theta,d\theta,d\varphi
$$

$$
= \frac{|A|^2}{2\mu_0 c}\cdot 2\pi\int_0^\pi\sin^3\theta,d\theta = \frac{|A|^2}{2\mu_0 c}\cdot 2\pi\cdot\frac{4}{3} = \frac{4\pi|A|^2}{3\mu_0 c}
$$

Dit is onafhankelijk van $r$ (energiebehoud) ✓. In de praktijk: identificeer $A$ uit de gegeven $r^{-1}$-termen en reken P uit.

---

## Vraag 12 (/5) — Potentiaal op een bol

**Gegeven:** $V(R,\theta) = V_0\frac{1}{2}(3\cos^2\theta - 1) = V_0,P_2(\cos\theta)$

### a) Potentiaal binnen en buiten de bol

Los de Laplace-vergelijking op in bolvormige coördinaten. De algemene oplossing:

$$
V(r,\theta) = \sum_\ell \left(A_\ell r^\ell + B_\ell r^{-(\ell+1)}\right)P_\ell(\cos\theta)
$$

* **Binnen** ($r<R$): regulariteit in $r=0$ vereist $B_\ell = 0$
* **Buiten** ($r>R$): $V\to 0$ voor $r\to\infty$ vereist $A_\ell = 0$

Grensvoorwaarde $V(R,\theta) = V_0,P_2(\cos\theta)$: enkel $\ell=2$ draagt bij.

$$
\textbf{Binnen:}\quad V_{\text{in}}(r,\theta) = V_0\left(\frac{r}{R}\right)^{!2}P_2(\cos\theta) = \frac{V_0}{2}\frac{r^2}{R^2}(3\cos^2\theta - 1)
$$

$$
\textbf{Buiten:}\quad V_{\text{out}}(r,\theta) = V_0\left(\frac{R}{r}\right)^{!3}P_2(\cos\theta) = \frac{V_0}{2}\frac{R^3}{r^3}(3\cos^2\theta - 1)
$$

*Controle:* Bij $r=R$: $V_{\text{in}} = V_{\text{out}} = V_0 P_2(\cos\theta)$ ✓

### b) Welke multipool domineert in het verre veld?

In de multipool-ontwikkeling geldt $V\sim r^{-(\ell+1)}$:

| $\ell$ | Naam       | $V\sim$  |
| -------- | ---------- | ---------- |
| 0        | monopool   | $r^{-1}$ |
| 1        | dipool     | $r^{-2}$ |
| 2        | kwadrupool | $r^{-3}$ |
| 3        | oktupool   | $r^{-4}$ |

Het buitenpotentiaal gaat als $r^{-3}$, dus de leidende bijdrage is een **kwadrupool** ($\ell=2$).

---

## Bonusvraag 1 (/2) — Eigentime van een foton

**Vraag:** Hoe lang doet een foton er in zijn *eigen* referentiekader over om de 4,7 lichtjaar naar Proxima Centauri af te leggen?

**Antwoord:** $\boxed{\tau = 0}$

Fotonen bewegen langs **nulgeodeten** (lichtkegels): $ds^2 = c^2 dt^2 - dx^2 = 0$.

De eigentime (proper time) is:

$$
d\tau^2 = dt^2 - \frac{dx^2}{c^2} = dt^2!\left(1-\frac{v^2}{c^2}\right) = 0 \quad\text{voor } v = c
$$

Een foton ervaart geen verstrijkende tijd. Er bestaat bovendien geen inertiaalstelsel dat meebeweegt met een foton (dat zou een Lorentz-boost vereisen met $\beta = 1$, waarvoor $\gamma \to \infty$). Voor het foton zijn begin en einde van de reis **gelijktijdig** — afstand en tijd zijn beide nul.

---

## Bonusvraag 2 (/2) — Susceptibiliteiten

**Elektrische susceptibiliteit van een perfecte geleider:**

Binnenin een perfecte geleider is $\mathbf{E} = 0$ (vrije ladingen schermen het veld volledig af). Voor het diëlektrisch model geldt $\mathbf{D} = \varepsilon_0(1+\chi_e)\mathbf{E}$. Omdat de vrije ladingen elk uitwendig veld compenseren, gedraagt de geleider zich alsof $\varepsilon_r\to\infty$, en dus:

$$
\boxed{\chi_e \to \infty}
$$

**Magnetische susceptibiliteit van een supergeleider:**

Het Meissner-effect zorgt voor $\mathbf{B} = 0$ binnenin de supergeleider, ook als het uitwendige veld niet nul is. Met $\mathbf{B} = \mu_0(\mathbf{H}+\mathbf{M}) = \mu_0(1+\chi_m)\mathbf{H} = 0$:

$$
\boxed{\chi_m = -1}
$$

De supergeleider is een  **perfecte diamagneet** : de magnetisatie $\mathbf{M} = -\mathbf{H}$ compenseert het uitwendige veld volledig. Dit is fundamenteel anders dan gewone diamagneten (waarbij $\chi_m \approx -10^{-5}$).

---

*Oplossingen opgesteld ter voorbereiding van het REM-examen 2024.*
