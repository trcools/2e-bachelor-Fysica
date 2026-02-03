# HC15 — Chemisch evenwicht (thermodynamische afleiding)

## 1) Kernprincipe (PT-controle)
Bij constante druk en temperatuur wil een systeem zijn Gibbs vrije energie minimaliseren:
$$
dG = \sum_i \mu_i\,dn_i \qquad (T,p\ \text{constant})
$$
Dus: **de enige “knop”** om $G$ te verlagen is de **samenstelling** $n_i$.

---

## 2) Fysisch evenwicht (fase-evenwicht) als prototype
Voor een zuivere stof $A$ die kan uitwisselen tussen vloeistof en gas:
$$
A(l)\rightleftharpoons A(g)
$$
Evenwicht wanneer uitwisseling $dn$ geen daling van $G$ meer geeft:
$$
\mu_A^{(l)}=\mu_A^{(g)}.
$$

### Ideaal gas + zuivere vloeistof
- Zuivere vloeistof: drukafhankelijkheid vaak verwaarloosbaar (klein molair volume) → neem $\mu_A^{(l)}\approx \mu_{A,l}^\circ(T)$.
- Ideaal gas: 
$$
\mu_A^{(g)}=\mu_{A,g}^\circ(T)+RT\ln\!\left(\frac{p_A}{p^\circ}\right).
$$
Dan volgt voor de dampdruk:
$$
\frac{p_A^{eq}}{p^\circ}=\exp\!\left(-\frac{\Delta G_{vap}^\circ(T)}{RT}\right),
\quad 
\Delta G_{vap}^\circ=\mu_{A,g}^\circ-\mu_{A,l}^\circ.
$$

**Temperatuurcorrectie (als tabellen bij 298 K staan):**
$$
\Delta G^\circ(T)=\Delta H^\circ(T)-T\Delta S^\circ(T),
$$
en vaak (over tientallen K) nemen we bij benadering:
$$
\Delta H^\circ(T)\approx \Delta H^\circ(298),\qquad
\Delta S^\circ(T)\approx \Delta S^\circ(298),
$$
maar **nooit** $\Delta G^\circ$ als constant nemen omdat $T\Delta S^\circ$ meespeelt.

---

## 3) Dynamisch evenwicht (interpretatie)
Evenwicht betekent niet “niets gebeurt”, maar:
- heen- en terugproces lopen nog,
- **flux heen = flux terug**,
- macroscopisch stationair (tijdonafhankelijk).

Dit geldt zowel voor fase-evenwichten als chemische reacties.

---

## 4) Chemisch evenwicht: reactie-extensie $\xi$ (vorderingsgraad)
Voor een algemene reactie:
$$
\sum_r \nu_r R_r \rightleftharpoons \sum_p \nu_p P_p
$$
Definieer reactie-extensie $\xi$ zodat:
$$
dn_i = \nu_i\,d\xi
$$
met $\nu_i>0$ voor producten en $\nu_i<0$ voor reagentia (tekenconventie).

Dan:
$$
dG=\sum_i \mu_i\,dn_i=\left(\sum_i \nu_i\mu_i\right)\,d\xi.
$$

### Reactie-Gibbsenergie
Definieer:
$$
\Delta_r G \equiv \sum_i \nu_i\mu_i.
$$
Dan:
- $\Delta_r G=0$  → evenwicht
- $\Delta_r G<0$ → reactie gaat **voorwaarts** (toename $\xi$ verlaagt $G$)
- $\Delta_r G>0$ → reactie gaat **achterwaarts**

Grafisch: $G(\xi)$ heeft een minimum waar $\Delta_r G=0$.

**Waarom is $G(\xi)$ soms krom en soms lineair?**
- Als $\mu_i$ afhangt van samenstelling (bv. gassen/opl.) → $G(\xi)$ typisch **krom** en minimum kan bij $\xi\neq 0,1$ liggen.
- Als alle betrokken stoffen zuiver vast/vloeibaar zijn ($a=1$) → $\mu_i$ ~ constant → $G(\xi)$ **lineair** → minimum op een uiteinde (alles reagentia of alles producten).

---

## 5) Chemische potentiaal en activiteit
Algemene splitsing:
$$
\mu_i=\mu_i^\circ(T)+RT\ln a_i.
$$

Invullen in $\Delta_r G$ geeft:
$$
\Delta_r G=\Delta_r G^\circ + RT\ln Q,
$$
waar
$$
\Delta_r G^\circ=\sum_i \nu_i\mu_i^\circ,
\qquad
Q=\frac{\prod_p a_p^{\nu_p}}{\prod_r a_r^{\nu_r}}
$$
het **reactiequotiënt** is.

Evenwicht betekent $\Delta_r G=0$, dus:
$$
0=\Delta_r G^\circ + RT\ln K
\quad\Rightarrow\quad
K=\exp\!\left(-\frac{\Delta_r G^\circ}{RT}\right),
$$
en bij evenwicht:
$$
Q=K.
$$

Interpretatie:
- $Q<K$ → voorwaarts
- $Q>K$ → achterwaarts
- $Q=K$ → evenwicht

---

## 6) Gasevenwichten: twee “modellen” (consistent blijven!)

### (A) Partiële-druk model ($K_p^\circ$)
Voor ideale gassen:
$$
a_i=\frac{p_i}{p^\circ}.
$$
Dan:
$$
K_p^\circ=\prod_i \left(\frac{p_i}{p^\circ}\right)^{\nu_i}
=\exp\!\left(-\frac{\Delta_r G^\circ}{RT}\right).
$$
Belangrijk: $K_p^\circ$ is **druk-onafhankelijk** (standaardtoestand $p^\circ=1$ bar zit erin).

### (B) Molfractie model ($K_x$)
Gebruik $p_i = x_i\,p_{tot}$ en schrijf:
$$
Q=\frac{\prod_p x_p^{\nu_p}}{\prod_r x_r^{\nu_r}}
\left(\frac{p_{tot}}{p^\circ}\right)^{\Delta\nu},
\quad 
\Delta\nu=\sum_p \nu_p-\sum_r \nu_r.
$$
Dus:
$$
K_p^\circ = K_x\left(\frac{p_{tot}}{p^\circ}\right)^{\Delta\nu}
\quad\Rightarrow\quad
K_x = K_p^\circ\left(\frac{p_{tot}}{p^\circ}\right)^{-\Delta\nu}.
$$

**Vuistregel:** kies $Q$ in termen van $p_i$ of in termen van $x_i$, maar gebruik dan ook de bijhorende constante.

---

## 7) Oplosstrategie met een massa-balans (en $\xi$ of $\alpha$)
1. Schrijf de reactie en $Q$ (met juiste machten $\nu_i$).
2. Introduceer reactie-extensie (in HC15 vaak “$\alpha$”):
   - $n_i = n_{i,0}+\nu_i\alpha$
3. Bereken molfracties:
$$
x_i=\frac{n_i}{\sum_j n_j}
$$
4. Vul in in $Q(\alpha)$ en los $Q(\alpha)=K$ op.
5. Gebruik $\alpha$ om alle $x_i$ of $p_i$ te bepalen.

Praktisch: vaak kan je met limieten werken:
- als $K$ extreem klein → productfractie $\approx 0$
- als $K$ extreem groot → reagentiafractie $\approx 0$
Dat is ook een snelle check tegen rekenfouten (bv. “factor 1000” door kJ i.p.v. J).

---

## 8) Temperatuurafhankelijkheid (Gibbs–Helmholtz / van ’t Hoff)
Uit
$$
K=\exp\!\left(-\frac{\Delta_r G^\circ}{RT}\right)
$$
volgt bij vaste druk:
$$
\frac{\partial \ln K}{\partial T}=\frac{\Delta_r H^\circ}{RT^2}.
$$
Dus:
- $\Delta_r H^\circ>0$ (endotherm) → $K$ **stijgt** met $T$ → meer producten bij hogere $T$.
- $\Delta_r H^\circ<0$ (exotherm) → $K$ **daalt** met $T$ → minder producten bij hogere $T$.

Dit is de nette, foutbestendige versie van “Le Châtelier”.

---

## 9) Drukafhankelijkheid (via $K_x$ en $\Delta\nu$)
Hoewel $K_p^\circ$ druk-onafhankelijk is, verschuift het evenwicht met $p_{tot}$ via:
$$
K_x = K_p^\circ\left(\frac{p_{tot}}{p^\circ}\right)^{-\Delta\nu}.
$$
- $\Delta\nu<0$ (minder mol gas aan productzijde) → hogere druk maakt $K_x$ groter → evenwicht naar **producten**.
- $\Delta\nu>0$ → hogere druk duwt naar **reagentia**.

Ook hier: dit is Le Châtelier, maar dan in rekenvorm.

---

## 10) Voorbeelden uit HC15 (conceptueel)
### Ozonvorming (typisch heel klein bij 298 K)
Afhankelijk van gekozen stoichiometrie verandert $K$ mee (kwadrateren, wortel nemen), maar **evenwichts-samenstelling blijft dezelfde** zolang je consistent blijft.

### Haber–Bosch (NH₃)
Thermodynamisch gunstig bij lage $T$ (grote $K$), maar kinetisch traag (N≡N sterk) → katalysator nodig.
Hogere $T$ versnelt kinetiek maar werkt tegen evenwicht als de reactie exotherm is; daarom combineer je met hoge druk als $\Delta\nu<0$.

---

## 11) Examenvallen (HC15-klassiekers)
- **kJ vs J:** in $K=\exp(-\Delta G^\circ/RT)$ moet $\Delta G^\circ$ in **J/mol**.
- $K_p^\circ$ is **niet** hetzelfde als $K_x$ als $\Delta\nu\neq 0$ en $p_{tot}\neq p^\circ$.
- Grote/kleine $K$ geeft vaak sterke intuïtie, maar let op: standaardtoestanden beïnvloeden absolute grootte; veranderingen van $K$ met $T$ zijn het betrouwbaarst.
