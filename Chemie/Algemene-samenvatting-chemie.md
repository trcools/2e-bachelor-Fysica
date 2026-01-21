# HC13–HC14 — Thermodynamica: F, G, μ, activiteiten, evenwicht (samenvatting)

## 1) Hoofdwetten + evolutieprincipe
**1e hoofdwet (conventie les):** warmte en arbeid positief als ze $U_{\text{sys}}$ doen stijgen
$$
dU = \delta q + \delta w
$$

**2e hoofdwet (reversibel):**
$$
dS = \frac{\delta q_{\text{rev}}}{T}
$$

**Evolutie/evenwicht (systeem + omgeving):**
$$
\Delta S_{\text{tot}}=\Delta S_{\text{sys}}+\Delta S_{\text{omg}} \ge 0
$$
Probleem: je moet de omgeving mee modelleren → opgelost door toestandfuncties $F$ en $G$.

---

## 2) Helmholtz vrije energie $F$ (of $A$)
**Definitie**
$$
F \equiv U - TS
$$

**Differentiaal**
$$
dF = dU - T\,dS - S\,dT
$$

**Temperatuurcontrole ($dT=0$):**
$$
dF = dU - T\,dS
$$

**Bij $T,V$-controle + enkel $pV$-arbeid:**
$$
\Delta F \le 0
$$
**Evenwicht:** $F$ minimaal (bij vaste $T,V$).

**Werk-interpretatie (reversibel = “gelijkheidsdenken”):** $\Delta F$ geeft de arbeidsgrens (maximaal haalbare arbeid onder de gekozen randvoorwaarden).

---

## 3) Gibbs vrije enthalpie $G$
**Definitie**
$$
G \equiv H-TS = U+PV-TS
$$

**Algemene differentiaal**
$$
dG = dU + P\,dV + V\,dP - T\,dS - S\,dT
$$

**Voor eenvoudig $pV$-werksysteem (vaste samenstelling):**
$$
dG = -S\,dT + V\,dP
$$

**Bij $T,P$-controle + enkel $pV$-arbeid:**
$$
\Delta G \le 0
$$
**Evenwicht:** $G$ minimaal (bij vaste $T,P$).

**Werkbetekenis (PT):** bij $T,P$ constant is $-\Delta G$ de maximaal haalbare *niet-$pV$* arbeid (bv. elektrische) voor een reversibel proces.

---

## 4) Afgeleiden: T- en P-afhankelijkheid van $G$
Uit
$$
dG=-S\,dT+V\,dP
$$
volgt
$$
\left(\frac{\partial G}{\partial T}\right)_P=-S,\qquad
\left(\frac{\partial G}{\partial P}\right)_T=V
$$
Dus: $S>0$ ⇒ $G$ daalt met stijgende $T$.

---

## 5) Drukafhankelijkheid: ideaal gas vs. gecondenseerde fasen
### Ideaal gas (molaire vorm)
$$
\left(\frac{\partial G_m}{\partial P}\right)_T = V_m = \frac{RT}{P}
$$
Integratie:
$$
G_m(T,P)=G_m^\circ(T)+RT\ln\!\left(\frac{P}{P^\circ}\right)
$$

**Ideaal gasmengsel:** vervang $P$ door partiële druk $P_i$:
$$
G_{m,i}=\mu_i=\mu_i^\circ(T)+RT\ln\!\left(\frac{P_i}{P^\circ}\right)
$$

### Gecondenseerde fasen (vloeistof/vaste stof)
$V_m$ klein ⇒ $V\,\Delta P$ klein ⇒ vaak:
$$
G_m(T,P)\approx G_m^\circ(T)
$$
(drukafhankelijkheid verwaarloosbaar in veel problemen)

---

## 6) Samenstelling → chemische potentiaal $\mu_i$
Bij $T,P$ vast is $G$ homogeen van graad 1 in $n_i$ (verdubbel alles ⇒ $G$ verdubbelt).

**Definitie**
$$
\mu_i \equiv \left(\frac{\partial G}{\partial n_i}\right)_{T,P,n_{j\neq i}}
$$

**Euler (homogeniteit):**
$$
G=\sum_i n_i\mu_i
$$

**Volledige differentiaal (met samenstelling):**
$$
dG=-S\,dT+V\,dP+\sum_i \mu_i\,dn_i
$$

---

## 7) Gibbs–Duhem relatie
Start met:
$$
G=\sum_i n_i\mu_i \Rightarrow dG=\sum_i\mu_i\,dn_i+\sum_i n_i\,d\mu_i
$$
Vergelijk met
$$
dG=-S\,dT+V\,dP+\sum_i\mu_i\,dn_i
$$
$\implies$
$$
\sum_i n_i\,d\mu_i=-S\,dT+V\,dP
$$

**Bij vaste $T,P$:**
$$
\sum_i n_i\,d\mu_i=0
$$
Dus $\mu_i$ zijn niet onafhankelijk.

---

## 8) Activiteiten: universele vorm voor $\mu$
$$
\mu_i=\mu_i^\circ+RT\ln a_i
$$
waar $a_i$ dimensieloos is.

**Gekende activiteiten:**
- **Zuivere gecondenseerde fase:** $a_i=1 \Rightarrow \mu_i=\mu_i^\circ$
- **Ideaal gas (zuiver of mengsel):** $a_i=P_i/P^\circ$

**Valkuil:** in $\ln$ mag enkel dimensieloos:
$$
\ln\!\left(\frac{P}{P^\circ}\right)\ \text{niet}\ \ln(P)
$$

---

## 9) Fase-evenwicht = gelijkheid van chemische potentialen
Voor component $A$ in twee fasen $\alpha,\beta$ bij $T,P$ constant:
$$
\mu_A^{(\alpha)}=\mu_A^{(\beta)}
$$
Materie gaat spontaan naar de fase met lagere $\mu$ tot gelijkheid.

---

## 10) Dampdruk (voorbeeld water) via $\mu$ en activiteit
Evenwicht vloeistof ↔ gas:
$$
\mu_g=\mu_l
$$
Invullen:
$$
\mu_g^\circ+RT\ln\!\left(\frac{P^*}{P^\circ}\right)=\mu_l^\circ
$$
Dus
$$
\ln\!\left(\frac{P^*}{P^\circ}\right)=-\frac{\mu_g^\circ-\mu_l^\circ}{RT}
=-\frac{\Delta G_{\text{vap}}^\circ}{RT}
$$
$$
P^*=P^\circ\exp\!\left(-\frac{\Delta G_{\text{vap}}^\circ}{RT}\right)
$$

---

## 11) Temperatuursafhankelijkheid van evenwichten
$$
\Delta G^\circ(T)=\Delta H^\circ(T)-T\Delta S^\circ(T)
$$

Via Hess-cyclus met opwarmen/afkoelen:
$$
\frac{d(\Delta H^\circ)}{dT}=\Delta C_p,\qquad
\frac{d(\Delta S^\circ)}{dT}=\frac{\Delta C_p}{T}
$$

Vaak eerste benadering (beperkt $T$-interval): $\Delta H^\circ,\Delta S^\circ$ ongeveer constant.

**Kookpunt bij 1 bar (evenwicht vloeistof–gas):**
$$
\Delta G^\circ(T_b)=0 \Rightarrow T_b\approx \frac{\Delta H^\circ}{\Delta S^\circ}
$$

---

## Mini-cheatsheet (must know)
- $F=U-TS$, bij $T,V$ const + enkel $pV$: $\Delta F\le 0$
- $G=H-TS=U+PV-TS$, bij $T,P$ const + enkel $pV$: $\Delta G\le 0$
- $dG=-S\,dT+V\,dP+\sum_i\mu_i dn_i$
- $\mu_i=(\partial G/\partial n_i)_{T,P}$ en $G=\sum_i n_i\mu_i$
- Gibbs–Duhem: $\sum_i n_i d\mu_i=-S dT + V dP$ (bij $T,P$ const: $=0$)
- $\mu_i=\mu_i^\circ+RT\ln a_i$ met $a_i$ dimensieloos
- Ideaal gas: $a_i=P_i/P^\circ$; zuivere condensed fase: $a=1$
- Fase-evenwicht: $\mu^{(\alpha)}=\mu^{(\beta)}$



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


# Hoorcollege (vermoedelijk 18 --> HC16 en mss 17) — Oplossingen, activiteiten, oplosbaarheid, Ksp, Henry

> Captions zijn automatisch gegenereerd; hieronder staat de **examenwaardige kern**.

---

## 1) Oplos-evenwicht: “suiker blijft niet eindeloos oplossen”
Voor een **moleculair** oplos-evenwicht (bv. suiker):
$$
A(s)\rightleftharpoons A(aq)
$$

Evenwichtsvoorwaarde in termen van **activiteiten**:
$$
K=\frac{a_{A(aq)}}{a_{A(s)}}.
$$

Voor een **zuivere gecondenseerde fase** (vaste stof, zuivere vloeistof) geldt:
$$
a=1 \quad \Rightarrow \quad K=a_{A(aq)}.
$$

Interpretatie: bij gegeven $T$ is er **één** activiteit (dus één concentratie) waarvoor vaste stof en opgeloste stof in evenwicht zijn → **verzadiging / oplosbaarheid**.

---

## 2) Waarom een roostermodel? (ideaal-gas-analogie voor vloeistoffen)
Doel: begrijpen wat de **activiteit** is van stoffen in oplossing, zodat je
$$
\mu_i=\mu_i^\circ+RT\ln a_i
$$
concreet kan gebruiken.

### 2.1 Roostermodel (A = oplosmiddel, B = opgeloste stof)
- Vloeibare oplossing voorgesteld als een rooster met $N$ cellen.
- Elke cel is bezet door **A** of **B**.
- Coördinatiegetal $z$ = aantal eerste naburen per cel.

### 2.2 Enthalpie: enkel **contactinteracties** (eerste buren)
Drie types contacten:
- AA met contact-enthalpie $h_{AA}$
- BB met contact-enthalpie $h_{BB}$
- AB met contact-enthalpie $h_{AB}$

Waarom enkel eerste buren?
- Van der Waals-interacties schalen typisch als $\propto 1/r^6$ → tweede buren zijn snel veel zwakker (bij $2r$ al $64\times$ kleiner).

**Fysische interpretatie bij mengen:**
- Je “breekt” een deel van AA en BB contacten
- Je “vormt” AB contacten
- Het menggedrag hangt af van hoe sterk AB is t.o.v. het gemiddelde van AA en BB.

Vaak wordt dit samengevat met een **uitwisselingsparameter** (in de les aangeduid als iets als $g$):
$$
g \sim \frac{h_{AB}-\frac{1}{2}(h_{AA}+h_{BB})}{kT}
$$
(dimensieloos: “hoe belangrijk is het bindingsverschil t.o.v. thermische energie $kT$?”)

### 2.3 Entropie van mengen (Boltzmann + combinatoriek)
Aantal configuraties:
$$
W=\frac{(N_A+N_B)!}{N_A!\,N_B!}
$$
Met Stirling:
$$
\Delta S_{\text{mix}}=-R\left(x_A\ln x_A+x_B\ln x_B\right).
$$

### 2.4 Chemische potentialen uit $G$
Vanuit $G=H-TS$:
$$
\mu_A=\left(\frac{\partial G}{\partial n_A}\right)_{T,p,n_B},\qquad
\mu_B=\left(\frac{\partial G}{\partial n_B}\right)_{T,p,n_A}.
$$

In de limieten:
- **Oplosmiddel** ($x_A\to 1$): $\mu_A \approx$ waarde van **zuivere vloeistof** ⇒ activiteit oplosmiddel vaak $\approx 1$.
- **Opgeloste stof** ($x_B\to 0$): logterm blijft belangrijk (activiteitenmodel nodig).

---

## 3) Activiteitsmodellen: molfracties vs concentraties (en waarom tabellen hier “strikt” zijn)
Altijd:
$$
\mu_i=\mu_i^\circ+RT\ln a_i
$$
Maar: je kan **kiezen** hoe je $a_i$ definieert → dan verandert ook wat je bedoelt met $\mu_i^\circ$.

### 3.1 Molfractiemodel
Voor (bijna) ideaal gedrag:
$$
a_i \approx x_i.
$$

### 3.2 Concentratiemodel (wat tabellen doorgaans gebruiken voor opgeloste soorten)
Bij sterke verdunning geldt ruwweg:
$$
x_B \approx \frac{n_B}{n_A} \propto c_B
$$
Daarom gebruikt men vaak:
$$
[B]\equiv \frac{c_B}{c^\circ},\qquad c^\circ=1~\text{mol L}^{-1}.
$$
Dan:
$$
\mu_B=\mu_B^{\circ,c}+RT\ln [B].
$$

**Cruciale examenvuistregel:**
> Als je $\mu^\circ$ of $\Delta G^\circ$ uit standaardtabellen haalt voor opgeloste ionen/stoffen, neem dan de activiteit als **onbenoemde concentratie** $[\,]$ (tenzij de opgave expliciet iets anders zegt).

---

## 4) Oplosbaarheid als evenwicht: $Q$ vs $K$
Voor het oplos-evenwicht geldt bij evenwicht:
$$
Q=K.
$$
- $Q<K$: systeem kan $G$ verlagen door **voorwaarts** te gaan → **meer oplossen** (onverzadigd)
- $Q>K$: systeem verlaagt $G$ door **achterwaarts** te gaan → **neerslag/kristallisatie** (oververzadigd)

---

## 5) Zouten in water: oplosbaarheidsproduct $K_{sp}$
Voor ionaire zouten is het oplos-evenwicht bv.:
$$
\text{AgCl}(s)\rightleftharpoons \text{Ag}^+(aq)+\text{Cl}^-(aq)
$$
Met $a(\text{vaste stof})=1$ en verdund $a\approx[\ ]$:
$$
K_{sp}=[\text{Ag}^+][\text{Cl}^-].
$$

### 5.1 Van $K_{sp}$ naar oplosbaarheid: **massabalans (vorderingsgraad)**
Definieer $s$ = molaire oplosbaarheid (mol/L).

**AgCl**
$$
[\text{Ag}^+]=s,\quad [\text{Cl}^-]=s
\Rightarrow K_{sp}=s^2 \Rightarrow s=\sqrt{K_{sp}}.
$$

**CaF$_2$**
$$
\text{CaF}_2(s)\rightleftharpoons \text{Ca}^{2+}+2\text{F}^-,
$$
$$
[\text{Ca}^{2+}]=s,\quad [\text{F}^-]=2s
\Rightarrow K_{sp}=s(2s)^2=4s^3
\Rightarrow s=\sqrt[3]{\frac{K_{sp}}{4}}.
$$

**Examenvalkuil:**
> $K_{sp}$ is een **product** van concentraties (met machten). Oplosbaarheid $s$ krijg je pas na de juiste wortel (√, ³√, …) via de stoichiometrie.

### 5.2 Gemeenschappelijk-ion-effect (oplosbaarheid is niet “één getal”)
$K_{sp}$ blijft hetzelfde, maar als je al een ion toevoegt (bv. AgNO$_3$ geeft Ag$^+$):
$$
[\text{Ag}^+]=[\text{Ag}^+]_0+s,\quad [\text{Cl}^-]=s
$$
$$
K_{sp}=([\text{Ag}^+]_0+s)\,s \approx [\text{Ag}^+]_0\,s
$$
Dus $s$ wordt **veel kleiner** → oplosbaarheid daalt sterk.

Interpretatie via $Q$:
- extra Ag$^+$ verhoogt $Q$
- als $Q>K_{sp}$ → neerslag totdat $Q=K_{sp}$.

Toepassing: zware metalen uit oplossing halen door een geschikt tegenion toe te voegen (precipitatie).

---

## 6) Temperatuursafhankelijkheid van oplosbaarheid (van ’t Hoff-idee)
Algemeen:
$$
\ln K = -\frac{\Delta G^\circ}{RT} \approx -\frac{\Delta H^\circ}{RT}+\frac{\Delta S^\circ}{R}.
$$
Dus bij (ongeveer) constante $\Delta H^\circ$:
- $\ln K$ (en vaak ook $\ln(\text{oplosbaarheid})$) is ~ lineair in $1/T$.
- Endotherm oplossen ($\Delta H^\circ>0$): oplosbaarheid ↑ bij hogere $T$.
- Exotherm oplossen ($\Delta H^\circ<0$): oplosbaarheid ↓ bij hogere $T$.

---

## 7) Opgeloste gassen: Henry
Voor een gas (bv. zuurstof):
$$
\text{O}_2(g)\rightleftharpoons \text{O}_2(aq)
$$
Activiteit gasfase:
$$
a_g=\frac{p}{p^\circ}.
$$
Voor verdund opgelost gas (vaak als concentratie):
$$
c = K_H\,p \quad \text{(Henry’s law)}.
$$

Interpretatie:
- hogere partiële druk $p$ → hogere opgeloste concentratie $c$.
- oplossen van gassen in water is vaak **exotherm** → bij hogere $T$ lost **minder** gas op.

Opmerking uit de les: databanken (bv. NIST) geven Henry-constanten soms in molaliteit (mol/kg), maar bij water en lage concentraties is dat numeriek dicht bij mol/L.

---

## 8) Examenrecept (altijd toepasbaar)
1. **Schrijf reactie** (oplossen, dissociatie, gas-oplossing).
2. **Zet activiteiten**:
   - vaste stof / zuivere vloeistof: $a=1$
   - verdunde opgeloste soorten: meestal $a\approx[i]=c_i/c^\circ$
   - gasfase: $a=p/p^\circ$
3. **Schrijf evenwichtsrelatie** ($K$, $K_{sp}$, Henry).
4. **Maak massabalans** (vorderingsgraad + beginconcentraties).
5. **Los op** (let op machten → juiste wortel).
6. **Interpreteer met $Q$ vs $K$** (richting, neerslag, over-/onderverzadigd).


# HC19–20 — Oplosbaarheid (Kₛ) & gemeenschappelijk-ion-effect (examengericht)

## 1) Wat de prof **echt** wil dat je snapt

- **Oplosbaarheid is een evenwichtsprobleem.** Een slecht oplosbaar zout staat in evenwicht tussen vaste stof en ionen in oplossing.
- **Het oplosbaarheidsproduct $K_s$ is een evenwichtsconstante** (dus thermodynamisch gelinkt aan $\Delta_r G^\circ$):

  $$
  K_s = \exp\!\left(-\frac{\Delta_r G^\circ}{RT}\right)
  $$

- **$K_s$ ≠ “oplosbaarheid”.** $K_s$ zegt iets over het evenwicht, maar de *molaire oplosbaarheid* $S$ (mol/L) hangt ook af van:
  - **stoichiometrie** van het zout  
  - **samenstelling** van de oplossing (gemeenschappelijke ionen, complexvorming, pH, …)  
  - **temperatuur**

**Cruciale examenvuistregel:**  
> Je mag $K_s$ van verschillende zouten **niet** zomaar onderling vergelijken om te zeggen welk zout het “meest oplosbaar” is, tenzij je eerst de **molaire oplosbaarheid $S$** uitrekent.

---

## 2) Basisopstelling: hoe je altijd start

Neem een algemeen zout:

$$
M_aX_b(s) \rightleftharpoons a\,M^{n+}(aq) + b\,X^{m-}(aq)
$$

Dan (in verdunde oplossingen, met activiteiten $\approx$ concentraties):

$$
K_s = [M^{n+}]^a [X^{m-}]^b
$$

### Molaire oplosbaarheid $S$ in zuiver water

Laat $S = x$ mol/L zout oplossen. Dan:

- $[M^{n+}] = a x$  
- $[X^{m-}] = b x$

Dus:

$$
K_s = (a x)^a (b x)^b
\quad\Rightarrow\quad
x = S
$$

**Mini-cheatsheet (komt rechtstreeks uit de slides/oefeningen):**

- $MX(s) \rightleftharpoons M^+ + X^-$  
  $K_s = x^2 \Rightarrow S = \sqrt{K_s}$

- $M_2X(s) \rightleftharpoons 2M^+ + X^{2-}$  
  $K_s = (2x)^2 x = 4x^3 \Rightarrow S = \left(\frac{K_s}{4}\right)^{1/3}$

- $M_2X_3(s) \rightleftharpoons 2M^{3+} + 3X^{2-}$  
  $K_s = (2x)^2 (3x)^3 = 108x^5 \Rightarrow S = \left(\frac{K_s}{108}\right)^{1/5}$

**Waarom dit belangrijk is:** in de slides staat expliciet een oefening waar zouten met verschillende stoichiometrie een $K_s$ krijgen en je *moet* rangschikken op oplosbaarheid. Dat kan enkel via $S$.

---

## 3) Gemeenschappelijk-ion-effect (dé klassieker)

Als er al een ion aanwezig is dat in het evenwicht voorkomt, dan verschuift het evenwicht **naar links** ⇒ **oplosbaarheid daalt**.

### Voorbeeldtype uit de slides: $Ag_2CrO_4$

$$
Ag_2CrO_4(s) \rightleftharpoons 2Ag^+ + CrO_4^{2-}
$$

$$
K_s = [Ag^+]^2 [CrO_4^{2-}]
$$

Stel: je hebt al $0{,}10\ \text{mol/L}$ $Ag^+$ in oplossing, en extra oplosbaarheid is $x$:

- $[Ag^+] = 0{,}10 + 2x$  
- $[CrO_4^{2-}] = x$

Dan:

$$
K_s = (0{,}10 + 2x)^2 x
$$

**De exametruc die de prof toont:**  
Als $0{,}10 \gg 2x$, dan:

$$
(0{,}10 + 2x) \approx 0{,}10
\quad\Rightarrow\quad
K_s \approx (0{,}10)^2 x
\Rightarrow x \approx \frac{K_s}{(0{,}10)^2}
$$

**Altijd doen op examen:** check achteraf of $2x \ll 0{,}10$ echt klopt.  
Als dat niet klopt → geen benadering, dan los je exact op (vaak een kubiek/kwadratisch, maar meestal is de benadering juist gekozen).

---

## 4) Stappenplan dat bijna elke oefening oplost

1. **Schrijf de oplosreactie** correct (met stoichiometrie).  
2. **Schrijf $K_s$** als product van ionconcentraties met machten.  
3. **ICE/massabalans**: zet beginconcentraties + verandering door oplossen ($x$).  
4. **Substitueer** alles in $K_s$.  
5. **Kies (indien mogelijk) een benadering** (bv. gemeenschappelijk ion domineert).  
6. **Los op voor $x$ (= $S$)**.  
7. **Consistentiecheck** van je benadering.  

---

## 5) Typische valkuilen (waar punten verdampen)

- $K_s$-waarden vergelijken zonder naar stoichiometrie te kijken.  
- Vergeten dat bij $M_2X$: $[M] = 2S$ en niet $S$.  
- Benadering maken (“$0{,}10 + 2x \approx 0{,}10$”) en **niet** checken.  
- Denken dat “groter $K_s$ altijd groter $S$” is — *soms*, maar niet universeel.

---

## 6) Wat je “paraat” wil hebben voor het examen

- Definitie: $K_s$ als evenwichtsconstante van oplossen.  
- Omzetting $K_s \leftrightarrow S$ via stoichiometrie (voor 1:1, 2:1, 2:3 moet je dit snel kunnen).  
- Gemeenschappelijk-ion-effect kunnen opstellen én benaderen.  
- Het mantra: **“$K_s$ is niet oplosbaarheid; $S$ is oplosbaarheid.”**  

---

# HC21 – Elektrochemie (examengerichte samenvatting)

## 1) De kernidee (wat je *altijd* moet kunnen)

Een **galvanische (volta-)cel** zet een **spontane redoxreactie** om in **elektrische arbeid** via elektronen die door een uitwendig circuit lopen.

Spontaniteit-criteria (cruciaal):

- **Spontaan:** $E_\text{cel} > 0 \;\Leftrightarrow\; \Delta G < 0$  
- **Evenwicht / cel “plat”:** $E_\text{cel} = 0 \;\Leftrightarrow\; \Delta G = 0$ en dan geldt $Q = K$

**Anode/Kathode (altijd examenvalkuil):**

- **Anode = oxidatie**  
- **Kathode = reductie**  
- In een **galvanische** cel: anode is typisch **negatief**, kathode **positief**.

---

## 2) Opbouw van een galvanische cel (Daniell/Zn–Cu is het archetype)

- 2 **halfcellen** (elk een elektrode + elektrolyt).  
- **Zoutbrug** (bv. KNO₃ of NH₄NO₃ in gel) zorgt voor ionenmigratie zodat elke halfcel elektrisch neutraal blijft.  
- **Uitwendig circuit**: elektronen lopen van anode → kathode.

**Actieve vs inerte elektroden**

- Actief: de elektrode doet mee in de redox (bv. Zn(s), Cu(s)).  
- Inert: geleidt enkel elektronen, doet niet mee (bv. Pt(s), grafiet).

---

## 3) Elektrodepotentialen en celspanning (rekenrecept)

**Standaardreductiepotentialen $E^\circ$** zijn gedefinieerd t.o.v. de **standaard-waterstofelektrode (SHE)** met $E^\circ = 0{,}000\ \text{V}$.

Belangrijkste formule:

$$
E^\circ_\text{cel} = E^\circ_{\text{kathode \,(reductie)}} - E^\circ_{\text{anode \,(reductie)}}
$$

**Examenschema om $E^\circ_\text{cel}$ te vinden**

1. Schrijf beide halfreacties als **reducties** (zoals in de tabel met $E^\circ$).  
2. De halfreactie met **grootste $E^\circ$** wordt de **kathode** (gaat effectief als reductie).  
3. De andere wordt **anode** (gaat effectief als oxidatie; teken keert om in je reactie, maar je gebruikt in de formule nog steeds het *reductie*-$E^\circ$).  
4. Bereken $E^\circ_\text{cel}$ met de formule hierboven.  
5. Balanceer elektronen voor de totale reactie, maar let op:
   - **Je vermenigvuldigt $E^\circ$ nooit met stoichiometrische factoren.**

---

## 4) Link met thermodynamica: arbeid, vrije energie en evenwicht

Elektrische arbeid en vrije energie hangen samen via:

$$
\Delta G = -n F E_\text{cel}
$$

waar $n$ = aantal uitgewisselde elektronen en $F$ = Faraday-constante.

Standaardrelaties (super-examenwaardig):

$$
\Delta G^\circ = -n F E^\circ_\text{cel}
$$

$$
\Delta G^\circ = -R T \ln K
$$

$$
E^\circ_\text{cel} = \frac{R T}{n F} \ln K
$$

Interpretatie in één oogopslag:

- $\Delta G^\circ < 0 \Rightarrow K > 1 \Rightarrow E^\circ_\text{cel} > 0$ (spontaan in standaardcondities)  
- $\Delta G^\circ = 0 \Rightarrow K = 1 \Rightarrow E^\circ_\text{cel} = 0$ (evenwicht)  
- $\Delta G^\circ > 0 \Rightarrow K < 1 \Rightarrow E^\circ_\text{cel} < 0$ (niet spontaan)

---

## 5) Nernstvergelijking: niet-standaard condities

Als concentraties/drukken afwijken van standaardcondities, gebruik je **Nernst**:

$$
E = E^\circ - \frac{R T}{n F} \ln Q
$$

Bij $25^\circ\text{C}$ vaak als:

$$
E = E^\circ - \frac{0{,}05916}{n} \log_{10} Q
$$

**Wat is $Q$?**  
Het reactiequotiënt: producten / reactanten met macht volgens stoichiometrie.

- Zuivere vaste stoffen en vloeistoffen komen **niet** in $Q$.

---

## 6) Concentratiecellen (klassiek examenvraagstuk)

Een **concentratiecel** heeft dezelfde halfreactie links en rechts, maar met **verschillende concentraties**.

Voorbeeldnotatie:

$$
Ag(s) \;|\; Ag^+(1{,}0\,\text{M}) \;||\; Ag^+(0{,}10\,\text{M}) \;|\; Ag(s)
$$

Belangrijk:

- De celspanning is typisch **klein**.  
- Je berekent elke halfcelpotentiaal met **Nernst**, en dan:

  $$
  E_\text{cel} = E_\text{rechts} - E_\text{links}
  $$

  (of consistent met kathode–anode)  

- De kant met **hogere reductiepotentiaal** fungeert als **kathode**.

---

## 7) Commerciële voltacellen (wat je moet herkennen + kernreacties)

Je hoeft dit meestal niet hyper-diep af te leiden, maar je moet het **type, de idee en vaak de halfreacties kunnen plaatsen**.

### Loodaccumulator (auto, ~12 V totaal)

- Ongeveer **2 V per cel**, typisch 6 cellen in serie.  
- Bij ontlading wordt **$H_2SO_4$ verbruikt** (dichtheid daalt).  
- Oplaadbaar: reacties omkeerbaar via externe stroombron.  
- Bij laden kan water-elektrolyse ($H_2/O_2$) optreden → veiligheidsaspect.

### Droge cellen (1,25–1,50 V)

- **Leclanché-element** (klassieke “zink-kool”-achtige cel).  
- **Alkalische batterij**: langere levensduur (zinkanode corrodeert trager in basisch milieu).  
- **Zilvercel**: gebruikt in kleine toestellen (uurwerken, pacemakers, hoorapparaten, …).

### Nikkel–cadmium (Ni–Cd, ~1,4 V)

- Oplaadbaar (producten blijven aan elektroden “kleven” volgens de cursuscontext).  
- Toepassingen: boormachines, scheerapparaten, …

### Brandstofcel (H$_2$/O$_2$, $E_\text{cel} \approx 1{,}2\ \text{V}$)

- Reagentia worden continu aangevoerd.  
- Nettoreactie:  
  $$
  2H_2 + O_2 \rightarrow 2H_2O
  $$
- Efficiëntie-idee: groot deel van theoretische $\Delta G$ → elektrische energie.  
- Nadelen: opslag van reagentia en dure elektroden.

### Li-ion (conceptueel herkennen)

- **Anode: grafiet**  
- **Kathode: LiCoO_2**  
- **Li⁺-transport** tussen elektroden (intercalatie/de-intercalatie als idee).

---

## 8) Examenvallen & mini-checklist

- **Anode ≠ altijd positief**: in een galvanische cel is de anode typisch **negatief** (oxidatie), kathode **positief** (reductie).  
- **$E^\circ$ nooit schalen met coëfficiënten.**  
- **Tabelwaarden zijn reductiepotentialen**: als je een oxidatie gebruikt, keer je de reactie om maar je werkt nog steeds met reductie-$E^\circ$ in  
  $$
  E^\circ_\text{cel} = E^\circ_\text{kath} - E^\circ_\text{an}.
  $$
- **Standaardcondities**: opgeloste species 1 M, gassen 1 bar, zuivere vaste stoffen/vloeistoffen activiteit 1.  
- Bij “cel plat”: **$E_\text{cel} = 0$ én $\Delta G = 0$ én $Q = K$.**

---

## 9) Snelle “rekenflow” (wat je op je kladpapier wil)

1. Identificeer halfreacties + $E^\circ$ uit tabel.  
2. Kies kathode = hoogste $E^\circ$.  
3. Bereken $E^\circ_\text{cel} = E^\circ_\text{kath} - E^\circ_\text{an}$.  
4. Balanceer de volledige reactie → bepaal $n$.  
5. Indien niet-standaard: gebruik  
   $$
   E = E^\circ - \frac{R T}{n F} \ln Q.
   $$  
6. Gebruik $\Delta G = -n F E$ en eventueel  
   $$
   E^\circ_\text{cel} = \frac{R T}{n F} \ln K
   $$  
   om $K$ te vinden.
