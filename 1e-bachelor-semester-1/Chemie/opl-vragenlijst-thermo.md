# Thermodynamica — Oplossingen (afleidingen)

Hieronder staan **model-afleidingen** (de “standaard ketting” die je op examen wil kunnen reproduceren).  
Notatie: $d$ voor exacte differentiaal (toestandsfunctie), $\delta$ voor pad-afhankelijke procesgrootheden.

---

## 1) Eerste hoofdwet + arbeid

### 1.1 Waarom $dU$ maar $\delta q,\delta w$?
- $U$ is een **toestandsfunctie**: $U=U(\text{toestand})$. Dus de verandering tussen twee toestanden hangt enkel af van begin/eind, niet van het pad. Daarom is $dU$ een **exacte differentiaal**.
- Warmte $q$ en arbeid $w$ zijn **energieoverdracht langs een procespad**. Hun waarden hangen af van het pad → **niet-exact** → $\delta q,\delta w$.

De eerste hoofdwet (gesloten systeem):
$$
dU=\delta q+\delta w
$$
waar het teken van $w$ afhangt van conventie. In jouw lijst gebruiken we de conventie waarin PV-arbeid negatief is bij expansie (werk door het systeem).

### 1.2 Splitsing van arbeid en PV-arbeid
Algemeen kan je arbeid splitsen:
$$
\delta w=\delta w_{PV}+\delta w_{\text{nuttig}}
$$
PV-arbeid: bij een zuiger met uitwendige druk $P_{\text{omg}}$ en volumeverandering $dV$.
- Kracht: $F=P_{\text{omg}}A$.
- Verplaatsing: $dx$, met $dV=A\,dx$.
- Werk **op het systeem**: $\delta w = F\,dx = P_{\text{omg}}A\,dx=P_{\text{omg}}\,dV$.
- Werk **door het systeem** (jouw conventie): $\delta w_{PV}=-P_{\text{omg}}\,dV$.

Dus:
$$
\delta w_{PV}=-P_{\text{omg}}\,dV.
$$

### 1.3 Combineer tot totale energiebalans
Invullen in de eerste hoofdwet:
$$
dU=\delta q+\delta w=\delta q - P_{\text{omg}}\,dV + \delta w_{\text{nuttig}}.
$$

**Impliciete assumpties**: gesloten systeem; macroscopische PV-arbeid via een goed gedefinieerde uitwendige druk; “nuttige arbeid” bundelt alle niet-PV bijdragen (elektrisch, oppervlaktespanning, aswerk, …).

### 1.4 Reversibel vs irreversibel (maximale arbeid)
- Reversibel betekent “quasi-statisch” en zonder dissipatie.
- Dissipatie (wrijving, turbulentie, eindige $\Delta T$ bij warmteoverdracht) produceert entropie en “verspilt” potentiële nuttige arbeid.
- Daarom levert een reversibel pad, bij gegeven begin/eindtoestand, **maximale** nuttige arbeid.

---

## 2) Warmtecapaciteiten $C_V$ en $C_P$

### 2.1 Afleiding bij constant volume
Neem gesloten systeem, geen nuttige arbeid: $\delta w_{\text{nuttig}}=0$, en $dV=0$:
$$
dU=\delta q - P_{\text{omg}}\,dV = \delta q.
$$
Definieer warmtecapaciteit bij constant volume:
$$
C_V \equiv \left(\frac{\delta q}{dT}\right)_V.
$$
Dus:
$$
dU = C_V\,dT \quad \Rightarrow \quad C_V=\left(\frac{\partial U}{\partial T}\right)_V.
$$

### 2.2 Afleiding bij constante druk
Bij $P=\text{const}$ en geen nuttige arbeid:
$$
dU=\delta q - P_{\text{omg}}\,dV.
$$
Herschik:
$$
\delta q = dU + P_{\text{omg}}\,dV.
$$
Definieer $C_P$:
$$
C_P \equiv \left(\frac{\delta q}{dT}\right)_P
\quad \Rightarrow \quad
\delta q = C_P\,dT.
$$
De extra term $P\,dV$ is de energie die nodig is om **ruimte te maken tegen de omgeving** tijdens expansie.

---

## 3) Enthalpie $H$ en differentiaal

### 3.1 Differentiaal van enthalpie
Definitie:
$$
H = U + P_{\text{omg}}V.
$$
Differentieer:
$$
dH = dU + P_{\text{omg}}\,dV + V\,dP_{\text{omg}}.
$$
Als $P_{\text{omg}}=\text{const}$, dan $dP_{\text{omg}}=0$ en:
$$
dH = dU + P_{\text{omg}}\,dV.
$$

### 3.2 Link met $C_P$
Uit sectie 2.2: bij $P=\text{const}$ en $\delta w_{\text{nuttig}}=0$ geldt:
$$
\delta q = dU + P_{\text{omg}}\,dV = dH.
$$
Dus:
$$
dH = C_P\,dT \quad \Rightarrow \quad C_P=\left(\frac{\partial H}{\partial T}\right)_P.
$$

---

## 4) Hess + standaardenthalpieën

### 4.1 Wet van Hess
Omdat $H$ een toestandsfunctie is:
$$
\Delta H_{\text{tot}} = H_{\text{eind}}-H_{\text{begin}}
$$
is **pad-onafhankelijk**. Elke reactie kan je opdelen in een som van deelreacties; de totale $\Delta H$ is de som van de deel-$\Delta H$’s:
$$
\Delta H = \sum_k \Delta H_k.
$$
Dat is Hess.

### 4.2 $\Delta_r H^\circ$ uit vormingsenthalpieën
Schrijf een reactie met stoichiometrische coëfficiënten $\nu_i$ (producten positief, reagentia negatief):
$$
\sum_i \nu_i A_i = 0.
$$
Dan:
$$
\Delta_r H^\circ = \sum_i \nu_i H_{m,i}^\circ.
$$
Als je expliciet “producten minus reagentia” wil:
$$
\Delta_r H^\circ=\sum_P \nu_P H_{m,P}^\circ-\sum_R \nu_R H_{m,R}^\circ.
$$

---

## 5) Entropie + tweede hoofdwet

### 5.1 $dS=\delta q_{\text{rev}}/T$ en ongelijkheid
Clausius-definitie (reversibel):
$$
dS \equiv \frac{\delta q_{\text{rev}}}{T}.
$$
Voor irreversibele processen geldt Clausius-ongelijkheid:
$$
dS \ge \frac{\delta q}{T},
$$
met strikt $>$ voor echte irreversibiliteit (entropieproductie).

### 5.2 Bij constant $P,T$ en geen nuttige arbeid: $dS=dH/T$
Bij reversibel en $P=\text{const}$:
$$
\delta q_{\text{rev}} = dH.
$$
Dus:
$$
dS = \frac{\delta q_{\text{rev}}}{T}=\frac{dH}{T}.
$$

---

## 6) Gibbs vrije energie $G$: definitie, spontaniteit, differentiaal

### 6.1 Differentiaal bij constant $T$
Definitie:
$$
G=H-TS.
$$
Differentieer:
$$
dG=dH - T\,dS - S\,dT.
$$
Bij $T=\text{const}$:
$$
dG = dH - T\,dS.
$$

### 6.2 Spontaniteitscriterium bij $P,T$ constant
Neem gesloten systeem. Schrijf de eerste en tweede hoofdwet samen.

Eerste hoofdwet:
$$
dU=\delta q - P_{\text{omg}}\,dV + \delta w_{\text{nuttig}}.
$$
Definieer $H=U+P_{\text{omg}}V$ (met $P_{\text{omg}}=\text{const}$):
$$
dH = dU + P_{\text{omg}}\,dV = \delta q + \delta w_{\text{nuttig}}.
$$
Tweede hoofdwet (algemeen):
$$
dS \ge \frac{\delta q}{T}
\quad \Rightarrow \quad
\delta q \le T\,dS.
$$
Combineer:
$$
dH = \delta q + \delta w_{\text{nuttig}} \le T\,dS + \delta w_{\text{nuttig}}.
$$
Herschik:
$$
dH - T\,dS \le \delta w_{\text{nuttig}}.
$$
Maar bij $P,T=\text{const}$ is $dG=dH-T\,dS$, dus:
$$
dG \le \delta w_{\text{nuttig}}.
$$
Als er **geen** nuttige arbeid mogelijk is ($\delta w_{\text{nuttig}}=0$), dan:
$$
dG \le 0,
$$
met $dG=0$ in evenwicht.

### 6.3 Fundamentele identiteit $dG=V\,dP-S\,dT$
Start van:
$$
G=U+PV-TS.
$$
Differentieer:
$$
dG=dU + P\,dV + V\,dP - T\,dS - S\,dT.
$$
Gebruik de fundamentele relatie voor een simpel comprimeerbaar systeem:
$$
dU = T\,dS - P\,dV.
$$
Invullen:
$$
dG = (T\,dS - P\,dV) + P\,dV + V\,dP - T\,dS - S\,dT
= V\,dP - S\,dT.
$$

### 6.4 Partiële afgeleide
Uit $dG = V\,dP - S\,dT$ volgt bij constant $P$:
$$
\left(\frac{\partial G}{\partial T}\right)_P = -S.
$$

---

## 7) Temperatuurafhankelijkheid via $G/T$

We willen:
$$
\left(\frac{\partial (G/T)}{\partial T}\right)_P=-\frac{H}{T^2}.
$$

Start:
$$
\frac{G}{T} = \frac{H-TS}{T} = \frac{H}{T} - S.
$$
Differentieer bij constant $P$:
$$
\left(\frac{\partial (G/T)}{\partial T}\right)_P
=
\left(\frac{\partial (H/T)}{\partial T}\right)_P
-
\left(\frac{\partial S}{\partial T}\right)_P.
$$
Werk de eerste term uit:
$$
\left(\frac{\partial (H/T)}{\partial T}\right)_P
=
\frac{1}{T}\left(\frac{\partial H}{\partial T}\right)_P
-\frac{H}{T^2}.
$$
Maar $\left(\frac{\partial H}{\partial T}\right)_P = T\left(\frac{\partial S}{\partial T}\right)_P$ (want bij $P=\text{const}$: $dS=dH/T$).
Dus:
$$
\frac{1}{T}\left(\frac{\partial H}{\partial T}\right)_P
=
\left(\frac{\partial S}{\partial T}\right)_P.
$$
Die twee termen cancellen, en je houdt over:
$$
\left(\frac{\partial (G/T)}{\partial T}\right)_P=-\frac{H}{T^2}.
$$

---

## 8) Drukafhankelijkheid van $G$

### 8.1 Algemeen: $(\partial G/\partial P)_T=V$
Uit $dG=V\,dP-S\,dT$:
bij constant $T$ is $dT=0$:
$$
dG = V\,dP
\quad \Rightarrow \quad
\left(\frac{\partial G}{\partial P}\right)_T = V.
$$

### 8.2 Ideaal gas: log-relatie
Voor een ideaal gas: $V=\frac{nRT}{P}$.
Dus:
$$
\left(\frac{\partial G}{\partial P}\right)_T=\frac{nRT}{P}.
$$
Integreer van $P^\circ$ naar $P$ bij constant $T$:
$$
G(P)-G(P^\circ) = \int_{P^\circ}^{P}\frac{nRT}{P'}\,dP'
= nRT\ln\!\left(\frac{P}{P^\circ}\right).
$$
Dus:
$$
G(P)=G(P^\circ)+nRT\ln\!\left(\frac{P}{P^\circ}\right).
$$

---

## 9) Helmholtz-energie $A$

### 9.1 Differentiaal $dA$
Definitie:
$$
A=U-TS.
$$
Differentieer:
$$
dA=dU - T\,dS - S\,dT.
$$
Gebruik $dU=T\,dS-P\,dV$:
$$
dA = (T\,dS-P\,dV) - T\,dS - S\,dT
= -P\,dV - S\,dT.
$$
Dus:
$$
dA=-S\,dT - P\,dV.
$$

### 9.2 Partiële afgeleiden
Uit $dA=-S\,dT-P\,dV$:
$$
\left(\frac{\partial A}{\partial T}\right)_V=-S,
\qquad
\left(\frac{\partial A}{\partial V}\right)_T=-P.
$$

---

## 10) Chemische potentiaal + Gibbs–Duhem

### 10.1 Definitie chemische potentiaal
Voor een multicomponent systeem:
$$
G=G(T,P,n_1,\dots,n_k).
$$
Differentiaal:
$$
dG=\left(\frac{\partial G}{\partial T}\right)_{P,n}dT
+\left(\frac{\partial G}{\partial P}\right)_{T,n}dP
+\sum_i \left(\frac{\partial G}{\partial n_i}\right)_{T,P,n_{j\ne i}}dn_i.
$$
Met de definitie:
$$
\mu_i \equiv \left(\frac{\partial G}{\partial n_i}\right)_{T,P,n_{j\ne i}}.
$$

### 10.2 Toon $G=\sum_i n_i\mu_i$ (extensiviteit)
Voor vaste $T,P$ is $G$ **extensief**: schaal alle $n_i \to \lambda n_i$, dan $G\to \lambda G$.
Dus $G$ is homogeen van graad 1 in $(n_i)$. Euler’s stelling:
$$
G=\sum_i n_i\left(\frac{\partial G}{\partial n_i}\right)_{T,P,n_{j\ne i}}
=\sum_i n_i\mu_i.
$$

### 10.3 Gibbs–Duhem
Start:
$$
G=\sum_i n_i\mu_i.
$$
Differentieer:
$$
dG = \sum_i \mu_i\,dn_i + \sum_i n_i\,d\mu_i.
$$
Maar uit de algemene vorm:
$$
dG=-S\,dT + V\,dP + \sum_i \mu_i\,dn_i.
$$
Vergelijk beide uitdrukkingen → de $\sum_i \mu_i dn_i$ vallen weg:
$$
\sum_i n_i\,d\mu_i = -S\,dT + V\,dP.
$$
Herschik:
$$
S\,dT - V\,dP + \sum_i n_i\,d\mu_i = 0.
$$

### 10.4 Speciale gevallen
- **Zuivere stof**: $G=n\mu$, dus molair:
  $$
  \mu = G_m.
  $$
- **Ideaal gas**: later (sectie 8/13) volgt $\mu=\mu^\circ + RT\ln(P/P^\circ)$.

---

## 11) Fase-evenwicht: Clausius–Clapeyron + fasenregel

### 11.1 Geïntegreerde Clausius–Clapeyron
Clapeyron:
$$
\frac{dP}{dT}=\frac{\Delta H_{\alpha\to\beta}}{T\,\Delta V}.
$$
Voor verdamping/sublimatie: $\Delta V \approx V_{\text{gas}}$ en ideaal gas:
$$
V_{\text{gas}}=\frac{RT}{P} \quad \Rightarrow \quad \Delta V \approx \frac{RT}{P}.
$$
Dan:
$$
\frac{dP}{dT}=\frac{\Delta H}{T}\frac{P}{RT}
=\frac{\Delta H}{R}\frac{P}{T^2}.
$$
Dus:
$$
\frac{1}{P}\,dP = \frac{\Delta H}{R}\frac{1}{T^2}\,dT.
$$
Integreer tussen $(P_1,T_1)$ en $(P_2,T_2)$ (met $\Delta H$ constant benaderd):
$$
\int_{P_1}^{P_2}\frac{1}{P}\,dP
=
\frac{\Delta H}{R}\int_{T_1}^{T_2}\frac{1}{T^2}\,dT.
$$
Links: $\ln(P_2/P_1)$. Rechts: $\frac{\Delta H}{R}\left[-\frac{1}{T}\right]_{T_1}^{T_2}
=\frac{\Delta H}{R}\left(-\frac{1}{T_2}+\frac{1}{T_1}\right)$.
Dus:
$$
\ln\!\left(\frac{P_2}{P_1}\right)
=
-\frac{\Delta H}{R}\left(\frac{1}{T_2}-\frac{1}{T_1}\right).
$$
Exponentieer:
$$
P_2 = P_1\exp\!\left[-\frac{\Delta H}{R}\left(\frac{1}{T_2}-\frac{1}{T_1}\right)\right].
$$

### 11.2 Fasenregel
Aantal vrijheidsgraden $f$ = aantal intensieve variabelen die je onafhankelijk kan kiezen zonder het aantal fasen te veranderen.
Gibbs fasenregel:
$$
f = 2 + c - p,
$$
met $c$ componenten en $p$ fasen.
- “2” komt van de typische intensieve variabelen $T$ en $P$.
- Elke extra component geeft een extra samenstellingsvrijheid, elke extra fase legt een evenwichtsvoorwaarde op.

---

## 12) Chemisch evenwicht: $\xi$ en $\Delta_r G$

Start met:
$$
dG=\sum_i \mu_i\,dn_i.
$$
Voor reactie $\sum_i \nu_i A_i =0$ (producten $\nu_i>0$, reagentia $\nu_i<0$).
Definieer vorderingsgraad $\xi$:
$$
dn_i = \nu_i\,d\xi.
$$
Invullen:
$$
dG=\sum_i \mu_i \nu_i\,d\xi
=\left(\sum_i \nu_i\mu_i\right)d\xi.
$$
Definieer:
$$
\Delta_r G \equiv \sum_i \nu_i\mu_i
\quad \Rightarrow \quad
dG=\Delta_r G\,d\xi.
$$
Evenwicht bij vast $T,P$ is minimum van $G$ → $dG=0$ bij variaties in $\xi$:
$$
\Delta_r G = 0 \quad (\text{evenwicht}).
$$
Teken:
- $\Delta_r G<0$: vooruitgaande reactie verlaagt $G$ → spontaan vooruit.
- $\Delta_r G>0$: spontaan achteruit.
- $\Delta_r G=0$: evenwicht.

---

## 13) Activiteiten, $Q$ en $K$

### 13.1 Van $\mu=\mu^\circ+RT\ln a$ naar $\Delta G=\Delta G^\circ+RT\ln Q$
Neem:
$$
\mu_i=\mu_i^\circ+RT\ln a_i.
$$
Dan:
$$
\Delta_r G = \sum_i \nu_i\mu_i
= \sum_i \nu_i\mu_i^\circ + RT\sum_i \nu_i \ln a_i.
$$
Definieer:
$$
\Delta_r G^\circ \equiv \sum_i \nu_i\mu_i^\circ.
$$
En gebruik $\sum_i \nu_i \ln a_i = \ln\left(\prod_i a_i^{\nu_i}\right)$:
$$
\Delta_r G = \Delta_r G^\circ + RT\ln\!\left(\prod_i a_i^{\nu_i}\right).
$$
Definieer reactiequotiënt:
$$
Q_a \equiv \prod_i a_i^{\nu_i}
= \frac{\prod_P a_P^{\nu_P}}{\prod_R a_R^{\nu_R}}.
$$
Dus:
$$
\Delta_r G=\Delta_r G^\circ + RT\ln Q_a.
$$

### 13.2 Evenwicht: $K=\exp(-\Delta_r G^\circ/RT)$
Bij evenwicht: $\Delta_r G=0$ en $Q_a=K_a^\circ$:
$$
0=\Delta_r G^\circ + RT\ln K_a^\circ
\quad \Rightarrow \quad
\ln K_a^\circ = -\frac{\Delta_r G^\circ}{RT}.
$$
Dus:
$$
K_a^\circ=\exp\!\left(-\frac{\Delta_r G^\circ}{RT}\right).
$$

---

## 14) Ideale gassen: $Q_P$, $K_P^\circ$ en $K_\chi$

### 14.1 Activiteit als $P/P^\circ$
Voor ideale gassen is de activiteit (standaardtoestand: $P^\circ$):
$$
a_i=\frac{P_i}{P^\circ}.
$$
Dus:
$$
Q_P=\prod_i \left(\frac{P_i}{P^\circ}\right)^{\nu_i}.
$$
Bij evenwicht: $Q_P=K_P^\circ$.

### 14.2 Naar molfracties
Met $P_i=\chi_i P_{\text{tot}}$:
$$
Q_P = \prod_i \left(\frac{\chi_i P_{\text{tot}}}{P^\circ}\right)^{\nu_i}
= \left(\frac{P_{\text{tot}}}{P^\circ}\right)^{\sum_i \nu_i}\prod_i \chi_i^{\nu_i}.
$$
Definieer:
$$
\Delta\nu \equiv \sum_i \nu_i
\quad (\text{som producten minus som reagentia}).
$$
En:
$$
Q_\chi \equiv \prod_i \chi_i^{\nu_i}
= \frac{\prod_P \chi_P^{\nu_P}}{\prod_R \chi_R^{\nu_R}}.
$$
Dan:
$$
Q_P = \left(\frac{P_{\text{tot}}}{P^\circ}\right)^{\Delta\nu} Q_\chi.
$$
Bij evenwicht: $Q_P=K_P^\circ$ en $Q_\chi=K_\chi$:
$$
K_P^\circ=\left(\frac{P_{\text{tot}}}{P^\circ}\right)^{\Delta\nu}K_\chi
\quad \Rightarrow \quad
K_\chi = K_P^\circ\left(\frac{P^\circ}{P_{\text{tot}}}\right)^{\Delta\nu}.
$$

---

## 15) Massabalans met $\alpha$ (gasreacties)

Neem reactie $\sum_i \nu_i A_i=0$ en beginhoeveelheden $n_{i,0}$.
Met extent $\xi$:
$$
n_i = n_{i,0}+\nu_i\xi.
$$
Soms gebruikt men een dimensieloze omzettingsgraad $\alpha$ met $\xi = \alpha n_0$ (waar $n_0$ een gekozen referentiehoeveelheid is, vaak de beginhoeveelheid van limiterend reagens).
Dan:
$$
n_i = n_{i,0}+\nu_i \alpha n_0.
$$
Voor reagentia (negatieve $\nu$) schrijft dit als “min”:
$$
n_R = n_{R,0} - |\nu_R|\alpha n_0.
$$
Molfracties:
$$
\chi_i(\alpha)=\frac{n_i(\alpha)}{\sum_j n_j(\alpha)}.
$$
En dan bouw je:
$$
Q_\chi(\alpha)=\prod_i \chi_i(\alpha)^{\nu_i}.
$$

---

## 16) Oplossingen: roostermodel → $S$, $H$, $G$, $\mu$

### 16.1 Mengentropie (ideale oplossing)
Neem $N$ roosterplaatsen, met $N_A$ deeltjes A en $N_B$ deeltjes B, $N=N_A+N_B$.
Aantal configuraties:
$$
\Omega=\frac{N!}{N_A!\,N_B!}.
$$
Boltzmann:
$$
S=k_B\ln\Omega.
$$
Gebruik Stirling: $\ln N!\approx N\ln N - N$:
$$
\ln\Omega \approx N\ln N - N_A\ln N_A - N_B\ln N_B.
$$
Schrijf $N_A=x_A N$, $N_B=x_B N$ met $x_A+x_B=1$:
$$
\ln\Omega \approx -N\left(x_A\ln x_A + x_B\ln x_B\right).
$$
Dus:
$$
\Delta S_{\text{mix}} = -Nk_B\left(x_A\ln x_A + x_B\ln x_B\right).
$$
Per mol (met $R=N_A^\text{(Av)}k_B$):
$$
\Delta S_{m,\text{mix}} = -R\left(\chi_A\ln\chi_A+\chi_B\ln\chi_B\right).
$$

### 16.2 Mengenthalpie (regelmatige oplossing – typische vorm)
In het simpelste “regular solution”-model:
$$
\Delta H_{\text{mix}} = \Omega\,n\,\chi_A\chi_B
$$
waar $\Omega$ een interactieparameter is (verschil in A–B bindingen t.o.v. A–A en B–B). De kern is: nul voor $\chi_A=0$ of $1$, maximum rond $\chi_A=\chi_B=1/2$.

### 16.3 Gibbs bij mengen
$$
\Delta G_{\text{mix}}=\Delta H_{\text{mix}}-T\Delta S_{\text{mix}}
= \Omega n\chi_A\chi_B + RT\,n\left(\chi_A\ln\chi_A+\chi_B\ln\chi_B\right).
$$
(Het entropiedeel is negatief, dus $-T\Delta S$ wordt positief maal log-termen.)

### 16.4 Chemische potentialen
Voor mengsels:
$$
\mu_i=\left(\frac{\partial G}{\partial n_i}\right)_{T,P,n_{j\ne i}}.
$$
Voor een ideale oplossing (waar $\Delta H_{\text{mix}}=0$) krijg je:
$$
\mu_i=\mu_i^\circ + RT\ln \chi_i,
$$
en in het regular-solution model komen er extra termen bij uit $\Delta H_{\text{mix}}$ (typisch lineair in $\chi$).

### 16.5 Activiteit en concentratie-benadering
In verdunde oplossingen gebruik je vaak:
$$
a_i \approx \frac{c_i}{c^\circ}
$$
(ideaal verdund; $c^\circ$ standaardconcentratie). Dan bouw je:
$$
Q_c=\prod_i \left(\frac{c_i}{c^\circ}\right)^{\nu_i}.
$$
Bij evenwicht:
$$
Q_c = K_c^\circ.
$$

---

## 17) Zuur–base: $K_a$, $K_w$, pH, buffer

### 17.1 Zuurconstante $K_a$
Reactie:
$$
HA + H_2O \rightleftharpoons H_3O^+ + A^-.
$$
Evenwichtsconstante (activiteiten):
$$
K_a=\frac{a_{H_3O^+}\,a_{A^-}}{a_{HA}\,a_{H_2O}}.
$$
Omdat water als oplosmiddel ~ constant is: $a_{H_2O}\approx 1$:
$$
K_a \approx \frac{a_{H_3O^+}\,a_{A^-}}{a_{HA}}.
$$
In ideale verdunde benadering $a_i\approx [i]/c^\circ$ en vaak absorbeer je $c^\circ$ in de definitie:
$$
K_a \approx \frac{[H_3O^+][A^-]}{[HA]}.
$$

### 17.2 Link met $\Delta_r G^\circ$
Uit sectie 13:
$$
K_a^\circ=\exp\!\left(-\frac{\Delta_r G^\circ}{RT}\right).
$$

### 17.3 Waterionproduct
Auto-ionisatie:
$$
2H_2O \rightleftharpoons H_3O^+ + OH^-.
$$
Met $a_{H_2O}\approx 1$:
$$
K_w = a_{H_3O^+}a_{OH^-}\approx [H^+][OH^-].
$$
Bij 25°C is $K_w\approx 10^{-14}$, dus in neutraal water:
$$
[H^+] = [OH^-] = 10^{-7}.
$$

### 17.4 “Snelle” pH-benaderingen
- **Sterk zuur** met beginconcentratie $c_0$:
  $$
  [H^+] \approx c_0 \Rightarrow pH=-\log c_0.
  $$
- **Zwak zuur** $HA$ met $c_0$:
  Stel $x=[H^+]$ in evenwicht, dan $[A^-]=x$, $[HA]=c_0-x$:
  $$
  K_a=\frac{x^2}{c_0-x}\approx \frac{x^2}{c_0}
  \Rightarrow x\approx \sqrt{K_a c_0}.
  $$
  Dus:
  $$
  pH \approx -\log\sqrt{K_a c_0}=\frac12\left(pK_a-\log c_0\right).
  $$

### 17.5 Henderson–Hasselbalch
Van:
$$
K_a=\frac{[H^+][A^-]}{[HA]}
\Rightarrow [H^+] = K_a\frac{[HA]}{[A^-]}.
$$
Neem $-\log$:
$$
pH = pK_a + \log\!\left(\frac{[A^-]}{[HA]}\right).
$$
Met concentraties van buffercomponenten:
$$
pH = pK_a + \log\!\left(\frac{c_{\text{base}}}{c_{\text{zuur}}}\right).
$$

---

## 18) Oplosbaarheid: $K_s$, $Q_s$, link met $S$

### 18.1 $Q_s$ en vaste stof activiteit
Voor:
$$
A_cB_d(s)\rightleftharpoons cA^+(aq)+dB^-(aq),
$$
$$
K_s=\frac{a_{A^+}^c a_{B^-}^d}{a_{A_cB_d(s)}}.
$$
Voor een zuivere vaste stof is $a_{solid}=1$, dus:
$$
K_s = a_{A^+}^c a_{B^-}^d \approx [A^+]^c[B^-]^d.
$$
Daarom:
$$
Q_s=[A^+]^c[B^-]^d.
$$

### 18.2 Link met $\Delta G^\circ$
Algemeen:
$$
K_s=\exp\!\left(-\frac{\Delta G^\circ}{RT}\right).
$$

### 18.3 Relatie met oplosbaarheid $S$
Als de molaire oplosbaarheid $S$ is (mol/L):
$$
[A^+] = cS,\qquad [B^-]=dS.
$$
Invullen:
$$
K_s = (cS)^c(dS)^d.
$$

---
