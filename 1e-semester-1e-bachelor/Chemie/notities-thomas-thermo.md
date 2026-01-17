# Notities thermo (net) — alle afleidingen uitgeschreven (markdown)

> Conventie (zoals typisch in chemie):  
> $\delta q>0$ als warmte *naar* het systeem gaat.  
> $\delta w>0$ als arbeid *op* het systeem wordt verricht.  
> (Dus expansiewerk door het systeem zelf is negatief.)

---

## 0) Warmtecapaciteit(en)

### Definitie
Warmtecapaciteit $C$ (van een systeem) is:
$$
C=\frac{\delta q}{dT}
$$
Je ziet vaak “molaire” warmtecapaciteit $C_m$ (per mol) en “soortelijke” (per massa).

---

## 1) pV-arbeid: $\delta w_{pV}=-p_\text{ext}\,dV$

Neem een zuiger met doorsnede $A$.

- Kracht door de omgeving: $F_\text{ext}=p_\text{ext}A$
- Verplaatsing: $dx$
- Volumeverandering: $dV=A\,dx$

Arbeid op het systeem door de omgeving:
$$
\delta w = -F_\text{ext}\,dx = -p_\text{ext}A\,dx=-p_\text{ext}\,dV
$$

Voor een eindig proces:
$$
w_{pV}=-\int_{V_1}^{V_2}p_\text{ext}\,dV
$$

**Reversibel** (quasi-statisch): $p_\text{ext}=p_\text{sys}$ op elk moment, dus
$$
w_{pV,\text{rev}}=-\int_{V_1}^{V_2}p\,dV
$$

---

## 2) Eerste hoofdwet: $dU=\delta q+\delta w$

De interne energie $U$ is een **toestandsfunctie** (pad-onafhankelijk), terwijl warmte en arbeid pad-afhankelijk zijn:
$$
dU=\delta q+\delta w
$$

Als de enige arbeid pV-arbeid is:
$$
dU=\delta q - p_\text{ext}\,dV
$$

Bij reversibel pV-proces:
$$
dU=\delta q_\text{rev}-p\,dV
$$

---

## 3) Warmtecapaciteiten $C_V$ en $C_p$

### 3.1) Isochoor: $dV=0$
Bij $V=\text{cte}$ en (typisch) enkel pV-arbeid:
$$
dU=\delta q_V
$$
Definieer:
$$
C_V=\left(\frac{\partial U}{\partial T}\right)_V
$$
Dan:
$$
\delta q_V=C_V\,dT
$$

### 3.2) Isobaar: $dp=0$ en enthalpie (vooruitblik)
Bij $p=\text{cte}$ is $\delta q_p$ meestal niet gelijk aan $dU$ omdat het systeem ook pV-arbeid levert/ontvangt.
Daarom introduceren we $H$.

---

## 4) Enthalpie: $H=U+pV$ en $q_p=\Delta H$

### 4.1) Definitie en differentiaal
$$
H:=U+pV
$$
Differentieer:
$$
dH=dU+p\,dV+V\,dp
$$

### 4.2) Enkel pV-arbeid
Met $dU=\delta q - p\,dV$ (bij isobaar proces is $p_\text{ext}=p$):
$$
dH=(\delta q-p\,dV)+p\,dV+V\,dp=\delta q+V\,dp
$$
Bij **constante druk** ($dp=0$):
$$
dH=\delta q_p
\quad\Rightarrow\quad
\Delta H=q_p
$$

### 4.3) Definitie van $C_p$
Bij $p=\text{cte}$:
$$
\delta q_p=C_p\,dT
$$
en omdat $\delta q_p=dH$:
$$
C_p=\left(\frac{\partial H}{\partial T}\right)_p
$$

---

## 5) Reactie-enthalpie en Hess

### 5.1) Reactie-enthalpie via molaire enthalpieën
Voor reactie $\sum_R \nu_R R \to \sum_P \nu_P P$:
$$
\Delta_r H=\sum_P \nu_P H_{m,P}-\sum_R \nu_R H_{m,R}
$$

### 5.2) Hess’ wet (bewijs-idee)
Omdat $H$ een toestandsfunctie is:
$$
\Delta H_{1\to 2}=H(2)-H(1)
$$
en voor tussenstap $A$:
$$
\Delta H_{1\to 2}=\Delta H_{1\to A}+\Delta H_{A\to 2}
$$
Dus reacties mogen in stappen opgesplitst worden: de totale $\Delta_rH$ is de som.

### 5.3) Standaard reactie-enthalpie uit vormingsenthalpieën
Met standaard (molaire) vormingsenthalpieën $ \Delta_f H_m^\circ$:
$$
\Delta_r H^\circ=\sum_P \nu_P\,\Delta_f H_{m,P}^\circ-\sum_R \nu_R\,\Delta_f H_{m,R}^\circ
$$

---

## 6) Tweede hoofdwet: entropie en Clausius

### 6.1) Definitie (reversibel)
$$
dS=\frac{\delta q_\text{rev}}{T}
$$

### 6.2) Clausius-ongelijkheid
Voor een willekeurig proces:
$$
dS \ge \frac{\delta q}{T}
$$
met gelijkheid enkel bij reversibel.

### 6.3) Geïsoleerd universum
Voor “systeem + omgeving” geïsoleerd:
$$
\Delta S_\text{univ}=\Delta S_\text{sys}+\Delta S_\text{omg}\ge 0
$$
Evenwicht (macro): $\Delta S_\text{univ}=0$.

---

## 7) Entropie van een ideaal gas — belangrijke afleiding

### 7.1) Isotherme reversibele expansie van 1 mol ideaal gas
Voor ideaal gas: $U=U(T)$, dus bij isotherm $dU=0$.

Eerste hoofdwet:
$$
0=dU=\delta q_\text{rev}+\delta w_\text{rev}
$$
met $\delta w_\text{rev}=-p\,dV$:
$$
\delta q_\text{rev}=p\,dV
$$
Ideale gaswet (1 mol): $p=\frac{RT}{V}$, dus
$$
\delta q_\text{rev}=\frac{RT}{V}\,dV
$$
Entropie:
$$
dS=\frac{\delta q_\text{rev}}{T}=\frac{R}{V}\,dV
$$
Integreren:
$$
\Delta S = R\ln\!\left(\frac{V_2}{V_1}\right)
$$
Voor $n$ mol:
$$
\Delta S = nR\ln\!\left(\frac{V_2}{V_1}\right)
$$

### 7.2) Temperatuurverandering (met warmtecapaciteit)
Als langs een reversibel pad $\delta q_\text{rev}=nC\,dT$:
$$
dS=\frac{nC\,dT}{T}
\Rightarrow
\Delta S=nC\ln\!\left(\frac{T_2}{T_1}\right)
$$
Kies $C=C_V$ (isochoor) of $C=C_p$ (isobaar) volgens proces.

---

## 8) Boltzmann en mengentropie

### 8.1) Boltzmann
$$
S=k_B\ln W
$$
waar $W$ het aantal microtoestanden is compatibel met de macrotoestand.

### 8.2) Ideale menging (twee componenten) — eindresultaat
Voor ideale menging van componenten met molfracties $x_i$:
$$
\Delta S_\text{mix}=-nR\sum_i x_i\ln x_i
$$
(De kern is combinatoriek: het aantal permutaties $W$ van de deeltjes over posities levert log-termen.)

---

## 9) Derde hoofdwet en absolute entropieën

### 9.1) Derde hoofdwet (Nernst)
Entropie van een perfect kristal bij $T=0$ K:
$$
S(0)=0
$$

### 9.2) Absolute entropie uit calorimetrie (schets van de opbouw)
Voor een stof:
$$
S(T)=\int_0^T \frac{C_p(T')}{T'}\,dT' + \sum_k \frac{\Delta H_k}{T_k}
$$
waar de som over faseovergangen $k$ loopt (smelten, verdampen, …) met overgangstemperatuur $T_k$.

---

## 10) Vrije enthalpie (Gibbs): afleiding en betekenis

### 10.1) Doel: een spontaniteitscriterium met alleen systeemeigenschappen bij $p,T$-controle
Neem “systeem + omgeving” als (bij benadering) geïsoleerd.

Voor de omgeving als warmtereservoir op temperatuur $T$:
$$
dS_\text{omg}=\frac{\delta q_\text{omg}}{T}
$$
en warmtebalans: $\delta q_\text{omg}=-\delta q_\text{sys}$.

Dus:
$$
dS_\text{univ}=dS_\text{sys}+dS_\text{omg}
=dS_\text{sys}-\frac{\delta q_\text{sys}}{T}\ge 0
$$
Vermenigvuldig met $-T$:
$$
-T\,dS_\text{univ} \le -T\,dS_\text{sys}+\delta q_\text{sys}
$$
Gebruik eerste hoofdwet $\delta q_\text{sys}=dU-\delta w$:
$$
-T\,dS_\text{univ} \le -T\,dS_\text{sys}+dU-\delta w
$$
Splits arbeid: $\delta w=\delta w_{pV}+\delta w_\text{nuttig}$ met $\delta w_{pV}=-p\,dV$ (bij $p$-controle):
$$
-T\,dS_\text{univ} \le -T\,dS_\text{sys}+dU+p\,dV-\delta w_\text{nuttig}
$$
Herschik:
$$
dU+p\,dV-T\,dS_\text{sys} \le \delta w_\text{nuttig}
$$
Definieer nu:
$$
G:=U+pV-TS
$$
Differentieer bij constante $p$ en $T$:
$$
dG=dU+p\,dV-T\,dS
$$
Dus:
$$
dG \le \delta w_\text{nuttig}
$$
Bij **geen nuttige arbeid** (enkel pV-arbeid): $\delta w_\text{nuttig}=0$, dus
$$
dG\le 0
$$
**Spontaan** bij $p,T$ constant $\Leftrightarrow \Delta G<0$.

### 10.2) Interpretatie
Bij $p,T$ constant is $-\Delta G$ de maximaal haalbare nuttige arbeid (reversibel):
$$
w_{\text{nuttig,max}} = -\Delta G
$$

---

## 11) Fundamentele differentiaal van $G$ en partiële afgeleiden

### 11.1) Startpunt
$$
G=U+pV-TS
$$
Differentieer:
$$
dG=dU+p\,dV+V\,dp-T\,dS-S\,dT
$$

Voor een reversibel proces met enkel pV-arbeid:
$$
dU=\delta q_\text{rev}+\delta w_\text{rev}=T\,dS-p\,dV
$$
Invullen:
$$
dG=(T\,dS-p\,dV)+p\,dV+V\,dp-T\,dS-S\,dT
$$
Dus:
$$
dG=-S\,dT+V\,dp
$$
(voor vaste samenstelling).

### 11.2) Conclusies
Omdat $G=G(T,p)$ (bij vaste samenstelling):
$$
dG=\left(\frac{\partial G}{\partial T}\right)_p dT+\left(\frac{\partial G}{\partial p}\right)_T dp
$$
Vergelijk met $dG=-S\,dT+V\,dp$:
$$
\left(\frac{\partial G}{\partial T}\right)_p=-S,
\qquad
\left(\frac{\partial G}{\partial p}\right)_T=V
$$

---

## 12) Drukafhankelijkheid van $G$ voor ideale gassen

Bij constante $T$:
$$
\left(\frac{\partial G}{\partial p}\right)_T=V
$$
Voor ideaal gas: $V=\frac{nRT}{p}$, dus
$$
dG=\frac{nRT}{p}\,dp
$$
Integreren van $p^\circ$ naar $p$:
$$
G(T,p)-G(T,p^\circ)=nRT\ln\!\left(\frac{p}{p^\circ}\right)
$$
Dus:
$$
G(T,p)=G^\circ(T)+nRT\ln\!\left(\frac{p}{p^\circ}\right)
$$
Molaire vorm (zuivere stof): $g=\mu$:
$$
\mu(T,p)=\mu^\circ(T)+RT\ln\!\left(\frac{p}{p^\circ}\right)
$$

---

## 13) Chemische potentiaal en samenstellingsafhankelijkheid

Voor variabele samenstelling:
$$
dG=-S\,dT+V\,dp+\sum_i \mu_i\,dn_i
$$
Definitie:
$$
\mu_i:=\left(\frac{\partial G}{\partial n_i}\right)_{T,p,n_{j\ne i}}
$$

### 13.1) Activiteit
Algemeen schrijf je:
$$
\mu_i=\mu_i^\circ+RT\ln a_i
$$
- ideaal gas: $a_i=\frac{p_i}{p^\circ}$
- zuivere vaste stof / zuivere vloeistof: $a_i=1$

---

## 14) Fase-evenwicht: $\mu^\alpha=\mu^\beta$ (korte variatie-redenering)

Neem één component $A$ in twee fasen $\alpha,\beta$ bij vaste $T,p$.
Laat een kleine hoeveelheid $dn$ van $\alpha$ naar $\beta$ gaan:
$$
dn^\alpha=-dn,\qquad dn^\beta=+dn
$$
Dan:
$$
dG=\mu_A^\alpha\,dn^\alpha+\mu_A^\beta\,dn^\beta
= -\mu_A^\alpha\,dn+\mu_A^\beta\,dn
=(\mu_A^\beta-\mu_A^\alpha)\,dn
$$
Evenwicht betekent $dG=0$ voor elke (toelaatbare) $dn$, dus:
$$
\mu_A^\alpha=\mu_A^\beta
$$

---

## 15) Clapeyron en Clausius–Clapeyron

### 15.1) Clapeyron
Langs de evenwichtslijn tussen $\alpha$ en $\beta$ geldt $\mu^\alpha=\mu^\beta$.
Differentieer:
$$
d\mu^\alpha=d\mu^\beta
$$
Voor een zuivere stof:
$$
d\mu=-\bar S\,dT+\bar V\,dp
$$
Dus:
$$
-\bar S^\alpha dT+\bar V^\alpha dp=-\bar S^\beta dT+\bar V^\beta dp
$$
Herschik:
$$
(\bar V^\beta-\bar V^\alpha)\,dp=(\bar S^\beta-\bar S^\alpha)\,dT
$$
Dus:
$$
\frac{dp}{dT}=\frac{\Delta \bar S}{\Delta \bar V}
$$
Met $\Delta \bar S=\frac{\Delta \bar H}{T}$:
$$
\frac{dp}{dT}=\frac{\Delta \bar H}{T\,\Delta \bar V}
$$

### 15.2) Clausius–Clapeyron (als $\beta$ gas is, ideaal)
Voor verdamping: $\Delta \bar V\approx \bar V_g$ en $\bar V_g=\frac{RT}{p}$.
Dan:
$$
\frac{dp}{dT}\approx \frac{\Delta \bar H}{T}\frac{p}{RT}
=\frac{p\,\Delta \bar H}{R\,T^2}
$$
Dus:
$$
\frac{1}{p}\,dp=\frac{\Delta \bar H}{R}\frac{1}{T^2}\,dT
$$
Integreer (neem $\Delta \bar H$ constant):
$$
\ln\!\left(\frac{p_2}{p_1}\right)
=-\frac{\Delta \bar H}{R}\left(\frac{1}{T_2}-\frac{1}{T_1}\right)
$$

---

## 16) Reactie-evolutie: $dG=\Delta_r G\,d\xi$

Voor reactie met stoichiometrische coëfficiënten $\nu_i$ (positief voor producten, negatief voor reagentia):
$$
dn_i=\nu_i\,d\xi
$$
Bij constante $T,p$:
$$
dG=\sum_i \mu_i\,dn_i=\sum_i \mu_i \nu_i\,d\xi
$$
Definieer reactie-Gibbs:
$$
\Delta_r G:=\sum_i \nu_i \mu_i
$$
Dan:
$$
dG=\Delta_r G\,d\xi
$$
Evenwicht (minimum van $G$) geeft:
$$
\Delta_r G=0
$$

---

## 17) $\Delta_r G=\Delta_r G^\circ+RT\ln Q$ en $K=\exp(-\Delta_rG^\circ/RT)$

Gebruik $\mu_i=\mu_i^\circ+RT\ln a_i$:
$$
\Delta_r G=\sum_i \nu_i\mu_i
=\sum_i \nu_i\mu_i^\circ+RT\sum_i \nu_i\ln a_i
$$
Definieer:
$$
\Delta_r G^\circ:=\sum_i \nu_i\mu_i^\circ
$$
en
$$
Q:=\prod_i a_i^{\nu_i}
\quad\Rightarrow\quad
\ln Q=\sum_i \nu_i\ln a_i
$$
Dus:
$$
\Delta_r G=\Delta_r G^\circ+RT\ln Q
$$
In evenwicht $\Delta_r G=0$ en $Q=K$:
$$
0=\Delta_r G^\circ+RT\ln K
\Rightarrow
K=\exp\!\left(-\frac{\Delta_r G^\circ}{RT}\right)
$$

---

## 18) (Handig) Link tussen $\Delta S_\text{univ}$ en $\Delta G_\text{sys}$ bij $p,T$ constant

Voor een proces bij $p,T$ constant met omgeving als reservoir:
$$
\Delta S_\text{omg}=-\frac{q_p}{T}
$$
en $q_p=\Delta H_\text{sys}$, dus:
$$
\Delta S_\text{univ}=\Delta S_\text{sys}-\frac{\Delta H_\text{sys}}{T}
$$
Vermenigvuldig met $-T$:
$$
-T\Delta S_\text{univ}=\Delta H_\text{sys}-T\Delta S_\text{sys}=\Delta G_\text{sys}
$$
Dus:
$$
\Delta S_\text{univ}=-\frac{\Delta G_\text{sys}}{T}
$$

---

## 19) Extra’s die vaak in dezelfde notities staan (zuur–base / oplosbaarheid)

### 19.1) Zuurconstante
Voor $HA+H_2O\rightleftharpoons A^-+H_3O^+$ en $a_{H_2O}\approx 1$:
$$
K_a=\frac{[H_3O^+][A^-]}{[HA]}
$$

### 19.2) Waterconstante
Voor $2H_2O\rightleftharpoons H_3O^++OH^-$:
$$
K_w=[H_3O^+][OH^-]
$$

### 19.3) Oplosbaarheidsproduct
Voor $A_mB_n(s)\rightleftharpoons mA^{z+}+nB^{z-}$ (vaste stof: activiteit 1):
$$
K_{sp}=[A^{z+}]^m[B^{z-}]^n
$$
Als de molaire oplosbaarheid $s$ is in zuiver water:
$$
[A^{z+}]=ms,\quad [B^{z-}]=ns
\Rightarrow
K_{sp}=(ms)^m(ns)^n
$$

---

# Einde
