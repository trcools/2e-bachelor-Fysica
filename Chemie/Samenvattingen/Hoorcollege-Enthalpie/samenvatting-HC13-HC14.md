# HC13–HC14 — Thermodynamica: F, G, μ, activiteiten, evenwicht (samenvatting)

## 1) Hoofdwetten + evolutieprincipe
**1e hoofdwet (conventie les):** warmte en arbeid positief als ze $U_{\text{sys}}$ doen stijgen
\[
dU = \delta q + \delta w
\]

**2e hoofdwet (reversibel):**
\[
dS = \frac{\delta q_{\text{rev}}}{T}
\]

**Evolutie/evenwicht (systeem + omgeving):**
\[
\Delta S_{\text{tot}}=\Delta S_{\text{sys}}+\Delta S_{\text{omg}} \ge 0
\]
Probleem: je moet de omgeving mee modelleren → opgelost door toestandfuncties $F$ en $G$.

---

## 2) Helmholtz vrije energie $F$ (of $A$)
**Definitie**
\[
F \equiv U - TS
\]

**Differentiaal**
\[
dF = dU - T\,dS - S\,dT
\]

**Temperatuurcontrole ($dT=0$):**
\[
dF = dU - T\,dS
\]

**Bij $T,V$-controle + enkel $pV$-arbeid:**
\[
\Delta F \le 0
\]
**Evenwicht:** $F$ minimaal (bij vaste $T,V$).

**Werk-interpretatie (reversibel = “gelijkheidsdenken”):** $\Delta F$ geeft de arbeidsgrens (maximaal haalbare arbeid onder de gekozen randvoorwaarden).

---

## 3) Gibbs vrije enthalpie $G$
**Definitie**
\[
G \equiv H-TS = U+PV-TS
\]

**Algemene differentiaal**
\[
dG = dU + P\,dV + V\,dP - T\,dS - S\,dT
\]

**Voor eenvoudig $pV$-werksysteem (vaste samenstelling):**
\[
dG = -S\,dT + V\,dP
\]

**Bij $T,P$-controle + enkel $pV$-arbeid:**
\[
\Delta G \le 0
\]
**Evenwicht:** $G$ minimaal (bij vaste $T,P$).

**Werkbetekenis (PT):** bij $T,P$ constant is $-\Delta G$ de maximaal haalbare *niet-$pV$* arbeid (bv. elektrische) voor een reversibel proces.

---

## 4) Afgeleiden: T- en P-afhankelijkheid van $G$
Uit
\[
dG=-S\,dT+V\,dP
\]
volgt
\[
\left(\frac{\partial G}{\partial T}\right)_P=-S,\qquad
\left(\frac{\partial G}{\partial P}\right)_T=V
\]
Dus: $S>0$ ⇒ $G$ daalt met stijgende $T$.

---

## 5) Drukafhankelijkheid: ideaal gas vs. gecondenseerde fasen
### Ideaal gas (molaire vorm)
\[
\left(\frac{\partial G_m}{\partial P}\right)_T = V_m = \frac{RT}{P}
\]
Integratie:
\[
G_m(T,P)=G_m^\circ(T)+RT\ln\!\left(\frac{P}{P^\circ}\right)
\]

**Ideaal gasmengsel:** vervang $P$ door partiële druk $P_i$:
\[
G_{m,i}=\mu_i=\mu_i^\circ(T)+RT\ln\!\left(\frac{P_i}{P^\circ}\right)
\]

### Gecondenseerde fasen (vloeistof/vaste stof)
$V_m$ klein ⇒ $V\,\Delta P$ klein ⇒ vaak:
\[
G_m(T,P)\approx G_m^\circ(T)
\]
(drukafhankelijkheid verwaarloosbaar in veel problemen)

---

## 6) Samenstelling → chemische potentiaal $\mu_i$
Bij $T,P$ vast is $G$ homogeen van graad 1 in $n_i$ (verdubbel alles ⇒ $G$ verdubbelt).

**Definitie**
\[
\mu_i \equiv \left(\frac{\partial G}{\partial n_i}\right)_{T,P,n_{j\neq i}}
\]

**Euler (homogeniteit):**
\[
G=\sum_i n_i\mu_i
\]

**Volledige differentiaal (met samenstelling):**
\[
dG=-S\,dT+V\,dP+\sum_i \mu_i\,dn_i
\]

---

## 7) Gibbs–Duhem relatie
Start met:
\[
G=\sum_i n_i\mu_i \Rightarrow dG=\sum_i\mu_i\,dn_i+\sum_i n_i\,d\mu_i
\]
Vergelijk met
\[
dG=-S\,dT+V\,dP+\sum_i\mu_i\,dn_i
\]
⇒
\[
\sum_i n_i\,d\mu_i=-S\,dT+V\,dP
\]

**Bij vaste $T,P$:**
\[
\sum_i n_i\,d\mu_i=0
\]
Dus $\mu_i$ zijn niet onafhankelijk.

---

## 8) Activiteiten: universele vorm voor $\mu$
\[
\mu_i=\mu_i^\circ+RT\ln a_i
\]
waar $a_i$ dimensieloos is.

**Gekende activiteiten:**
- **Zuivere gecondenseerde fase:** $a_i=1 \Rightarrow \mu_i=\mu_i^\circ$
- **Ideaal gas (zuiver of mengsel):** $a_i=P_i/P^\circ$

**Valkuil:** in $\ln$ mag enkel dimensieloos:
\[
\ln\!\left(\frac{P}{P^\circ}\right)\ \text{niet}\ \ln(P)
\]

---

## 9) Fase-evenwicht = gelijkheid van chemische potentialen
Voor component $A$ in twee fasen $\alpha,\beta$ bij $T,P$ constant:
\[
\mu_A^{(\alpha)}=\mu_A^{(\beta)}
\]
Materie gaat spontaan naar de fase met lagere $\mu$ tot gelijkheid.

---

## 10) Dampdruk (voorbeeld water) via $\mu$ en activiteit
Evenwicht vloeistof ↔ gas:
\[
\mu_g=\mu_l
\]
Invullen:
\[
\mu_g^\circ+RT\ln\!\left(\frac{P^*}{P^\circ}\right)=\mu_l^\circ
\]
Dus
\[
\ln\!\left(\frac{P^*}{P^\circ}\right)=-\frac{\mu_g^\circ-\mu_l^\circ}{RT}
=-\frac{\Delta G_{\text{vap}}^\circ}{RT}
\]
\[
P^*=P^\circ\exp\!\left(-\frac{\Delta G_{\text{vap}}^\circ}{RT}\right)
\]

---

## 11) Temperatuursafhankelijkheid van evenwichten
\[
\Delta G^\circ(T)=\Delta H^\circ(T)-T\Delta S^\circ(T)
\]

Via Hess-cyclus met opwarmen/afkoelen:
\[
\frac{d(\Delta H^\circ)}{dT}=\Delta C_p,\qquad
\frac{d(\Delta S^\circ)}{dT}=\frac{\Delta C_p}{T}
\]

Vaak eerste benadering (beperkt $T$-interval): $\Delta H^\circ,\Delta S^\circ$ ongeveer constant.

**Kookpunt bij 1 bar (evenwicht vloeistof–gas):**
\[
\Delta G^\circ(T_b)=0 \Rightarrow T_b\approx \frac{\Delta H^\circ}{\Delta S^\circ}
\]

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
