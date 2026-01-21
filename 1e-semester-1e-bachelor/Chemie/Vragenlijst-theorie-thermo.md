# Thermodynamica — Vragenlijst (afleidingen)

## 1) Eerste hoofdwet + arbeid
- [ ] Leid de **eerste hoofdwet** af/verklaar de notatie: waarom $dU$ maar $\delta q,\delta w$? Start van $U$ als toestandsfunctie en $q,w$ als procesgrootheden.
- [ ] Leid de vorm af:
  $$
  \delta w=\delta w_{PV}+\delta w_{\text{nuttig}}
  $$
  en toon waarom bij P–V arbeid:
  $$
  \delta w_{PV}=-P_{\text{omg}}\,dV.
  $$
  (Begin bij $-F\,dx$ en koppel naar $P\,dV$.)
- [ ] Combineer dit tot:
  $$
  dU=\delta q-P_{\text{omg}}\,dV+\delta w_{\text{nuttig}}.
  $$
  Welke assumpties zitten hier impliciet?
- [ ] Reversibel vs irreversibel: geef een thermodynamische reden waarom reversibel “maximale arbeid” kan opleveren.

## 2) Warmtecapaciteiten $C_V$ en $C_P$
- [ ] Toon dat in een **gesloten systeem** met $\delta w_{\text{nuttig}}=0$ en $V=$ const:
  $$
  dU=\delta q=C_V\,dT,
  \qquad
  C_V=\left(\frac{\partial U}{\partial T}\right)_V.
  $$
- [ ] Toon dat bij $P=$ const:
  $$
  dU+P_{\text{omg}}\,dV=\delta q=C_P\,dT.
  $$
  Waar komt die extra term fysisch vandaan?

## 3) Enthalpie $H$ en differentiaal
- [ ] Vertrek van de definitie:
  $$
  H=U+P_{\text{omg}}V
  $$
  en leid $dH$ af. Welke term valt weg als $P_{\text{omg}}$ constant is?
- [ ] Leid af:
  $$
  C_P=\left(\frac{\partial H}{\partial T}\right)_P.
  $$

## 4) Hess + standaardenthalpieën
- [ ] Bewijs (kort maar strak) waarom de **wet van Hess** volgt uit “enthalpie is toestandsfunctie”.
- [ ] Leid de formule voor standaard reactie-enthalpie af uit vormingsenthalpieën:
  $$
  \Delta_r H^\circ=\sum_P \nu_P H_{m,P}^\circ-\sum_R \nu_R H_{m,R}^\circ.
  $$

## 5) Entropie + tweede hoofdwet
- [ ] Leid voor reversibele processen:
  $$
  dS=\frac{\delta q_\text{rev}}{T}
  $$
  en formuleer de irreversibele ongelijkheid:
  $$
  dS>\frac{\delta q}{T}.
  $$
  Wat is de *inhoud* van dat “>”?
- [ ] Toon dat bij constante $P$ en $T$ en zonder nuttige arbeid:
  $$
  dS=\frac{dH}{T}.
  $$

## 6) Gibbs vrije energie $G$: definitie, spontaniteit, differentiaal
- [ ] Vertrek van
  $$
  G=H-TS
  $$
  en leid (bij constante $T$) af:
  $$
  dG=dH-T\,dS.
  $$
- [ ] Toon het spontaniteitscriterium bij $P,T$ constant:
  - waarom $dG\le \delta w_{\text{nuttig}}$?
  - en waarom bij enkel P–V arbeid $dG\le 0$?
- [ ] Leid de “fundamentele identiteit” af:
  $$
  dG=V\,dP-S\,dT
  $$
  door te starten van $G=U+PV-TS$ en $dU=T\,dS-P\,dV$.
- [ ] Toon de partiële afgeleide:
  $$
  \left(\frac{\partial G}{\partial T}\right)_P=-S.
  $$

## 7) Temperatuurafhankelijkheid via $G/T$
- [ ] Vertrek van $G=U+PV-TS$ en leid af:
  $$
  \left(\frac{\partial (G/T)}{\partial T}\right)_P=-\frac{H}{T^2}.
  $$
  Let op: je moet expliciet de algebra met $G/T$ uitschrijven.

## 8) Drukafhankelijkheid van $G$
- [ ] Toon eerst algemeen:
  $$
  \left(\frac{\partial G}{\partial P}\right)_T=V.
  $$
- [ ] Gebruik de ideale gaswet om af te leiden:
  $$
  \left(\frac{\partial G}{\partial P}\right)_T=\frac{nRT}{P}
  \Rightarrow
  G(P)=G(P^\circ)+nRT\ln\!\left(\frac{P}{P^\circ}\right).
  $$
  (Schrijf de integraalstappen uit.)

## 9) Helmholtz-energie $A$
- [ ] Vertrek van
  $$
  A=U-TS
  $$
  en leid af:
  $$
  dA=-S\,dT-P\,dV
  $$
  (combineer met $dU=T\,dS-P\,dV$.)
- [ ] Toon:
  $$
  \left(\frac{\partial A}{\partial T}\right)_V=-S,
  \qquad
  \left(\frac{\partial A}{\partial V}\right)_T=-P.
  $$

## 10) Chemische potentiaal + Gibbs–Duhem
- [ ] Definieer chemische potentiaal en leid af waarom:
  $$
  \mu_i=\left(\frac{\partial G}{\partial n_i}\right)_{P,T,n_{j\ne i}}.
  $$
- [ ] Toon dat (voor mengsels):
  $$
  G=\sum_i n_i\mu_i.
  $$
  Welke aanname over “extensiviteit” gebruik je?
- [ ] Leid **Gibbs–Duhem** af door $G=\sum_i n_i\mu_i$ te differentiëren en te vergelijken met
  $$
  dG=-S\,dT+V\,dP+\sum_i \mu_i\,dn_i:
  \qquad
  S\,dT - V\,dP + \sum_i n_i\,d\mu_i = 0.
  $$
- [ ] Leid de “gevallen” af/verklaar:
  - zuivere stof: $\mu_A=\mu_A^\circ=G_{m,A}$
  - ideale gassen: koppel $\mu$ aan $\ln P$ (activiteit via $P/P^\circ$).

## 11) Fase-evenwicht: Clausius–Clapeyron + fasenregel
- [ ] Herleid de geïntegreerde **Clausius–Clapeyron**-vorm:
  $$
  \ln\!\left(\frac{P_2}{P_1}\right)= -\frac{\Delta_{\alpha\to\beta}H}{R}\left(\frac{1}{T_2}-\frac{1}{T_1}\right),
  $$
  en herschrijf naar $P_2$ expliciet. (Schrijf de integraalstappen uit.)
- [ ] Leid/verklaar de **fasenregel**:
  $$
  f = 2 + c - p,
  $$
  en geef een fysische interpretatie van $f$.

## 12) Chemisch evenwicht: vorderingsgraad $\xi$ en $\Delta_r G$
- [ ] Start van:
  $$
  dG=\sum_P \mu_P\,dn_P + \sum_R \mu_R\,dn_R
  $$
  en definieer $\xi$ zodat $dn_P=\nu_P\,d\xi$ en $dn_R=-\nu_R\,d\xi$. Leid af:
  $$
  dG=\Delta_r G\,d\xi,
  \qquad
  \Delta_r G=\sum_P \nu_P \mu_P-\sum_R \nu_R \mu_R.
  $$
- [ ] Toon dat evenwicht $\Delta_r G=0$ betekent, en koppel de tekenregels $(<0,=0,>0)$ aan richting van spontaniteit.

## 13) Activiteiten, reactiequotiënt $Q$ en evenwichtsconstante $K$
- [ ] Vertrek van:
  $$
  \mu_A=\mu_A^\circ+RT\ln(a_A)
  $$
  en leid af:
  $$
  \Delta G=\Delta G^\circ + RT\ln(Q_a),
  \quad
  Q_a=\frac{\prod_P a_P^{\nu_P}}{\prod_R a_R^{\nu_R}}.
  $$
- [ ] Zet $\Delta_r G=0$ bij evenwicht en leid af:
  $$
  K_a^\circ=\exp\!\left(-\frac{\Delta_r G^\circ}{RT}\right)
  \quad\text{en}\quad
  Q_a=K_a^\circ\ \text{bij evenwicht}.
  $$

## 14) Ideale gassen: $Q_P$, $K_P^\circ$ en $K_\chi$
- [ ] Toon waarom in ideale gasmengsels:
  $$
  a=\frac{P}{P^\circ}
  $$
  en leid de vorm van $Q_P$ af met partiële drukken.
- [ ] Leid de omzetting af naar molfracties:
  $$
  Q_\chi=\frac{\prod_P \chi_P^{\nu_P}}{\prod_R \chi_R^{\nu_R}},
  \qquad
  K_\chi = K_P^\circ\left(\frac{P^\circ}{P_{\text{tot}}}\right)^{\Delta\nu}.
  $$
  (Belangrijk: toon waar $\Delta\nu$ precies vandaan komt.)

## 15) Massabalans met $\alpha$ (gasreacties)
- [ ] Definieer $\alpha$ zoals in jouw samenvatting en leid de algemene molbalans af:
  $$
  n_R=n_{R,0}-\alpha\,n_0\,\nu_R,
  \qquad
  n_P=\alpha\,n_0\,\nu_R \ (\text{met juiste stoichiometrie}).
  $$
- [ ] Leid de molfracties af als functie van $\alpha$ en toon hoe je $Q_\chi(\alpha)$ opbouwt.

## 16) Oplossingen: roostermodel → $S$, $H$, $G$, $\mu$
- [ ] Entropie in oplossing: start van “aantal schikkingsmogelijkheden” in een rooster met $N$ plaatsen en leid de vorm af:
  $$
  S \sim -R(\chi_A\ln\chi_A+\chi_B\ln\chi_B)
  $$
  (toon waar de log-termen vandaan komen).
- [ ] Enthalpie in oplossing: leid de (molaire) vorm af die in je samenvatting staat met interacties en $\chi_A\chi_B$-term.
- [ ] Combineer naar Gibbs:
  $$
  G=H-TS
  $$
  en reproduceer de expliciete $G(\chi_A,\chi_B)$-vorm uit je samenvatting.
- [ ] Leid de chemische potentialen af via partiële afgeleiden:
  $$
  \mu_A=\left(\frac{\partial G}{\partial n_A}\right)_{P,T,n_B},\quad
  \mu_B=\left(\frac{\partial G}{\partial n_B}\right)_{P,T,n_A},
  $$
  en toon dat je log-termen $\ln(\chi_A)$, $\ln(\chi_B)$ krijgt.
- [ ] Activiteit in oplossing: leid/verklaar waarom:
  $$
  a_B \approx [B]=\frac{c_B}{c^\circ}\approx \frac{\rho_A\chi_B}{M_A}.
  $$
- [ ] Leid $Q_c$ af:
  $$
  Q_c=\frac{\prod_P [P]^{\nu_P}}{\prod_R [R]^{\nu_R}},
  \qquad
  Q_c=K_c^\circ\ \text{bij evenwicht}.
  $$

## 17) Zuur–base: $K_a$, $K_w$, pH-formules, buffer
- [ ] Leid de evenwichtsrelatie voor een zuur in water af en verklaar waarom $a_{\mathrm{H_2O}}=1$:
  $$
  K_a^\circ=\frac{[\mathrm{H_3O^+}][\mathrm{A^-}]}{[HA]}.
  $$
- [ ] Leid af:
  $$
  K_a^\circ=\exp\!\left(-\frac{\Delta_r G^\circ}{RT}\right).
  $$
- [ ] Auto-ionisatie water: leid de $K_w$-relatie af en link naar neutraal water:
  $$
  K_w=[\mathrm{H^+}][\mathrm{OH^-}]=10^{-14},\quad [\mathrm{H^+}]=[\mathrm{OH^-}]=10^{-7}.
  $$
- [ ] Leid de “snelle” pH-benaderingen af voor:
  - sterke zuren/basen: $pH=-\log c_0$
  - zwakke zuren/basen: de $\tfrac12$-formules (toon de approximatie-stap!)
- [ ] Buffer: leid Henderson–Hasselbalch af:
  $$
  pH=pK_a+\log\!\left(\frac{c_{\text{base}}}{c_{\text{zuur}}}\right),
  $$
  vanuit $K_a$ en massabalans/elektroneutraliteit.

## 18) Oplosbaarheid: $K_s$, $Q_s$, relatie met oplosbaarheid $S$
- [ ] Voor $A_cB_d(s)\rightleftharpoons cA^+(aq)+dB^-(aq)$: leid af
  $$
  Q_s=[A^+]^c[B^-]^d
  $$
  en leg uit waarom vaste stof activiteit $=1$ is.
- [ ] Toon bij evenwicht:
  $$
  Q_s=K_s=\exp\!\left(-\frac{\Delta G^\circ}{RT}\right).
  $$
- [ ] Leid de relatie tussen $K_s$ en oplosbaarheid $S$ af (met stoichiometrie!):
  $$
  [A^+]=cS,\quad [B^-]=dS \Rightarrow K_s=(cS)^c(dS)^d.
  $$
