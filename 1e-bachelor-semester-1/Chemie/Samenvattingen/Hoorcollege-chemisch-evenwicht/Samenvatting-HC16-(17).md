# Hoorcollege (vermoedelijk 18) — Oplossingen, activiteiten, oplosbaarheid, Ksp, Henry

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
