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
E^\circ_\text{cel} = E^\circ_{\text{kathode (reductie)}} - E^\circ_{\text{anode (reductie)}}
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

### Brandstofcel (H₂/O₂, $E_\text{cel} \approx 1{,}2\ \text{V}$)

- Reagentia worden continu aangevoerd.  
- Nettoreactie:  
  $$
  2H_2 + O_2 \rightarrow 2H_2O
  $$
- Efficiëntie-idee: groot deel van theoretische $\Delta G$ → elektrische energie.  
- Nadelen: opslag van reagentia en dure elektroden.

### Li-ion (conceptueel herkennen)

- **Anode: grafiet**  
- **Kathode: LiCoO₂**  
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
