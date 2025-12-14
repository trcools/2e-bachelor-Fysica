# HC21 – Elektrochemie (examengerichte samenvatting)

## 1) De kernidee (wat je *altijd* moet kunnen)
Een **galvanische (volta-)cel** zet een **spontane redoxreactie** om in **elektrische arbeid** via elektronen die door een uitwendig circuit lopen.

Spontaniteit-criteria (cruciaal):
- **Spontaan:** $E_\text{cel} > 0 \;\Leftrightarrow\; \Delta G < 0$
- **Evenwicht / cel “plat”:** $E_\text{cel}=0 \;\Leftrightarrow\; \Delta G=0$ en dan geldt $Q=K$

**Anode/Kathode (altijd examenvalkuil):**
- **Anode = oxidatie**
- **Kathode = reductie**
- In een **galvanische** cel: anode is typisch **negatief**, kathode **positief**.

## 2) Opbouw van een galvanische cel (Daniell/Zn–Cu is het archetype)
- 2 **halfcellen** (elk een elektrode + elektrolyt).
- **Zoutbrug** (bv. KNO₃ of NH₄NO₃ in gel) zorgt voor ionenmigratie zodat elke halfcel elektrisch neutraal blijft.
- **Uitwendig circuit**: elektronen lopen van anode → kathode.

**Actieve vs inerte elektroden**
- Actief: de elektrode doet mee in de redox (bv. Zn(s), Cu(s)).
- Inert: geleidt enkel e⁻, doet niet mee (bv. Pt(s), grafiet).

## 3) Elektrodepotentialen en celspanning (rekenrecept)
**Standaardreductiepotentialen $E^\circ$** zijn gedefinieerd t.o.v. de **standaard-waterstofelektrode (SHE)** met $E^\circ = 0.000\ \text{V}$.

Belangrijkste formule:
\[
E^\circ_\text{cel} = E^\circ_{\text{kathode (reductie)}} - E^\circ_{\text{anode (reductie)}}
\]

**Examenschema om $E^\circ_\text{cel}$ te vinden**
1. Schrijf beide halfreacties als **reducties** (zoals in de tabel met $E^\circ$).
2. De halfreactie met **grootste $E^\circ$** wordt de **kathode** (gaat effectief als reductie).
3. De andere wordt **anode** (gaat effectief als oxidatie; teken keert om in je reactie, maar je gebruikt in de formule nog steeds het *reductie*-$E^\circ$).
4. Bereken $E^\circ_\text{cel}$ met de formule hierboven.
5. Balanceer elektronen voor de totale reactie, maar let op:
   - **Je vermenigvuldigt $E^\circ$ nooit met stoichiometrische factoren.**

## 4) Link met thermodynamica: arbeid, vrije energie en evenwicht
Elektrische arbeid en vrije energie hangen samen via:
\[
\Delta G = -nF E_\text{cel}
\]
waar $n$ = aantal uitgewisselde elektronen en $F$ = Faraday-constante.

Standaardrelaties (super-examenwaardig):
\[
\Delta G^\circ = -nF E^\circ_\text{cel}
\]
\[
\Delta G^\circ = -RT\ln K
\]
\[
E^\circ_\text{cel}=\frac{RT}{nF}\ln K
\]

Interpretatie in één oogopslag:
- $\Delta G^\circ<0 \Rightarrow K>1 \Rightarrow E^\circ_\text{cel}>0$ (spontaan in standaardcondities)
- $\Delta G^\circ=0 \Rightarrow K=1 \Rightarrow E^\circ_\text{cel}=0$ (evenwicht)
- $\Delta G^\circ>0 \Rightarrow K<1 \Rightarrow E^\circ_\text{cel}<0$ (niet spontaan)

## 5) Nernstvergelijking: niet-standaard condities
Als concentraties/drukken afwijken van standaardcondities, gebruik je **Nernst**:
\[
E = E^\circ - \frac{RT}{nF}\ln Q
\]
Bij $25^\circ\text{C}$ vaak als:
\[
E = E^\circ - \frac{0.05916}{n}\log_{10}Q
\]

**Wat is $Q$?**
Het reactiequotiënt: producten/reactanten met macht volgens stoichiometrie.
- Zuivere vaste stoffen en vloeistoffen komen **niet** in $Q$.

## 6) Concentratiecellen (klassiek examenvraagstuk)
Een **concentratiecel** heeft dezelfde halfreactie links en rechts, maar met **verschillende concentraties**.
Voorbeeldnotatie:
\[
\text{Ag(s)}\;|\;\text{Ag}^+(1.0\,\text{M})\;||\;\text{Ag}^+(0.1\,\text{M})\;|\;\text{Ag(s)}
\]

Belangrijk:
- De celspanning is typisch **klein**.
- Je berekent elke halfcelpotentiaal met **Nernst**, en dan:
\[
E_\text{cel}=E_\text{rechts}-E_\text{links}\quad(\text{of consistent met kathode–anode})
\]
- De kant met “hogere reductiepotentiaal” fungeert als **kathode**.

## 7) Commerciële voltacellen (wat je moet herkennen + kernreacties)
Je hoeft dit meestal niet hyper-diep af te leiden, maar je moet het **type, de idee en vaak de halfreacties kunnen plaatsen**.

### Loodaccumulator (auto, ~12 V totaal)
- Ongeveer **2 V per cel**, typisch 6 cellen in serie.
- Bij ontlading wordt **H₂SO₄ verbruikt** (dichtheid daalt).
- Oplaadbaar: reacties omkeerbaar via externe stroombron.
- Bij laden kan water-elektrolyse (H₂/O₂) optreden → veiligheidsaspect.

### Droge cellen (1.25–1.50 V)
- **Leclanché-element** (klassieke “zink-kool”-achtige cel).
- **Alkalische batterij**: langere levensduur (zinkanode corrodeert trager in basisch milieu).
- **Zilvercel**: gebruikt in kleine toestellen (uurwerken, pacemakers, hoorapparaten, …).

### Nikkel–cadmium (Ni–Cd, ~1.4 V)
- Oplaadbaar (producten blijven aan elektroden “kleven” volgens de cursuscontext).
- Toepassingen: boormachines, scheerapparaten, …

### Brandstofcel (H₂/O₂, $E_\text{cel}\approx 1.2$ V)
- Reagentia worden continu aangevoerd.
- Nettoreactie: $2\text{H}_2 + \text{O}_2 \rightarrow 2\text{H}_2\text{O}$
- Efficiëntie-idee: groot deel van theoretische $\Delta G$ → elektrische energie; nadelen: opslag en dure elektroden.

### Li-ion (conceptueel herkennen)
- **Anode: grafiet**
- **Kathode: LiCoO₂**
- **Li⁺-transport** tussen elektroden (intercalatie/de-intercalatie als idee).

## 8) Examenvallen & mini-checklist
- **Anode ≠ altijd positief**: in galvanische cel is anode typisch **negatief** (oxidatie), kathode positief (reductie).
- **$E^\circ$ nooit schalen met coëfficiënten**.
- **Tabelwaarden zijn reductiepotentialen**: als je een oxidatie gebruikt, keer je de reactie om maar je werkt nog steeds met reductie-$E^\circ$ in $E^\circ_\text{cel}=E^\circ_\text{kath}-E^\circ_\text{an}$.
- **Standaardcondities**: opgeloste species 1 M, gassen 1 bar, zuivere vaste stoffen/vloeistoffen activiteit 1.
- Bij “cel plat”: **$E_\text{cel}=0$ én $\Delta G=0$ én $Q=K$**.

## 9) Snelle “rekenflow” (wat je op je kladpapier wil)
1. Identificeer halfreacties + $E^\circ$ uit tabel.
2. Kies kathode = hoogste $E^\circ$.
3. $E^\circ_\text{cel}=E^\circ_\text{kath}-E^\circ_\text{an}$.
4. Balanceer reactie → bepaal $n$.
5. Indien niet-standaard: $E=E^\circ-\frac{RT}{nF}\ln Q$.
6. $\Delta G=-nFE$ en eventueel $K$ via $E^\circ_\text{cel}=\frac{RT}{nF}\ln K$.
