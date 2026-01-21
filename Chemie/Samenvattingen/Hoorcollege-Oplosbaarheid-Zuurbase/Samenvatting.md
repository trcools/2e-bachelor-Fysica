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
\mathrm{M_aX_b(s) \rightleftharpoons a\,M^{n+}(aq) + b\,X^{m-}(aq)}
$$

Dan (in verdunde oplossingen, met activiteiten $\approx$ concentraties):
$$
K_s = [\mathrm{M^{n+}}]^a[\mathrm{X^{m-}}]^b
$$

### Molaire oplosbaarheid $S$ in zuiver water
Laat $S = x$ mol/L zout oplossen. Dan:
- $[\mathrm{M^{n+}}] = a x$
- $[\mathrm{X^{m-}}] = b x$

Dus:
$$
K_s = (a x)^a (b x)^b
\quad\Rightarrow\quad
x = S
$$

**Mini-cheatsheet (komt rechtstreeks uit de slides/oefeningen):**

- $\text{MX(s)} \rightleftharpoons \text{M}^+ + \text{X}^-$  
  $K_s = x^2 \Rightarrow S = \sqrt{K_s}$

- $\text{M}_2 \text{X(s)} \rightleftharpoons 2\text{M}^+ + \text{X}^{2-}$  
  $K_s = (2x)^2 x = 4x^3 \Rightarrow S = \left(\dfrac{K_s}{4}\right)^{1/3}$

- $\text{M}_2 \text{X}_3(s) \rightleftharpoons 2\text{M}^{3+} + 3\text{X}^{2-}$  
  $K_s = (2x)^2 (3x)^3 = 108x^5 \Rightarrow S = \left(\dfrac{K_s}{108}\right)^{1/5}$


**Waarom dit belangrijk is:** in de slides staat expliciet een oefening waar zouten met verschillende stoichiometrie een $K_s$ krijgen en je *moet* rangschikken op oplosbaarheid. Dat kan enkel via $S$.

---

## 3) Gemeenschappelijk-ion-effect (dé klassieker)
Als er al een ion aanwezig is dat in het evenwicht voorkomt, dan verschuift het evenwicht **naar links** ⇒ **oplosbaarheid daalt**.

### Voorbeeldtype uit de slides: $\mathrm{Ag_2CrO_4}$
$$
\mathrm{Ag_2CrO_4(s) \rightleftharpoons 2Ag^+ + CrO_4^{2-}}
$$
$$
K_s = [\mathrm{Ag^+}]^2[\mathrm{CrO_4^{2-}}]
$$

Stel: je hebt al $0{,}10\ \mathrm{mol/L}$ $\mathrm{Ag^+}$ in oplossing, en extra oplosbaarheid is $x$:
- $[\mathrm{Ag^+}] = 0{,}10 + 2x$
- $[\mathrm{CrO_4^{2-}}] = x$

Dan:
$$
K_s = (0{,}10+2x)^2(x)
$$

**De exametruc die de prof toont:**
Als $0{,}10 \gg 2x$, dan:
$$
(0{,}10+2x)\approx 0{,}10
\quad\Rightarrow\quad
K_s \approx (0{,}10)^2 x
\Rightarrow x\approx \frac{K_s}{(0{,}10)^2}
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
- Vergeten dat bij $\mathrm{M_2X}$: $[\mathrm{M}]=2S$ en niet $S$.
- Benadering maken (“$0{,}10+2x \approx 0{,}10$”) en **niet** checken.
- Denken dat “groter $K_s$ altijd groter $S$” is — *soms*, maar niet universeel.

---

## 6) Wat je “paraat” wil hebben voor het examen
- Definitie: $K_s$ als evenwichtsconstante van oplossen.
- Omzetting $K_s \leftrightarrow S$ via stoichiometrie (voor 1:1, 2:1, 2:3 moet je dit snel kunnen).
- Gemeenschappelijk-ion-effect kunnen opstellen én benaderen.
- Het mantra: **“$K_s$ is niet oplosbaarheid; $S$ is oplosbaarheid.”**
