# Projectieve Representaties — volledige uitleg (§3.7)

*Gebaseerd op Group Theory for Physicists (Verstraete & Vancraeynest-De Cuiper), met uitgebreide toelichting per stap.*

---

## 0. Waarom bestaat dit hoofdstuk? — de fysische motivatie

Dit is het stuk dat de meeste studenten missen, en zonder dit klopt niets van de rest.

In de kwantummechanica is de golffunctie $|\psi\rangle$ **geen** fysisch observabele grootheid. Wat je effectief meet is:

- de kansverdeling $P(x) = \bar\psi(x)\psi(x)$,
- matrixelementen $\langle\psi|\hat O|\psi\rangle$.

Beide zijn **invariant** onder $|\psi\rangle \to e^{i\varphi}|\psi\rangle$. Met andere woorden: $|\psi\rangle$ en $e^{i\varphi}|\psi\rangle$ zijn **fysisch dezelfde toestand**. Een golffunctie is dus eigenlijk geen vector, maar een **straal** (ray) — een vector op een globale fase na.

**Het probleem:** stel dat je systeem transformeert onder een symmetriegroep $G$:
$$|\psi\rangle \mapsto U(g)|\psi\rangle.$$

Omdat fase fysisch niet telt, *moet* het evenveel gelden als $|\psi\rangle$ getransformeerd wordt door $e^{i\varphi(g)}U(g)$ — exact dezelfde fysica. Ga je nu twee symmetrieën na elkaar toepassen, dan hoeft er **geen reden** te zijn waarom de fasefactoren netjes wegvallen. Het meest algemene wat je kan eisen is:

$$\boxed{U(g)U(h) = e^{i\omega(g,h)}\,U(gh)}$$

Dit is **geen** gewone representatie meer (daar zou $U(g)U(h)=U(gh)$ moeten gelden, zonder fase) — dit noemen we een **projectieve representatie**.

> **Kernidee in één zin:** een projectieve representatie is precies wat je krijgt als je representatietheorie toepast op stralen (toestanden op een fase na) in plaats van op vectoren.

---

## 1. De 2-cocycle voorwaarde — waarom $\omega$ niet willekeurig is

$\omega$ kan niet zomaar elke functie zijn: de $U(g)$'s moeten nog steeds **associatief** vermenigvuldigen, want matrixvermenigvuldiging is altijd associatief:

$$\big(U(g)U(h)\big)U(k) = U(g)\big(U(h)U(k)\big)$$

Werk beide kanten uit:

**Links:** $U(g)U(h)U(k) = e^{i\omega(g,h)}U(gh)U(k) = e^{i\omega(g,h)}e^{i\omega(gh,k)}U(ghk)$

**Rechts:** $U(g)U(h)U(k) = U(g)e^{i\omega(h,k)}U(hk) = e^{i\omega(h,k)}e^{i\omega(g,hk)}U(ghk)$

Gelijkstellen geeft de **2-cocycle conditie**:

$$\boxed{\omega(g,h) + \omega(gh,k) = \omega(g,hk) + \omega(h,k) \pmod{2\pi}}$$

**Waarom dit de naam "2-cocycle" draagt:** dit is exact dezelfde vergelijking die je tegenkwam bij groepscohomologie (Appendix 2.B) in de context van centrale extensies. Hier duikt diezelfde wiskundige structuur weer op, maar nu in een totaal andere toepassing — representatietheorie i.p.v. groepsuitbreidingen. **Dat hergebruik is precies waarom de cursus dit "group cohomology" stuk zo benadrukt.**

---

## 2. Gauge-vrijheid: niet elke $\omega$ is fysisch verschillend

Stel je herdefinieert $U(g) \to e^{i\varphi(g)}U(g)$ (een herkeuze van fase per groepselement — dit heet **gauge fixing**). Hoe verandert $\omega$?

$$e^{i\varphi(g)}U(g)\cdot e^{i\varphi(h)}U(h) = e^{i(\varphi(g)+\varphi(h))}e^{i\omega(g,h)}U(gh) = e^{i(\varphi(g)+\varphi(h)+\omega(g,h)-\varphi(gh))}\cdot e^{i\varphi(gh)}U(gh)$$

dus

$$\boxed{\omega(g,h) \to \omega(g,h) + \varphi(g)+\varphi(h)-\varphi(gh)}$$

De term $\varphi(g)+\varphi(h)-\varphi(gh)$ heet een **coboundary**. Twee cocycles die enkel van elkaar verschillen door een coboundary, geven **fysisch dezelfde** projectieve representatie — je hebt gewoon een andere fase-conventie gekozen.

**Gevolg:** de fysisch relevante informatie zit niet in $\omega$ zelf, maar in de **equivalentieklasse**

$$[\omega] \in H^2(G, U(1)) := \frac{\{\text{2-cocycles}\}}{\{\text{coboundaries}\}}$$

de **tweede cohomologiegroep**. Dit is exact dezelfde $H^2$ die je kent uit centrale extensies — het feit dat hij hier terugkomt is geen toeval, maar een diepe structurele connectie in de wiskunde.

> **Het idee om te onthouden:** $H^2(G,U(1))$ classificeert alle "écht verschillende" manieren waarop een groep projectief kan representeren. Is $[\omega]=0$ (triviale klasse), dan kan je via een gauge-keuze altijd terug naar een gewone lineaire representatie. Is $[\omega]\neq 0$, dan zit je vast aan een intrinsiek projectief effect dat je nooit kan wegtransformeren.

### Formele definitie

> **Definitie 3.14.** Een **projectieve representatie** van $G$ is een verzameling matrices $X_g$ die voldoen aan
> $$X_g X_h = e^{i\omega(g,h)} X_{gh}$$
> waarbij $\omega$ een representant is van een **niet-triviale** cohomologieklasse $[\omega] \in H^2(G,U(1))$.

(Niet-triviaal is belangrijk: als $[\omega]=0$ kan je het via een gauge-transformatie wegwerken tot een gewone representatie — dan is er niets "echt projectiefs" aan de hand.)

---

## 3. Een tweede motivatie: tensorproducten van representaties

Onafhankelijk van quantum-rays duikt dezelfde structuur op bij tensorproducten. Stel $U(g)\otimes\bar U(g)$ vormt een gewone representatie (dit gebeurt vaak, bv. bij dichtheidsmatrices $\rho \to U\rho U^\dagger$). Betekent dit dat $U(g)$ zelf een gewone representatie is?

**Niet noodzakelijk!** Het volstaat dat $U(g)$ een **projectieve** representatie is:

$$\big(U(g)\otimes\bar U(g)\big)\big(U(h)\otimes\bar U(h)\big) = e^{i\omega(g,h)}U(gh)\otimes e^{-i\omega(g,h)}\bar U(gh) = U(gh)\otimes\bar U(gh)$$

De fasefactoren $e^{i\omega}$ en $e^{-i\omega}$ **vallen exact tegen elkaar weg** in het tensorproduct, ongeacht of $\omega$ triviaal is of niet! Dit toont aan dat projectieve representaties geen exotische wiskundige curiositeit zijn, maar **natuurlijk verschijnen** telkens als je met fysische grootheden werkt die zelf al fase-ongevoelig zijn (zoals dichtheidsmatrices).

---

## 4. Uitgewerkt voorbeeld: $\mathbb{Z}_2 \oplus \mathbb{Z}_2$

Dit is **de kleinste groep met een niet-triviale 2-cocycle**: $H^2(\mathbb{Z}_2\oplus\mathbb{Z}_2, U(1)) \cong \mathbb{Z}_2$.

Representant van de niet-triviale klasse (tabel met rijen $g$, kolommen $h$):

$$e^{i\omega(g,h)} = \begin{pmatrix} & (0,0) & (1,0) & (0,1) & (1,1) \\ (0,0) & 1&1&1&1 \\ (1,0)&1&1&1&1\\ (0,1)&1&1&-1&-1\\(1,1)&1&1&-1&-1\end{pmatrix}$$

Een projectieve representatie die hierbij hoort, gebruikt de **Pauli-matrices**:

$$\rho(0,0)=\mathbb{1}, \quad \rho(0,1) = i\sigma_Y, \quad \rho(1,0)=-\sigma_X, \quad \rho(1,1)=\sigma_Z$$

**Controle (waarom dit werkt):**
$$\rho(0,1)\rho(1,0) = (i\sigma_Y)(-\sigma_X) = -i\sigma_Y\sigma_X = -i(-i\sigma_Z) = -\sigma_Z = -\rho(1,1)$$

Vergelijk met wat je verwacht: $\rho(0,1)\rho(1,0) \stackrel{!}{=} e^{i\omega((0,1),(1,0))}\rho(1,1) = (-1)\rho(1,1)$. ✓ Klopt exact met de tabel.

**Waarom dit voorbeeld zo leerrijk is:** het laat zien dat Pauli-matrices (die je toch al kent uit $\mathfrak{su}(2)$) van nature een **projectieve** representatie van een eenvoudige discrete groep vormen. Het anticommuteren van Pauli-matrices ($\sigma_X\sigma_Y = -\sigma_Y\sigma_X$) ís precies het niet-triviale teken in de cocycle.

---

## 5. Projectieve representaties construeren: de "projectieve reguliere representatie"

Herinner je de **gewone** reguliere representatie: $R(g) = \sum_h |gh\rangle\langle h|$ (permuteert de basisvectoren via linksvermenigvuldiging met $g$).

Gegeven een niet-triviale cocycle $\omega$, construeer je een projectieve versie door simpelweg een fase toe te voegen:

$$\boxed{R_\omega(g) = \sum_h e^{i\omega(g,h)}\,|gh\rangle\langle h|}$$

**Check dat dit werkt** (dit is de berekening die je op een examen zou moeten kunnen reproduceren):

$$R_\omega(g)R_\omega(h) = \sum_{x,y} e^{i\omega(g,x)}e^{i\omega(h,y)}|gx\rangle\langle x|hy\rangle\langle y|$$

De factor $\langle x|hy\rangle = \delta_{x,hy}$ dwingt $x=hy$:

$$= \sum_y e^{i\omega(g,hy)}e^{i\omega(h,y)}\,|ghy\rangle\langle y|$$

Pas nu de 2-cocycle conditie toe in de vorm $\omega(g,hy)+\omega(h,y) = \omega(gh,y)+\omega(g,h)$:

$$= e^{i\omega(g,h)}\sum_y e^{i\omega(gh,y)}|ghy\rangle\langle y| = e^{i\omega(g,h)} R_\omega(gh)$$

$\Rightarrow$ Dit **bewijst** dat $R_\omega$ inderdaad een geldige projectieve representatie is met exact cocycle $\omega$.

---

## 6. Normalisatie: $R_\omega(e) = \mathbb{1}$

Een technisch maar belangrijk detail: kan je altijd een gauge kiezen zodat het neutraal element gewoon de identiteit blijft (zoals je zou verwachten)?

**Definitie:** een cocycle is **genormaliseerd** als $\omega(g,e)=\omega(e,h)=0$ voor alle $g,h$.

**Propositie 3.12 — elke cocycle is equivalent met een genormaliseerde:**

*Stap 1:* Uit de cocycle-conditie met $h=g^{-1}$, $k=e$: $\omega(g,g^{-1})+\omega(e,e) = \omega(g,g^{-1})+\omega(g',e)$ voor een geschikte $g'$ — uitwerken toont aan dat $\omega(g,e)$ **onafhankelijk is van $g$**, dus $\omega(g,e) = \omega(e,e)$ voor alle $g$.

*Stap 2:* Gauge-transformeer met $\varphi(e) = \omega(e,e)$ en $\varphi(g)=0$ voor $g\neq e$. Dan $\omega(e,e) \to \omega(e,e) - \varphi(e) = 0$, en dus $\omega(g,e) \to 0$ voor alle $g$.

**Propositie 3.13 — je kan ook $\omega(g,g^{-1})=0$ forceren:**

Met een genormaliseerde cocycle geeft de cocycle-conditie $\omega(g,g^{-1}) = \omega(g^{-1},g)$. Kies vervolgens $\varphi(g) = \tfrac12\omega(g,g^{-1})$ (symmetrisch in $g \leftrightarrow g^{-1}$) om dit naar nul te transformeren.

**Gevolg (Corollary 3.13.1):** je kan *altijd* een gauge kiezen waarin $R_\omega(e) = \mathbb{1}$. Dit lijkt triviaal, maar vereist dus wel degelijk een bewijs — het is niet automatisch.

> **Waarom dit ertoe doet op een examen:** dit soort propositie test of je begrijpt dat "gauge-vrijheid" niet zomaar een woord is — je kan effectief concrete eigenschappen (zoals $R_\omega(e)=\mathbb 1$) afdwingen door een geschikte $\varphi$ te kiezen, en dat moet je kunnen *construeren*, niet enkel beweren.

---

## 7. Geen 1-dimensionale invariante deelruimtes — het belangrijkste structurele gevolg

Dit is misschien wel **de belangrijkste stelling** van de hele sectie, want ze verklaart waarom projectieve representaties fysisch zo'n groot verschil maken.

**Stelling:** als $\omega$ niet-triviaal is, dan heeft de projectieve representatie $\{X_g\}$ **geen** 1-dimensionale invariante deelruimte.

**Bewijs (contradictie):** stel er bestaat zo'n invariante deelruimte, opgespannen door $|v\rangle$, met
$$X_g|v\rangle = e^{i\alpha(g)}|v\rangle.$$

Pas nu twee keer toe:
$$X_gX_h|v\rangle = e^{i\alpha(g)}e^{i\alpha(h)}|v\rangle$$

maar ook, via de definitie van projectieve representatie:
$$X_gX_h|v\rangle = e^{i\omega(g,h)}X_{gh}|v\rangle = e^{i\omega(g,h)}e^{i\alpha(gh)}|v\rangle$$

Gelijkstellen:
$$\omega(g,h) = \alpha(g)+\alpha(h)-\alpha(gh)$$

Maar dit is precies de vorm van een **coboundary**! Dus $\omega$ zou triviaal zijn — **contradictie** met de aanname dat $\omega$ niet-triviaal is. $\blacksquare$

**Wat dit fysisch betekent:** voor een niet-triviale projectieve representatie kan je dus nooit "ontsnappen" naar een eendimensionale (= scalaire, fase-achtige) representatie. Het kleinst mogelijke irrep heeft minstens dimensie 2 — en dit is exact de wiskundige kern van waarom spin-$\tfrac12$ deeltjes (2-dimensionale projectieve irrep van $SO(3)$!) niet kunnen bestaan als 1-dimensionale ("scalaire") objecten.

**Aanvullend gevolg:** omdat 1-dimensionale invariante deelruimtes uitgesloten zijn, en omdat (analoog aan het gewone geval) invariante deelruimtes in de reguliere representatie altijd voorkomen met multipliciteit gelijk aan hun dimensie, volgt dat een groep van orde $\#G = \sum_i d_i^2$ enkel een niet-triviale cocycle kan hebben als minstens één $d_i \geq 2$.

---

## 8. De veralgemening van bekende stellingen naar het projectieve geval

Het mooie (en het examenrelevante!) van projectieve representaties is dat **bijna de hele machinerie van gewone representatietheorie overleeft**, zolang je twee irreps vergelijkt die tot **dezelfde** cohomologieklasse $[\omega]$ behoren:

| Stelling (gewoon geval) | Geldt nog steeds voor projectief, mits zelfde $[\omega]$? |
|---|---|
| Schur's lemma 1 (invariante deelruimte = hele ruimte, want unitair) | ✓ Ja |
| Schur's lemma 2 ($[X^{(\alpha)},B]=0 \,\forall g \Rightarrow B = c\mathbb{1}$) | ✓ Ja |
| Great Orthogonality Theorem | ✓ Ja, met een aangepast bewijs (zie hieronder) |
| Karakter-orthogonaliteit / klassenbegrip | ✓ Ja, met $\text{Tr}(X(x)) = \text{Tr}(X(y)X(x)X(y^{-1}))$ |
| Multipliciteit van irrep in reguliere representatie = $\dim$ | ✓ Ja |

**Schets van het bewijs van de veralgemeende orthogonaliteitsstelling:** start van dezelfde som als in het gewone bewijs, maar nu met extra fasefactoren $e^{i\omega(h,g)}$ uit beide irreps $\alpha$ en $\beta$. Omdat beide irreps **dezelfde** cohomologieklasse $[\omega]$ hebben, vallen de fasefactoren tegen elkaar weg:
$$e^{i\omega(h,g)}\cdot e^{-i\omega(h,g)} = 1$$
waarna de som exact terugvalt op de structuur van het gewone bewijs, en Schur's lemma (dat nog steeds geldt, zie boven) levert $B = c\cdot\mathbb{1}$ met $c=\text{Tr}(B)$ — vanwaar de gebruikelijke orthogonaliteitsrelaties volgen.

> **Cruciale voorwaarde om te onthouden:** dit werkt enkel als je **twee irreps van dezelfde cohomologieklasse** vergelijkt. Vergelijk je irreps van *verschillende* klassen $[\omega] \neq [\omega']$, dan vallen de fasefactoren niet tegen elkaar weg en breekt de stelling.

---

## 9. Waarom dit fysisch zo belangrijk is

De cursus noemt expliciet twee gevolgen, en het is goed om beide te kunnen reproduceren:

1. **Het bestaan van half-integer elektronspin.** Zoals je weet uit hoofdstuk 4: $SU(2)$ heeft zowel integer als half-integer spin-irreps, maar enkel de integer irreps zijn *echte* representaties van $SO(3)$ — de half-integer irreps zijn enkel **projectieve** representaties van $SO(3)$ (met cocycle gegeven door het teken $-1$ bij een rotatie over $2\pi$). Omdat golffuncties stralen zijn (een fase-vrijheid hebben), is dit perfect toelaatbaar: spin-$\tfrac12$ deeltjes (fermionen) "leven" precies in deze projectieve sector. Zonder de toelating van projectieve representaties zouden fermionen, en dus de stabiliteit van materie, niet kunnen bestaan.

2. **Symmetry-protected topologische fasen van materie** — in de moderne kwantummaterie-fysica classificeert men fasen die niet door lokale ordeparameters te onderscheiden zijn, maar wel door welke (niet-triviale) cohomologieklasse $[\omega] \in H^2(G,U(1))$ de symmetrie projectief representeert op de rand van het systeem.

---

## 10. Samenvattend stappenplan — "hoe pak ik een examenvraag over dit onderwerp aan?"

1. **Herken het probleem:** golffuncties zijn stralen (fase-onbepaald) $\Rightarrow$ symmetrierepresentaties moeten enkel sluiten op een fase na: $U(g)U(h)=e^{i\omega(g,h)}U(gh)$.
2. **Associativiteit afdwingen** $\Rightarrow$ leidt tot de 2-cocycle conditie $\omega(g,h)+\omega(gh,k)=\omega(g,hk)+\omega(h,k)$.
3. **Gauge-vrijheid herkennen:** $U(g)\to e^{i\varphi(g)}U(g)$ verschuift $\omega$ met een coboundary $\Rightarrow$ enkel $[\omega]\in H^2(G,U(1))$ is fysisch.
4. **Triviaal vs. niet-triviaal:** $[\omega]=0$ $\Rightarrow$ je kan terug naar een gewone representatie. $[\omega]\neq0$ $\Rightarrow$ écht nieuw fenomeen.
5. **Geen 1D irreps** als $[\omega]\neq 0$ — leid dit zelf af via de contradictie met de coboundary-vorm (sectie 7 hierboven, dit is een populaire examenvraag).
6. **Constructiemethode:** ken de projectieve reguliere representatie $R_\omega(g)=\sum_h e^{i\omega(g,h)}|gh\rangle\langle h|$ en kan de afleiding van $R_\omega(g)R_\omega(h)=e^{i\omega(g,h)}R_\omega(gh)$ reproduceren via de cocycle-conditie.
7. **Connecteer met hoofdstuk 4:** spin-$\tfrac12$ = 2-dimensionale projectieve irrep van $SO(3)$, mogelijk gemaakt doordat $SU(2)$ de dekgroep is (zie eerdere samenvatting over topologie).

---

## 11. Veelgemaakte denkfouten (en hoe je ze vermijdt)

- **"Een projectieve representatie is gewoon een representatie met een extra fase."** Niet correct als verklaring — de fase $\omega(g,h)$ hangt af van *beide* groepselementen $g$ én $h$, niet enkel van het product. Dat is precies wat het tot een 2-cocycle maakt i.p.v. een simpele globale herschaling.
- **"Je kan altijd $\omega$ wegtransformeren."** Enkel als $[\omega]=0$ in $H^2(G,U(1))$. De hele essentie van de sectie is dat dit voor sommige groepen (zoals $\mathbb{Z}_2\oplus\mathbb{Z}_2$, of $SO(3)$) niet kan.
- **"Schur's lemma geldt niet meer voor projectieve representaties."** Het geldt wél nog — zolang je binnen dezelfde cohomologieklasse blijft vergelijken. Dat is een vaak vergeten randvoorwaarde.
- **Verwarring tussen $\omega$ en $[\omega]$.** $\omega$ is een specifieke functie (een gauge-keuze); $[\omega]$ is de equivalentieklasse modulo coboundaries — enkel die laatste is fysisch/groep-theoretisch fundamenteel.
