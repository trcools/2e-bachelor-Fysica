# Sterrenstelsels — volledig herschreven cursus (vlotte versie)

Deze versie herschrijft de volledige inhoud in duidelijke taal, met minder herhaling en met focus op fysisch inzicht.
Je krijgt per hoofdstuk:

- wat je conceptueel moet begrijpen,
- welke formules echt belangrijk zijn,
- hoe je de theorie aan observaties koppelt.

---

## Hoe je dit document gebruikt

1. Lees eerst de uitleg van een hoofdstuk zonder formules te memoriseren.
2. Kijk daarna pas naar de formulelijst en check of je elke term fysisch kan uitleggen.
3. Sluit af met één korte testvraag: “kan ik dit zonder notities uitleggen?”.

Doel: **snappen > blokken**.

### Universele symbolen (komen overal terug)

- $c$: lichtsnelheid
- $h$: constante van Planck
- $k$: constante van Boltzmann
- $G$: gravitatieconstante
- $m_p$, $m_H$: protonmassa en waterstofmassa
- $M_\odot$: zonsmassa

Als een formule niet klopt in een oefening, check dan eerst:

1) **eenheden**, 2) **aannames** (isotropie, evenwicht, dun/dik medium), 3) **grootorde**.

---

## Hoofdstuk 1 — Historiek van sterrenstelsels

### Waarover gaat dit?

Dit hoofdstuk geeft de wetenschappelijke context: hoe astronomen van “diffuse nevels” naar het moderne beeld van externe melkwegstelsels gingen.
De centrale vraag was lang: zijn die nevels objecten binnen de Melkweg of zelfstandige stelsels ver daarbuiten?

De echte waarde van dit hoofdstuk is niet rekenwerk, maar het begrijpen van de methodologie:

- welke observaties beschikbaar waren,
- waarom eerdere modellen fout of onvolledig waren,
- hoe betere metingen (afstanden, helderheid, structuur) de discussie beslisten.

### Wat je moet onthouden

- Extragalactische sterrenkunde is historisch jong.
- Interpretatie van waarnemingen evolueerde mee met instrumentele vooruitgang.
- Het “grote debat” was essentieel om het schaalbeeld van het heelal te corrigeren.

### Formules

- Geen kernformules in dit hoofdstuk; het is vooral conceptueel.

---

## Hoofdstuk 2 — Morfologie van melkwegstelsels en observaties

### Waarover gaat dit?

Dit hoofdstuk leert je de basistaal van observationele sterrenkunde:

- wat je meet (flux, spectra, beelden),
- hoe je dat omzet naar fysica (lichtkracht, kleur, structuur),
- en hoe je stelsels classificeert.

De grote lijn is: **meten -> kalibreren -> interpreteren**.

### Belangrijkste ideeën

#### 2.1 Waarnemen van melkwegstelsels

In praktijk is “licht” in de cursus breder dan zichtbaar licht: optisch, NIR, UV, radio, röntgen...
Elke golflengte vertelt iets anders over dezelfde galaxie (oude sterren, jong gas, stof, warme/hete componenten).

Fotometrie geeft totale flux per band, spectroscopie geeft verdeling over golflengte, en IFU-data koppelt beide ruimtelijk.

#### 2.2 Classificatie

De Hubble-sequentie is een beschrijvend schema: elliptisch, lenticulair, spiraal, balkspiraal, irregulair.
Belangrijk: dit is geen volledige “evolutievolgorde”, maar een morfologische ordening.

#### 2.3 Projectie en deprojectie

Wat je ziet is een projectie op de hemel. Hellingshoek beïnvloedt vorm, oppervlaktehelderheid en afgeleide parameters.

#### 2.4 Modellen

Lichtprofielen vatten structuur compact samen:

- de Vaucouleurs voor klassieke bulges/elliptischen,
- exponentiële wet voor schijven,
- Sérsic als algemene vorm.

Bulge-schijf decompositie is essentieel om fysische componenten apart te analyseren.

### Kernformules

$$
\lambda = \frac{c}{\nu}
$$

$$
F_\lambda = \frac{L_\lambda}{4\pi d^2}
$$

$$
m_X = -2.5\log_{10}\left(\frac{F_X}{F_{X,0}}\right)
$$

$$
X-Y = m_X - m_Y
$$

$$
j = 10\left(1-\frac{b}{a}\right)
$$

$$
\mu = -2.5\log_{10}I + C
$$

$$
I(R)=I_e\exp\left[-7.669\left(\left(\frac{R}{R_e}\right)^{1/4}-1\right)\right]
$$

$$
I(R)=I_0\exp\left[-b_m\left(\frac{R}{R_e}\right)^{1/m}\right]
$$

$$
I(R)\propto e^{-R/h}
$$

### Vanwaar komen deze formules?

- $\lambda=c/\nu$ komt uit de golfeigenschap snelheid = golflengte × frequentie.
- $F_\lambda=L_\lambda/(4\pi d^2)$ is de inverse-kwadratenwet: dezelfde totale energie wordt over een boloppervlak verdeeld.
- Magnitudes zijn logaritmisch omdat helderheidsperceptie en historische fotometrie log-schaal gebruiken.
- de Vaucouleurs- en Sérsic-profielen zijn **empirische fits** op waargenomen helderheidsprofielen; ze zijn dus data-gedreven modellen.
- Het exponentiële schijfprofiel volgt uit de waarneming dat schijflicht in veel spiralen ongeveer exponentieel afneemt met straal.

### Symbolen in dit hoofdstuk

- $F_\lambda$: gemeten flux
- $L_\lambda$: intrinsieke lichtkracht
- $d$: afstand
- $m_X$: magnitude in band $X$
- $R_e$: effectieve straal (helft van het licht binnen $R_e$)
- $h$: schijfschaallengte
- $m$ (in Sérsic): concentratie-index

---

## Hoofdstuk 3 — Interstellair gas

### Waarover gaat dit?

Gas is de ruwe bouwstof voor stervorming en bepaalt sterk de evolutie van een sterrenstelsel.
Je leert de belangrijkste gasfasen en hoe observaties in verschillende golflengtes die fasen zichtbaar maken.

### Belangrijkste ideeën

#### 3.1–3.2 Gascomponenten

- Neutraal atomair gas (H I)
- Moleculair gas (stervorming)
- Warm/heet geïoniseerd gas

Elke fase heeft een eigen temperatuur-, dichtheids- en emissieregime.

#### 3.3 Waarnemingen

- Radio: o.a. 21 cm, cruciaal voor H I-verdeling en gasmassa.
- Röntgen: traceert heet diffuse component.

#### 3.4 Energiebudget

Interstellair gas wordt beïnvloed door stralingsvelden, magnetische velden en kosmische straling.
Dat bepaalt ionisatie, druksteun, koeling en dynamica.

### Kernformules

$$
\frac{M_{\mathrm{HI}}}{M_\odot}=2.356\times10^5\,D^2\int S_\nu(v_{\mathrm{rad}})\,dv_{\mathrm{rad}}
$$

$$
E=h\nu=\frac{hc}{\lambda}
$$

### Vanwaar komen deze formules?

- De H I-massaformule komt uit radiotransfer + kalibratie van de 21 cm-lijnflux naar aantal H I-atomen.
- De factor $2.356\times10^5$ bevat eenhedenconversies (Jy, km/s, Mpc) en fysieke constanten.
- $E=h\nu$ is de quantumrelatie voor fotonen en verbindt direct waargenomen frequentie met energie-inhoud.

### Symbolen in dit hoofdstuk

- $D$: afstand in Mpc
- $S_\nu(v_{\mathrm{rad}})$: fluxdichtheid als functie van radiële snelheid
- $v_{\mathrm{rad}}$: radiële snelheid
- $M_{\mathrm{HI}}$: neutrale waterstofmassa

---

## Hoofdstuk 4 — Interstellair stof

### Waarover gaat dit?

Stof verandert observaties systematisch: objecten lijken zwakker en roder.
Wie stof niet corrigeert, maakt snel fouten in afstanden, leeftijden, metalliciteit of SFR.

### Belangrijkste ideeën

#### 4.1–4.2 Observaties van stof

Stof zie je indirect via extinctie in UV/optisch en direct via thermische emissie in IR/sub-mm.

#### 4.2.1 Extinctie en reddening

Totale verzwakking is niet gelijk aan kleurverandering; beide worden apart gekwantificeerd.

#### 4.3 Oorsprong

Stof ontstaat in evoluerende sterren en supernova-omgevingen, en wordt verwerkt in het ISM.

### Kernformules

$$
A_\lambda=-2.5\log_{10}\left(\frac{F_\lambda}{F^0_\lambda}\right)
$$

$$
I_\lambda=I^0_\lambda e^{-\tau_\lambda}
$$

$$
A_\lambda\approx1.086\,\tau_\lambda
$$

$$
R_V=\frac{A_V}{A_B-A_V}=\frac{A_V}{E(B-V)}
$$

$$
E(B-V)=(B-V)_{\mathrm{obs}}-(B-V)_{\mathrm{int}}
$$

$$
\frac{N_H}{E(B-V)}\approx5.8\times10^{21}\;\mathrm{cm^{-2}\,mag^{-1}}
$$

$$
\frac{A_V}{N_H}\approx5.3\times10^{-22}\;\mathrm{mag\,cm^2\,H^{-1}}
$$

### Vanwaar komen deze formules?

- $A_\lambda$-definitie volgt uit de magnitude-definitie toegepast op verzwakte versus intrinsieke flux.
- $I_\lambda=I^0_\lambda e^{-\tau_\lambda}$ komt uit exponentiële absorptie in radiatieve overdracht.
- De factor $1.086$ in $A_\lambda\approx1.086\tau_\lambda$ is puur de omzetting van natuurlijke log naar logaritme basis 10.
- Relaties met $N_H$ zijn empirisch bepaald voor typische zichtlijnen en kunnen variëren per omgeving.

### Symbolen in dit hoofdstuk

- $A_\lambda$: extinctie (in magnitudes)
- $\tau_\lambda$: optische diepte
- $E(B-V)$: kleur-exces (reddening)
- $R_V$: totale/selectieve extinctie
- $N_H$: waterstofkolomdichtheid

---

## Hoofdstuk 5 — Stralingsprocessen

### Waarover gaat dit?

Dit hoofdstuk is de fysische basis voor spectra:

- continuümprocessen,
- lijnvorming,
- lijnverbreding,
- diagnostische grootheden zoals equivalente breedte.

### Belangrijkste ideeën

#### 5.1 Historiek en quantumkoppeling

Foto-elektrisch en Compton-effect tonen dat licht deeltjeskarakter heeft met impuls en energie.

#### 5.2 Continuümstraling

De zwarte straler is de referentie. Veel astrofysische bronnen benader je hiermee om temperatuur/energie-output te schatten.

#### 5.3–5.4 Lijnstraling en moleculaire overgangen

Energieniveaus bepalen mogelijke lijnen; botsingen en stralingsvelden bepalen populaties en lijnsterkte.

#### 5.5 Lijnprofielen

Natuurlijke, Doppler- en botsingsverbreding bepalen de vorm van lijnen. Die vorm bevat fysische info over temperatuur, dichtheid en turbulentie.

### Kernformules

$$
E=h\nu=\frac{hc}{\lambda}=pc
$$

$$
\Delta\lambda=\frac{h}{m_ec}(1-\cos\theta)
$$

$$
B_\nu(T)=\frac{2h\nu^3}{c^2}\frac{1}{e^{h\nu/kT}-1}
$$

$$
L=4\pi R^2\sigma_{\mathrm{SB}}T^4
$$

$$
T_{\mathrm{eff}}=\left(\frac{L}{4\pi\sigma_{\mathrm{SB}}R^2}\right)^{1/4}
$$

$$
\lambda_{\max}T=b\approx2.898\times10^{-3}\;\mathrm{m\,K}
$$

$$
\frac{1}{\lambda}=R_H\left(\frac{1}{m^2}-\frac{1}{n^2}\right)
$$

$$
n_{\mathrm{crit}}\sim\frac{\sum A_{ul}}{\sum q_{ul}}
$$

$$
P_L(\nu)=\frac{1}{\pi}\frac{\gamma}{(\nu-\nu_0)^2+\gamma^2}
$$

$$
P_D(\nu)=\frac{1}{\sqrt{2\pi}\,\sigma\nu_0}\exp\left[-\frac{(\nu-\nu_0)^2}{2\sigma^2\nu_0^2}\right]
$$

$$
\sigma=\sqrt{\frac{kT}{mc^2}}
$$

$$
W_\lambda=\int\left(1-\frac{F_\lambda}{F_0}\right)d\lambda
$$

### Vanwaar komen deze formules?

- Planck, Stefan-Boltzmann en Wien komen uit thermodynamica + quantumfysica voor straling in evenwicht.
- Comptonverschuiving volgt uit behoud van energie en impuls bij foton-elektronbotsing.
- Rydberg-formule beschrijft energieniveaus in waterstof (overgangen tussen kwantumniveaus).
- $n_{\mathrm{crit}}$ komt uit competitie tussen radiatieve de-excitatie en botsingen.
- Lijnprofielen (Lorentz, Gauss) zijn modelvormen voor verschillende verbredingsmechanismen.

### Symbolen in dit hoofdstuk

- $B_\nu(T)$: Planckfunctie
- $\sigma_{\mathrm{SB}}$: Stefan-Boltzmannconstante
- $R_H$: Rydbergconstante
- $A_{ul}$: Einsteincoëfficiënt (spontane emissie)
- $q_{ul}$: botsingscoëfficiënt
- $W_\lambda$: equivalente breedte

---

## Hoofdstuk 6 — Stervorming

### Waarover gaat dit?

Stervorming is een fysisch evenwichtsprobleem: wanneer kan een gaswolk haar eigen zwaartekracht niet meer weerstaan?
Je volgt het traject van instabiele wolk tot jonge sterren en de observabele signatures daarvan.

### Belangrijkste ideeën

#### 6.1 Instabiliteit en collapse

Viriaalbalans en Jeans-criterium bepalen of collapse mogelijk is.

#### 6.2–6.4 Protosterren en H II-gebieden

Massieve jonge sterren ioniseren hun omgeving: Strömgren-sferen geven een schaal voor deze ionisatiezones.

#### 6.5 Tracers van recente stervorming

H$\alpha$, UV en IR meten recente SFR, elk met eigen gevoeligheid voor stof en tijdsschaal.

#### 6.6–6.7 SFH en schaalrelaties

Stervormingsgeschiedenis wordt vaak met parametrische modellen beschreven.
De Kennicutt-Schmidt-relatie linkt gasdichtheid aan SFR-dichtheid.

#### 6.8 Toepassing op de Melkweg

Viriale massabepalingen en wolkschaalrelaties verbinden theorie met meetbare grootheden.

### Kernformules

viriaal theorema:

$$
2K+U=0
$$

$$
U\approx-\frac{A\,GM^2}{R}\quad(A\sim3/5)
$$

$$
K=\frac{3}{2}NkT=\frac{3}{2}\frac{MkT}{\mu m_H}
$$

Jeansmassa:

$$
M>M_J
$$

$$
M_J\approx\left(\frac{5kT}{G\mu m_H}\right)^{3/2}\left(\frac{3}{4\pi\rho_0}\right)^{1/2}
$$

**Ionisatietempo**:

$$
R_S=\left(\frac{3Q_0}{4\pi n_H^2\alpha_B(T)}\right)^{1/3}
$$

Gasmassa binnen de Strömgrensfeer (gas nagenoeg compleet geïoniseerd):

$$
M_{\mathrm{HII}}=\frac{4\pi}{3}R_S^3m_pn_H
$$

Energieetoestanden bepaald door de Boltzmann vergelijking:

$$
\frac{n_u}{n_l}=\frac{g_u}{g_l}e^{-h\nu_{ul}/kT}
$$

Star formation rate (SFR) voor de $H_\alpha$ lijn gecorrigeerd voor extinctie door interstellair stof:

$$
\mathrm{SFR}=5.37\times10^{-35}\left(\frac{L_{\mathrm{H}\alpha}}{\mathrm{W}}\right)M_\odot\,\mathrm{yr^{-1}}
$$

Een typische schatter voor de stervormingsgraad:

$$
\mathrm{SFR}=4.47\times10^{-37}\left(\frac{L_{\mathrm{FUV}}}{\mathrm{erg\,s^{-1}}}\right)M_\odot\,\mathrm{yr^{-1}}
$$

met $L_{\text{FUV}}$ de lichtkracht in de GALEX FUV band. Deze is handig om galaxieën te selecteren met abnormaal hoge UV-flux.

Stervorming met obscured (verduisterde) en unobscured (onverduisterde) term:

$$
\mathrm{SFR}=5.37\times10^{-35}\left[\frac{L_{\mathrm{H}\alpha}}{\mathrm{W}}+0.020\frac{L_{24\mu m}}{\mathrm{W}}\right]M_\odot\,\mathrm{yr^{-1}}
$$

Populaire vorm van parametrisatie "het tau model":

$$
\mathrm{SFR}(t)\propto e^{-(t-T_0)/\tau}
$$

met $T_0$ een bepaald tijdstip waarna de SFR exponentieeel afneemt.

"Delayed tau model":

$$
\mathrm{SFR}(t)\propto (t-T_0)e^{-(t-T_0)/\tau}
$$

(Laat een toename van de SFR over tijd toe.)

De trend tussen de stervormingsgraad en de dichtheid van de totale koude gasmassa, waarbij beiden vaak in eenheden van oppervlaktedichtheid worden uitgedrukt:

$$
\Sigma_{\mathrm{SFR}}=A\,\Sigma_{\mathrm{gas}}^n
$$

Wordt ook wel de Kennicutt-Schmidt relatie of stervormingswet genoemd.

$$
M_{\mathrm{vir}}\approx\frac{5R\sigma_r^2}{G},\qquad \alpha_{\mathrm{vir}}=\frac{5\sigma_r^2R}{GM}
$$

$\alpha_{vir} \equiv \frac{2K}{U}$ is de viriaalparameter, de verhouding tussen de kinetische en potentiele energie uit het viriaaltheorema. als deze $\approx$ 1, dan is er evenwicht tussen de druk en zwaartekracht.

### Vanwaar komen deze formules?

- Jeans- en viriaalrelaties komen uit energie-evenwicht tussen thermische beweging en zwaartekracht.
- $R_S$ (Strömgren) volgt uit ionisatie-evenwicht: aantal ionisaties per seconde = aantal recombinaties per seconde.
- SFR-kalibraties (H$\alpha$, FUV, IR) zijn empirisch geijkt op populatiesynthesemodellen + waarnemingen.
- Kennicutt-Schmidt is een geobserveerde schaalrelatie tussen gasdichtheid en stervormingsactiviteit.

### Symbolen in dit hoofdstuk

- $M_J$: Jeansmassa
- $Q_0$: aantal ioniserende fotonen per seconde
- $\alpha_B(T)$: case-B recombinatiecoëfficiënt
- $\Sigma_{\mathrm{SFR}}$, $\Sigma_{\mathrm{gas}}$: oppervlaktesnelheden/dichtheden
- $\sigma_r$: (radiële) snelheidsdispersie

---

## Hoofdstuk 7 — Levenscyclus van sterren

### Waarover gaat dit?

De massa van een ster bepaalt bijna alles: levensduur, energiebudget, en eindtoestand.
Voor sterrenstelsels is dit hoofdstuk cruciaal omdat stellair evolutieproces gas verrijkt en terugkoppeling veroorzaakt.

### Belangrijkste ideeën

- Lage-massasterren evolueren traag en eindigen rustiger (witte dwergen).
- Hoge-massasterren evolueren snel en eindigen explosief (supernovae, compacte resten).
- Supernovae en nucleosynthese bepalen de chemische evolutie van galaxieën.

### Kernformules

Ideale gaswet:

$$
P=nkT=\frac{\rho kT}{\mu m_H}
$$

De coulomb-barriere: - de elektrostatische energie die moet overwonnen worden zodat 2 atoomkernen dicht genoeg zijn om aan kernfussie te doen -

$$
U_{\mathrm{coul}}=\frac{kZ_1Z_2e^2}{r}
$$

Energie dat vrijkomt bij het instorten van een sferisch symmetrische ster met constante dichtheid (bij benadering)

$$
E\approx-\frac{3}{10}\frac{GM^2}{R}
$$

### Vanwaar komen deze formules?

- De ideale gaswet koppelt microscopische thermische energie aan macroscopische druk in sterinterieur.
- Coulombbarrière beschrijft hoeveel elektrostatische afstoting moet overwonnen worden voor fusie.
- De gravitatie-energieformule is een orde-groottebenadering uit het viriaalidee voor zelf-graviterende bollen.

### Symbolen in dit hoofdstuk

- $\rho$: massadichtheid
- $\mu$: gemiddelde moleculaire massa
- $Z_1, Z_2$: atoomnummers van reagerende kernen

---

## Hoofdstuk 8 — Stellaire populaties

### Waarover gaat dit?

Je observeert geïntegreerd licht, niet individuele sterren (behalve lokaal).
Daarom modelleer je een stelsel als superpositie van populaties met verschillende leeftijden, massa’s en metalliciteiten.

### Belangrijkste ideeën

#### 8.1 Populatiesynthese

Single Stellar Population (SSP) is de bouwsteen: één leeftijd, één metalliciteit, één IMF.

#### 8.2 Populaties in verschillende stelseltypes

Elliptischen worden vaak gedomineerd door oude sterren; spiraalstelsels tonen sterkere menging van oud en jong.

#### 8.3 Melkweg

Populaties I/II/III en metalliciteitspatronen helpen de vormingsgeschiedenis reconstrueren.

### Kernformules


De geïntegreerde fluxformule is letterlijk de som van bijdragen van alle stertypes, geschaald met afstand:

$$
F_\lambda=\frac{1}{4\pi d^2}\sum_i N_iL_{\lambda,i}
$$

Initiele massafunctie (IMF):

IMF-vormen (Salpeter, Kroupa) zijn empirische/semiepirische beschrijvingen van stervormingsuitkomsten.

$$
\xi(M)\propto M^{-2.35}\quad\text{(Salpeter)}
$$

$$
\xi(M)\propto
\begin{cases}
M^{-1.2}, & M\lesssim0.5M_\odot\\
M^{-2.2}, & 0.5M_\odot\lesssim M\lesssim1.0M_\odot\\
M^{-4.5}, & M\gtrsim1.0M_\odot
\end{cases}
\quad\text{(Kroupa, zoals in de cursusnotatie)}
$$

$[\mathrm{Fe}/\mathrm{H}]$ is een logaritmische maat omdat abundantieverhoudingen over grote dynamische bereiken lopen:

$$
[\mathrm{Fe}/\mathrm{H}]=\log_{10}\left[\frac{(N_{\mathrm{Fe}}/N_H)_\star}{(N_{\mathrm{Fe}}/N_H)_\odot}\right]
$$

### Symbolen in dit hoofdstuk

- $N_i$: aantal sterren van type $i$
- $L_{\lambda,i}$: spectraalluminositeit van type $i$
- $\xi(M)$: initiële massafunctie
- $[\mathrm{Fe}/\mathrm{H}]$: metalliciteitsindicator t.o.v. de zon

---

## Hoofdstuk 9 — Stellaire kinematica en dynamica

### Waarover gaat dit?

Hier gebruik je bewegingen om de massa- en potentiaalstructuur van stelsels te bepalen.
De lijn is: spectrum -> snelheid -> snelheidsverdeling -> dynamisch model.

### Belangrijkste ideeën

#### 9.1 LOSVD

Uit Dopplerverschuivingen haal je line-of-sight snelheden. Hun distributie (LOSVD) draagt info over rotatie, dispersie en anisotropie.

#### 9.2 Fundamenteel vlak

Relaties zoals Faber-Jackson en het fundamenteel vlak koppelen observabele grootheden aan dynamische toestand en structuur.

#### 9.3 Ruimtelijk opgeloste kinematica

IFU-data laat lokale dynamica zien: anisotropie, opwarming van schijven, kinematisch ontkoppelde componenten.

#### 9.4 Melkwegkinematica

Oort-constanten en lokale benaderingen zijn praktisch om de nabij-zonneomgeving te karakteriseren.

### Kernformules

Dopplerformules volgen uit golfkinematica in de limiet $v\ll c$.

$$
\Delta\lambda\approx\lambda_0\frac{v_{\mathrm{los}}}{c}
$$

$$
u\equiv c\ln\lambda,\qquad \Delta u\approx v_{\mathrm{los}}
$$

De convolutie $G=F\otimes S$ zegt: waargenomen spectrum = templatespectrum “verbreed” door snelheidsverdeling.

$$
G(u)\propto F(u)\otimes S(u)
$$

Momenten ($\bar v_{\mathrm{los}}$, $\sigma_{\mathrm{los}}$) komen uit kansrekening op de LOSVD.

$$
\bar v_{\mathrm{los}}=\int F(v_{\mathrm{los}})v_{\mathrm{los}}\,dv_{\mathrm{los}}
$$

$$
\sigma^2_{\mathrm{los}}=\int F(v_{\mathrm{los}})(v_{\mathrm{los}}-\bar v_{\mathrm{los}})^2\,dv_{\mathrm{los}}
$$

Gaussische verdelingsfunctie:

$$
F(v_{\mathrm{los}})=\frac{1}{\sqrt{2\pi}\sigma_{\mathrm{los}}}\exp\left[-\frac{(v_{\mathrm{los}}-\bar v_{\mathrm{los}})^2}{2\sigma^2_{\mathrm{los}}}\right]
$$

$$
v_r\approx A\,d\,\sin(2l)
$$

$$
A\approx14.8\;\mathrm{km\,s^{-1}\,kpc^{-1}},\quad B\approx-12.4\;\mathrm{km\,s^{-1}\,kpc^{-1}}
$$

### Vanwaar komen deze formules?

- 
- Oort-relaties zijn lokale linearisaties van de galactische rotatie rond de positie van de zon.

### Symbolen in dit hoofdstuk

- $v_{\mathrm{los}}$: snelheid langs de gezichtslijn
- $\bar v_{\mathrm{los}}$: gemiddelde lijn-of-sight snelheid
- $\sigma_{\mathrm{los}}$: lijn-of-sight snelheidsdispersie
- $A,B$: Oort-constanten
- $l$: galactische lengte

---

## Hoofdstuk 10 — Donkere materie en alternatieven

### Waarover gaat dit?

Dit hoofdstuk vertrekt van één observatief feit: rotatiekrommen blijven vaak vlak waar je een Kepler-afname verwacht.
Dat dwingt je naar twee interpretaties:

- extra massa (donkere materie),
- of aangepaste zwaartekracht in zwakke-versnellingsregimes.

### Belangrijkste ideeën

#### 10.1 Spiraalstelsels

Rotatiekrommen geven directe massa-informatie als functie van straal.

#### 10.2 Elliptische stelsels

Massa-informatie komt uit stellaire kinematica, planetaire nevels/bolhopen en X-ray halo’s.

#### 10.3 Aard van donkere materie

Baryonische kandidaten, WIMPs en alternatieve frameworks worden vergeleken op observabele consequenties.

#### 10.3.3 Alternatieven (MOND)

In het lage-versnellingsregime verandert de effectieve dynamische wet zodat vlakke rotatiekrommen natuurlijk ontstaan.

### Kernformules

$$
v_{\mathrm{circ}}(R)=\sqrt{\frac{GM(R)}{R}},\qquad M(R)=\frac{Rv_{\mathrm{circ}}^2}{G}
$$

$$
v_{\mathrm{circ}}\propto R^{-1/2}\quad\text{(Kepler verwacht zonder extra massa)}
$$

$$
L\propto v_{\max}^{\alpha},\quad \alpha\approx4
$$

$$
M_b=M_\star+M_{\mathrm{gas}}\approx M_\star+1.4M_{\mathrm{HI}}
$$

$$
\frac{dP}{dr}=-\frac{GM_{\mathrm{tot}}(<r)\rho_{\mathrm{gas}}(r)}{r^2}
$$

$$
P=n_{\mathrm{gas}}kT=\frac{\rho_{\mathrm{gas}}}{\mu m_p}kT
$$

$$
\mathbf{F}=m\,\mu\!\left(\frac{a}{a_0}\right)\mathbf{a},\qquad \frac{GM(R)}{R^2}=\mu\!\left(\frac{a}{a_0}\right)a
$$

$$
a\approx\frac{\sqrt{GM(R)a_0}}{R},\qquad v_{\mathrm{circ}}\approx\left(GM(R)a_0\right)^{1/4}
$$

$$
L\propto v_{\mathrm{circ}}^4
$$

### Vanwaar komen deze formules?

- $v_{\mathrm{circ}}=\sqrt{GM(R)/R}$ volgt uit centripetale balans in een sferisch potentiaal.
- Kepler-afname $v\propto R^{-1/2}$ geldt wanneer buiten een straal de ingesloten massa quasi constant is.
- Tully-Fisher is een empirische relatie tussen rotatiesnelheid en lichtkracht/baryonmassa.
- Hydrostatisch evenwicht in X-ray gas komt uit krachtenbalans tussen drukgradiënt en zwaartekracht.
- MOND-formules wijzigen effectief de dynamica in het lage-versnellingsregime $a\ll a_0$.

### Symbolen in dit hoofdstuk

- $v_{\mathrm{circ}}$: cirkelsnelheid
- $M(R)$: ingesloten massa binnen straal $R$
- $M_b$: baryonische massa
- $\rho_{\mathrm{gas}}$: gasdichtheid
- $a_0$: MOND-versnellingsschaal

---

## Rode draad door de volledige cursus

Als je alles reduceert tot één keten, dan is het dit:

1. **Waarnemen**: flux, spectra, beelden, snelheden.
2. **Corrigeren**: afstand, extinctie, projectie, instrumenteffecten.
3. **Afleiden**: massa, temperatuur, dichtheid, SFR, metalliciteit, kinematica.
4. **Interpreteren**: vorming en evolutie van sterrenstelsels.

Als je bij elke oefening expliciet deze 4 stappen benoemt, ga je sneller én met minder stress leren.

---

## Snelle examenstrategie (zonder overkill)

- Begin met hoofdstuk 2 + 5: dit zijn je “taal” en “fysica” van observaties.
- Voeg hoofdstuk 4 toe vóór je veel toepassingen doet (stofcorrecties vermijden fouten).
- Pak dan hoofdstuk 6 + 10: grote fysische thema’s met veel examenvragen.
- Eindig met 8 + 9: interpretatie van spectra en dynamica op hoger niveau.

Vuistregel: als je een formule niet in woorden kan uitleggen, ken je ze nog niet echt.
