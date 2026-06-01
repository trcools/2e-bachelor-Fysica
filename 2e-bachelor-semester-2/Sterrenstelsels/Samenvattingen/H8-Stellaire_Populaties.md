# H8 - Stellaire populaties

## Kernidee

In extragalactische sterrenkunde observeer je meestal niet individuele sterren, maar het geintegreerde licht van een volledig stelsel. Dat licht is een mengsel van sterren met verschillende leeftijden, massa's en metalliciteiten. Om een melkwegstelsel fysisch te begrijpen, moet je dus de totale waarneming ontbinden in stellare populaties en hun evolutie modelleren.

Het belangrijkste inzicht is dat het spectrum van een stelsel niet enkel afhangt van hoeveel sterren er zijn, maar vooral van welke sterren domineren in lichtbijdrage. Jonge, massieve sterren zijn extreem helder maar kortlevend; oude populaties bevatten veel meer sterren, maar zijn per ster doorgaans veel zwakker. Daardoor wordt het geintegreerde licht sterk gewogen naar de meest lichtkrachtige componenten.

## Populatiesynthese

De standaardbouwsteen is een Single Stellar Population (SSP): een ideale populatie waarin alle sterren tegelijk gevormd zijn, met dezelfde initiele metalliciteit en dezelfde initiele massafunctie (IMF). Zo'n SSP is geen letterlijk waargenomen object, maar een theoretisch basismodel waarmee men complexere stelsels kan opbouwen.

Voor een samengesteld stelsel wordt het geintegreerde spectrum geschreven als som van bijdragen van verschillende stertypes:

$$
F_\lambda = \frac{1}{4\pi d^2} \sum_i N_i L_{\lambda,i}
$$

Hierin is $N_i$ het aantal sterren van type $i$, $L_{\lambda,i}$ de spectraalluminositeit per stertype en $d$ de afstand. In continue vorm kan men dit uitbreiden tot een integraal over massa, leeftijd en metalliciteit:

$$
F_\lambda(t) = \frac{1}{4\pi d^2} \int \xi(M)\,L_\lambda(M,t,Z)\,dM
$$

Deze relatie maakt duidelijk waarom populatiesynthese essentieel is: uit het geobserveerde spectrum probeer je terug te redeneren welke mix van SSP's het stelsel vormt.

## De initiele massafunctie

De IMF beschrijft hoeveel sterren er gevormd worden per massabereik. Ze is een van de belangrijkste invoerparameters in populatiesynthese, omdat ze bepaalt hoe sterk lage- en hoge-massasterren bijdragen aan het licht, de massa en de chemische evolutie.

Een klassieke empirische benadering is de wet van Salpeter:

$$
\xi(M) \propto M^{-2.35}
$$

Voor moderne toepassingen wordt vaak een gebroken machtwet gebruikt, zoals in de Kroupa-achtige notatie van de cursus:

$$
\xi(M) \propto
\begin{cases}
M^{-1.2}, & M \lesssim 0.5\,M_\odot \\
M^{-2.2}, & 0.5\,M_\odot \lesssim M \lesssim 1.0\,M_\odot \\
M^{-4.5}, & M \gtrsim 1.0\,M_\odot
\end{cases}
$$

De precieze vorm van de IMF is niet uit eerste principes af te leiden, maar wordt empirisch bepaald. Toch heeft die keuze grote gevolgen: een topzware IMF levert veel meer UV-licht, meer supernovae en snellere verrijking van het interstellair medium.

## Leeftijd, kleur en spectra

De leeftijd van een populatie bepaalt in sterke mate haar kleur en spectrale eigenschappen. Jonge populaties bevatten nog O- en B-sterren en zijn dus blauw en UV-rijk. Oude populaties missen die massieve sterren en worden roder, met een spectrum dat sterker wordt gedragen door koelere sterren op de hoofdreeks, subreuzen en reuzen.

Daarom wordt de waargenomen kleur van een stelsel gebruikt als ruwe indicator voor zijn sterpopulatie. Een blauw stelsel wijst vaak op recente of lopende stervorming, terwijl een rood stelsel meestal wordt gedomineerd door oude sterren.

Belangrijk is wel dat kleur niet uniek is: ouderdom, metalliciteit en stofverduistering kunnen grotendeels dezelfde waargenomen kleur produceren. Dat is het klassieke leeftijd-metalliciteit-stof-degeneratieprobleem.

## Metalliciteit en chemische evolutie

Metalliciteit geeft aan hoeveel elementen zwaarder dan helium aanwezig zijn in sterren of gas. In de sterrenkunde wordt die informatie vaak samengevat via:

$$
[\mathrm{Fe}/\mathrm{H}] = \log_{10}\left[\frac{(N_{\mathrm{Fe}}/N_H)_\star}{(N_{\mathrm{Fe}}/N_H)_\odot}\right]
$$

Een positieve waarde betekent dat het object metaalrijker is dan de Zon, een negatieve waarde dat het metaalarmer is. Metalliciteit is belangrijk omdat ze informatie draagt over eerdere generaties sterren: zware elementen worden immers opgebouwd in sterren en verspreid via sterwinden en supernovae.

In de vroegste populaties verwacht men dus zeer lage metalliciteiten. Naarmate stervorming en terugkoppeling doorgaan, stijgt de metaalinhoud van het gas en dus van latere generaties sterren.

## Populaties in verschillende stelseltypes

Verschillende soorten stelsels tonen verschillende populatiesamenstellingen. Elliptische stelsels worden typisch gedomineerd door oude sterren, met weinig of geen recente stervorming. Hun licht wordt dus vooral geleverd door een oude, relatief uniforme populatie.

Spiraalstelsels bevatten een mengsel van oude en jonge sterren. Hun schijf herbergt doorgaans lopende stervorming, terwijl de bolcomponent en de buitenste delen eerder oude populaties bevatten.

Voor de Melkweg is het klassieke onderscheid tussen populaties I, II en III historisch en fysisch belangrijk:

- Populatie I: metaalrijke, jongere sterren in de schijf en armen.
- Populatie II: metaalarme, oudere sterren in de halo en bolhopen.
- Populatie III: hypothetische eerste generatie sterren, vrijwel zonder metalen.

Die indeling is nuttig als eerste ordeningsprincipe, maar in werkelijkheid bestaat er een continuum van leeftijden en metalliciteiten. De vormingsgeschiedenis van een stelsel reconstrueer je dus uit patronen in leeftijd, metalliciteit, ruimtelijke verdeling en kinematica.

## Waarom dit hoofdstuk belangrijk is

Stellaire populaties vormen de brug tussen waarneming en geschiedenis. Uit een geintegreerd spectrum probeer je af te leiden hoe een stelsel is opgebouwd, wanneer zijn sterren gevormd zijn, hoe snel het gas verrijkt raakte en welke evolutiekanalen dominant waren. Dat maakt populatiesynthese een fundamenteel instrument in de extragalactische sterrenkunde.

## Kernformules

$$
F_\lambda = \frac{1}{4\pi d^2} \sum_i N_i L_{\lambda,i}
$$

$$
\xi(M) \propto M^{-2.35}
$$

$$
\xi(M) \propto
\begin{cases}
M^{-1.2}, & M \lesssim 0.5\,M_\odot \\
M^{-2.2}, & 0.5\,M_\odot \lesssim M \lesssim 1.0\,M_\odot \\
M^{-4.5}, & M \gtrsim 1.0\,M_\odot
\end{cases}
$$

$$
[\mathrm{Fe}/\mathrm{H}] = \log_{10}\left[\frac{(N_{\mathrm{Fe}}/N_H)_\star}{(N_{\mathrm{Fe}}/N_H)_\odot}\right]
$$

## Samenvatting in een zin

Hoofdstuk 8 leert dat een melkwegstelsel moet worden begrepen als een samengestelde verzameling van stellare populaties, waarvan leeftijd, IMF en metalliciteit samen het geobserveerde licht en de evolutie van het stelsel bepalen.
