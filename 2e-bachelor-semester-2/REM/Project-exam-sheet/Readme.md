# Inleiding

Voor het examen relativiteit en elektromagnetisme mogen we een A4-blad meenemen met daarop op de voor- en achterkant formules, bewijzen en definities.

# Doel

Het doel van dit project is om zo'n A4-blad te maken, omzetbaar naar een pdf-bestand, waarop op de voor- en achterkant alle definities, bewijzen en formules staan die gekend moeten zijn voor het examen.

# Te kennen leerstof

De te kennen leerstof staat niet als PDF in deze repo, omdat PDF-bestanden hier niet standaard mee geversioneerd worden.
Plaats het document daarom lokaal, bijvoorbeeld als:

- `_lokaal/REM_2024_chapter-summary.pdf`

Als je dit bestand nog niet hebt, haal het dan via het officiële cursusmateriaal, Toledo of via de docent.

# Hulplijn

Het handboek dat in de les gevolgd werd is *Introduction to Electrodynamics* van David J. Griffiths.
Gebruik hiervoor een legale bron, bijvoorbeeld de universiteitsbibliotheek, het officiële cursusplatform of een persoonlijke legale kopie.

Het formuleblaadje staat ook niet in deze repo. Plaats het lokaal, bijvoorbeeld als:

- `_lokaal/RelEM-2026_formula_sheet.pdf`

Als je dit document niet hebt, vraag het na bij medestudenten, de docent of controleer het cursusplatform.

# Samenvatting

Er staan [hier](https://github.com/paa-bachelor/semester-2/tree/19fc942d88c8d307774d6a29201e2fe955be93c9/relativity-and-electromagnetism/summary) samenvattingen van de lessen van een mede collega die gebruikt kunnen worden.

# Uitvoering in deze map

Er is nu een printklare exam-sheet gemaakt in dit project:

- `rem_exam_sheet.tex` (A4, 2 pagina's: voor- en achterkant)

Inhoud:

- Kernformules elektrostatische en magnetostatische velden
- Maxwellvergelijkingen en golfrelaties
- Randvoorwaarden, potentiaalformulering en Poynting-vector
- Speciale relativiteit: Lorentztransformatie, energie-impuls, veldtransformaties

Compileren naar PDF (vanuit deze map):

```bash
pdflatex rem_exam_sheet.tex
pdflatex rem_exam_sheet.tex
```

Dan krijg je `rem_exam_sheet.pdf` die je dubbelzijdig kan afdrukken.
