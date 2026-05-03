# Inleiding

Voor het examen relativiteit en electromagnetisme mogen we een A4 blad meenemen met daarop voor en achterkant aan formules, bijwijzen en definities.

# Doel

Het doel van dit project is om zo een A4 blad te maken omzetbaar naar  een pdf-bestand, waar er op de voor en achterkant alles op staat van definities, bewijzen en formules die te kennen zijn voor het examen.

# Te kennen leerstof

De te kennen leerstof staat [hier](https://github.com/trcools/2e-bachelor-Fysica/blob/sub-branch/2e-bachelor-semester-2/REM/REM_2024_chapter-summary.pdf). REM/REM_2024_chapter-summary.pdf

# Hulplijn

Dit is het [handboek](https://github.com/trcools/2e-bachelor-Fysica/blob/sub-branch/2e-bachelor-semester-2/REM/David%20J.%20Griffiths-Introduction%20to%20Electrodynamics-Addison-Wesley%20(2012).pdf) die we hebben gevolgd tijdens de les. Natuurlijk is niet alles van het handboek te kennen. REM/David J. Griffiths-Introduction to Electrodynamics-Addison-Wesley (2012).pdf

Dit is een [fomuleblaadje](https://github.com/trcools/2e-bachelor-Fysica/blob/sub-branch/2e-bachelor-semester-2/REM/RelEM-2026_formula_sheet.pdf) die we tijdens de oefeningen hebbeen gekregen, ik ben niet zeker of we die ook oop het examen krijigen. REM/RelEM-2026_formula_sheet.pdf

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
