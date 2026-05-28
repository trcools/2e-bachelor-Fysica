# latex

Uitvoer en hulpbestanden van LaTeX‑builds (PDF, .aux, .log, etc.).

- Doel: bewaar gecompileerde PDF's en LaTeX intermediate files hier.
- Voorbeeld: `my_cheatsheet.pdf`, `.aux`, `.log` bestanden.
- Tip: voeg `.gitignore` regels toe om grote of tijdelijke bestanden uit te sluiten.

### LaTeX build artifacts are safe to delete once you're done—they're regeneratable byproducts.

Use latexmk -c to clean build artifacts without deleting the PDF
Use latexmk -C to clean everything (including PDF)
