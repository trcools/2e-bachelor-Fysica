# PDF naar Markdown Converter

Dit document beschrijft de algemene workflow om PDF's naar markdown om te zetten met formules als tekst en figuren als afbeeldingen.

## Script-locatie

`pdf-naar-md-bestand/pdf_naar_md.py`

## Wat doet het?

Het script zet PDF's om naar markdown met:

- Kleine formuleblokken naar `$$ ... $$`
- Figuren en grafieken als `.png` in een assets-map
- Per PDF een `.md` bestand met dezelfde bestandsnaam

## Installatie (eenmalig)

Navigeer naar de map `2e-bachelor-semester-2` in de repo en maak een venv aan (of gebruik een bestaande):

```bash
cd <repo-root>/2e-bachelor-semester-2
python -m venv .venv
source .venv/bin/activate        # macOS/Linux
# .venv\Scripts\activate         # Windows
pip install pymupdf rapidocr_onnxruntime
```

Je mag een bestaande venv gebruiken. Belangrijk is enkel dat de twee packages (`pymupdf` en `rapidocr_onnxruntime`) daarin geïnstalleerd zijn.

## Gebruik:

### Basis

```bash
cd <repo-root>/2e-bachelor-semester-2
source .venv/bin/activate        # macOS/Linux
python pdf-naar-md-bestand/pdf_naar_md.py --input-dir REM/Samenvatting-REM
```

### Andere map

```bash
python pdf-naar-md-bestand/pdf_naar_md.py --input-dir PAD/NAAR/JE/PDFMAP
```

### Met aparte outputmap

```bash
python pdf-naar-md-bestand/pdf_naar_md.py --input-dir PAD/NAAR/PDFS --output-dir PAD/NAAR/OUTPUT
```

### Recursief (ook submappen)

```bash
python pdf-naar-md-bestand/pdf_naar_md.py --input-dir PAD/NAAR/PDFS --recursive
```

## Handige opties

- `--formula-max-width 240`
- `--formula-max-height 90`
- `--render-scale 4`
- `--assets-suffix _assets`

Voorbeeld met fijnregeling:

```bash
python pdf-naar-md-bestand/pdf_naar_md.py \
	--input-dir REM/Samenvatting-REM \
	--formula-max-width 260 \
	--formula-max-height 110 \
	--render-scale 4
```

## Git ignore tip

Als je gegenereerde output niet naar GitHub wil sturen, voeg dit toe in de relevante `.gitignore`:

```gitignore
*.md
*_assets/
```

Of specifieker per map (voorbeeld):

```gitignore
/Samenvatting-REM/*.md
/Samenvatting-REM/*_assets/
```

## Troubleshooting

**Script is traag**
OCR kost tijd op grote PDF's. Dat is normaal.

**Formules missen of zijn slordig**
Verhoog `--render-scale` of pas de formulebreedte/-hoogte aan.

**Afbeeldingen worden als formule gezien**
Verlaag `--formula-max-width` en `--formula-max-height`.
