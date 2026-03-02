# _lokaal — Privé bestanden (alleen lokaal)

Deze map is bedoeld voor bestanden die je **alleen lokaal** wilt bewaren — ze worden **nooit** door Git bijgehouden of naar GitHub gepusht.

## Wat kun je hier bewaren?

- Professorslides (PDFs van hoorcolleges)
- Cursusboeken en syllabi
- Examens en oefeningen die je van Ufora hebt gedownload
- Elk ander bestand dat je privé wilt houden

## Hoe werkt het?

De `.gitignore` bevat de regel `_lokaal/**`, wat betekent:
- Alles wat je in deze map plaatst blijft **alleen op jouw computer**
- Git zal deze bestanden nooit tracken of uploaden naar GitHub
- Alleen dit `README.md` bestand is zichtbaar op GitHub

## Aanbevolen mappenstructuur

Je kunt submappen aanmaken die overeenkomen met de vakken in de repository:

```
_lokaal/
├── README.md                      ← dit bestand (zichtbaar op GitHub)
├── Biophysics/
│   ├── slides/                    ← weekelijkse slides van de prof
│   └── Biophysics_2025_2026-1.pdf ← examenopgave
├── Thermische-Fysica/
│   ├── cursusthermische.pdf
│   └── slides/
├── Sterrenstelsels/
│   ├── Cursus_sterrenstelsels_2025_2026.pdf
│   └── slides/
├── Chemie-2/
│   └── Powerpoints/               ← hoorcolleges HC1-HC22
├── Groepen-en-Representaties/
│   └── group-representation-theory.pdf
└── REM/
    └── Griffiths-Electrodynamics.pdf
```

## Wil je toch cloud-back-up?

Als je ook online toegang wil tot deze privébestanden, heb je twee opties:

### Optie 1: Maak deze repository privé

Ga naar **GitHub → Settings → Danger Zone → Change visibility → Make private**.
Dan kun je alles in de repository bewaren zonder dat het publiek zichtbaar is.
Je kunt daarna via GitHub nog steeds je eigen bestanden bekijken.

### Optie 2: Maak een aparte privérepository

Maak een tweede, privérepository aan op GitHub (bijv. `2e-bachelor-Fysica-lokaal`).
Zet daarin alleen de privébestanden. De publieke repository blijft dan publiek voor je eigen werk.

```bash
# Maak een nieuwe lokale repo van de _lokaal map:
cd _lokaal
git init
git remote add origin git@github.com:<JOUW-GITHUB-GEBRUIKERSNAAM>/2e-bachelor-Fysica-lokaal.git
git add .
git commit -m "Voeg privébestanden toe"
git push -u origin main
```
