# Beveiligingsrichtlijnen / Security Guidelines

## Nederlands

### Bescherming van Cursusmateriaal

Deze repository bevat studiemateriaal voor de 2e bachelor Fysica. Om de repository goed te beschermen en academische integriteit te waarborgen, volg deze richtlijnen:

#### Wat NIET te committen:

1. **Persoonlijke/Gevoelige informatie**

   - Studentnummers
   - Persoonlijke notities met gevoelige informatie
   - API keys of wachtwoorden
2. **Auteursrechtelijk beschermd materiaal**

   - Professorslides (hoorcolleges)
   - Cursusboeken en syllabi van docenten
   - Commercieel gepubliceerde leerboeken
   - Examenopgaven van het lopende academiejaar
3. **Grote binaire bestanden**

   - Video's van colleges
   - Grote datasets (gebruik externe hosting)

#### Wat WEL te committen:

1. **Eigen werk**

   - Zelfgemaakte samenvattingen
   - Eigen code en scripts
   - Eigen notities en uitwerkingen
2. **Open source materiaal**

   - Jupyter notebooks met eigen analyses
   - Python scripts voor berekeningen
   - Documentatie en README bestanden

---

### PDFs en documenten privé bewaren

Wil je PDFs (slides, syllabi, leerboeken) bewaren zonder ze publiek te maken?
Je hebt drie opties:

#### Optie 1 — `_lokaal/` map (aanbevolen voor lokaal gebruik)

Er is een speciale map `_lokaal/` aangemaakt in deze repository. Alles wat je daarin plaatst wordt **nooit** door Git bijgehouden.

```bash
# Kopieer een PDF naar de lokale privémap:
cp ~/Downloads/cursusthermische.pdf _lokaal/Thermische-Fysica/

# Git zal dit bestand negeren:
git status  # _lokaal/*.pdf verschijnt NIET in de output
```

Zie [`_lokaal/README.md`](_lokaal/README.md) voor een aanbevolen mappenstructuur.

#### Optie 2 — Maak de repository privé

Ga naar **GitHub → Settings → Danger Zone → Change visibility → Make private**.
Dan kun je alles in de repository bewaren en zijn je bestanden alleen zichtbaar voor jou (en uitgenodigde medewerkers).

#### Optie 3 — Aparte privérepository

Maak een tweede, privé GitHub-repository aan (bijv. `2e-bachelor-Fysica-lokaal`) en bewaar de privébestanden daarin. De publieke repository blijft dan publiek voor je eigen werk.

---

### .gitignore Configuratie

De `.gitignore` blokkeert automatisch:

- `*.pdf`, `*.docx`, `*.pptx` — documentbestanden
- `_lokaal/**` — alles in de privémap

Om één specifieke PDF toch te committen (bijv. een eigen samenvatting):

```bash
git add -f mijn_samenvatting.pdf
```

### Pre-commit Hooks

De repository gebruikt pre-commit hooks die automatisch controleren op:

- Grote bestanden (>500KB wordt gewaarschuwd)
- Private keys
- Merge conflicten
- Python syntax errors

Installeer ze met:

```bash
pip install pre-commit
pre-commit install
```

### Wat te doen als je per ongeluk gevoelig materiaal hebt gecommit:

1. **Verwijder het bestand uit de tracking (maar bewaar het lokaal):**

   ```bash
   git rm --cached bestand.pdf
   ```
2. **Voeg het toe aan .gitignore:**

   ```bash
   echo "bestand.pdf" >> .gitignore
   ```
3. **Commit de wijziging:**

   ```bash
   git add .gitignore
   git commit -m "Remove sensitive file from tracking"
   ```
4. **Voor reeds gepushte gevoelige data:** Neem contact op met de repository eigenaar of gebruik `git filter-branch` of BFG Repo-Cleaner.
