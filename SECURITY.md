# Beveiligingsrichtlijnen / Security Guidelines

## Nederlands

### Bescherming van Cursusmateriaal

Deze repository bevat studiemateriaal voor de 2e bachelor Fysica. Om de repository goed te beschermen en academische integriteit te waarborgen, volg deze richtlijnen:

#### Wat NIET te committen:

1. **Persoonlijke/Gevoelige informatie**
   - Studentnummers
   - Persoonlijke notities met gevoelige informatie
   - API keys of wachtwoorden

2. **Grote binaire bestanden**
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

### .gitignore Configuratie

De `.gitignore` file bevat commentaar bij regels die je kunt uitcommentariëren om PDFs en andere documenten te blokkeren:

```gitignore
# Verwijder de # om PDFs te blokkeren:
# *.pdf
# *.docx
# *.pptx
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

1. **Verwijder het bestand uit de staging area:**
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

4. **Voor reeds gepushte gevoelige data:** Neem contact op met de repository eigenaar of gebruik `git filter-branch` of BFG Repo-Cleaner (gevorderd).

### Best Practices

1. ✅ Deel alleen je eigen samenvattingen en code
2. ✅ Gebruik duidelijke commit messages
3. ✅ Check altijd `git status` voor je commit
4. ✅ Review je wijzigingen met `git diff`
5. ❌ Commit geen persoonlijke gegevens


