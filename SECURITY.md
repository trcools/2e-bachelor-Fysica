# Beveiligingsrichtlijnen / Security Guidelines

## Nederlands

### Bescherming van Cursusmateriaal

Deze repository bevat studiemateriaal voor de 2e bachelor Fysica. Om de repository goed te beschermen en academische integriteit te waarborgen, volg deze richtlijnen:

#### Wat NIET te committen:

1. **Auteursrechtelijk beschermde materialen**
   - PDF's van colleges en powerpoints van docenten
   - Examens en oude tentamens
   - Gescande boekpagina's of handouts

2. **Persoonlijke/Gevoelige informatie**
   - Studentnummers
   - Persoonlijke notities met gevoelige informatie
   - API keys of wachtwoorden

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
5. ❌ Commit geen copyrighted materiaal van docenten
6. ❌ Deel geen examens of antwoordmodellen
7. ❌ Commit geen persoonlijke gegevens

---

## English

### Protecting Course Materials

This repository contains study materials for 2nd year Physics Bachelor. To properly protect the repository and maintain academic integrity, follow these guidelines:

#### What NOT to commit:

1. **Copyrighted materials**
   - PDF slides from professors
   - Exams and past papers
   - Scanned book pages or handouts

2. **Personal/Sensitive information**
   - Student IDs
   - Personal notes with sensitive info
   - API keys or passwords

3. **Large binary files**
   - Lecture videos
   - Large datasets (use external hosting)

#### What TO commit:

1. **Your own work**
   - Self-made summaries
   - Your own code and scripts
   - Your own notes and solutions

2. **Open source materials**
   - Jupyter notebooks with your analyses
   - Python scripts for calculations
   - Documentation and README files

### .gitignore Configuration

The `.gitignore` file contains commented lines that you can uncomment to block PDFs and other documents:

```gitignore
# Remove the # to block PDFs:
# *.pdf
# *.docx
# *.pptx
```

### Pre-commit Hooks

The repository uses pre-commit hooks that automatically check for:
- Large files (>500KB warning)
- Private keys
- Merge conflicts
- Python syntax errors

Install them with:
```bash
pip install pre-commit
pre-commit install
```

### What to do if you accidentally committed sensitive material:

1. **Remove the file from staging area:**
   ```bash
   git rm --cached file.pdf
   ```

2. **Add it to .gitignore:**
   ```bash
   echo "file.pdf" >> .gitignore
   ```

3. **Commit the change:**
   ```bash
   git add .gitignore
   git commit -m "Remove sensitive file from tracking"
   ```

4. **For already pushed sensitive data:** Contact the repository owner or use `git filter-branch` or BFG Repo-Cleaner (advanced).

### Best Practices

1. ✅ Only share your own summaries and code
2. ✅ Use clear commit messages
3. ✅ Always check `git status` before committing
4. ✅ Review your changes with `git diff`
5. ❌ Don't commit copyrighted materials from professors
6. ❌ Don't share exams or answer keys
7. ❌ Don't commit personal information

---

## Reporting Security Issues

If you discover a security issue (accidentally committed passwords, API keys, etc.), please:
1. Do NOT create a public issue
2. Contact the repository owner directly
3. Delete sensitive data from your local copy

