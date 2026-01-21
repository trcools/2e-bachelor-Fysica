# Quick Reference: Repository Bescherming

## Dagelijkse Workflow

### Voor je iets commit:

```bash
# 1. Check wat er veranderd is
git status

# 2. Bekijk de wijzigingen
git diff

# 3. Voeg alleen veilige bestanden toe
git add <bestand>

# 4. Commit met duidelijke message
git commit -m "Beschrijving van wijziging"
```

## Veelgebruikte Commando's

### Een bestand uit tracking verwijderen (maar lokaal behouden):

```bash
git rm --cached bestand.pdf
echo "bestand.pdf" >> .gitignore
git commit -m "Stop tracking bestand.pdf"
```

### Alle PDFs uit tracking verwijderen:

```bash
git rm --cached '*.pdf'
echo "*.pdf" >> .gitignore
git commit -m "Stop tracking all PDFs"
```

### Pre-commit hooks installeren:

```bash
pip install pre-commit
pre-commit install
```

### Pre-commit hooks testen:

```bash
pre-commit run --all-files
```

## Checklist voor Nieuwe Bestanden

Voordat je een nieuw bestand commit, vraag jezelf af:

- [ ] Is dit mijn eigen werk?
- [ ] Bevat het geen copyrighted materiaal van docenten?
- [ ] Bevat het geen persoonlijke/gevoelige informatie?
- [ ] Is het bestand niet te groot (<1MB)?
- [ ] Zou ik dit publiek willen delen?

Als je bij één van deze vragen twijfelt, commit het bestand **niet**.

## .gitignore Aanpassen

Om PDFs te blokkeren, verwijder de `#` in `.gitignore`:

```gitignore
# Van:
# *.pdf

# Naar:
*.pdf
```

Dan:
```bash
git add .gitignore
git commit -m "Block PDF files from being committed"
```

## Hulp Nodig?

Zie [SECURITY.md](SECURITY.md) voor uitgebreide richtlijnen.

