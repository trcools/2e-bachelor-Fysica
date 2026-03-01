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

## Privé bestanden bewaren (slides, PDFs, leerboeken)

### Gebruik de `_lokaal/` map

Zet bestanden die je privé wilt houden in de `_lokaal/` map.
Git negeert deze map volledig — niets daarin wordt ooit gepusht naar GitHub.

```bash
# Kopieer een PDF naar de privémap:
cp ~/Downloads/cursusthermische.pdf _lokaal/Thermische-Fysica/

# Git zal dit bestand NIET zien:
git status   # _lokaal/*.pdf staat er NIET bij
```

Zie [`_lokaal/README.md`](_lokaal/README.md) voor een aanbevolen mapstructuur.

## Veelgebruikte Commando's

### Een bestand uit tracking verwijderen (maar lokaal bewaren):

```bash
git rm --cached bestand.pdf
echo "bestand.pdf" >> .gitignore
git commit -m "Stop tracking bestand.pdf"
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
Zet het in `_lokaal/` als je het toch wilt bewaren.

## Hulp Nodig?

Zie [SECURITY.md](SECURITY.md) voor uitgebreide richtlijnen.
