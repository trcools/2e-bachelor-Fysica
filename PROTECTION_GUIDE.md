### Comitten

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

### Privé bestanden

__Gebruik de `_lokaal/` map__

Zet bestanden die je privé wilt houden in de `_lokaal/` map.
Git negeert deze map volledig — niets daarin wordt ooit gepusht naar GitHub.

```bash
# Kopieer een PDF naar de privémap:
cp ~/Downloads/cursusthermische.pdf _lokaal/Thermische-Fysica/

# Git zal dit bestand NIET zien:
git status   # _lokaal/*.pdf staat er NIET bij
```

Zie [`_lokaal/README.md`](_lokaal/README.md) voor een aanbevolen mapstructuur.

### Veelgebruikte Commando's

__Een bestand uit tracking verwijderen (maar lokaal bewaren):__

```bash
git rm --cached bestand.pdf
echo "bestand.pdf" >> .gitignore
git commit -m "Stop tracking bestand.pdf"
```

__Pre-commit hooks installeren:__

```bash
pip install pre-commit
pre-commit install
```

__Pre-commit hooks testen:__

```bash
pre-commit run --all-files
```

__Nog hulp Nodig?__

Zie [SECURITY.md](SECURITY.md) voor uitgebreide richtlijnen.
