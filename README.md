# 2e Bachelor Fysica - Studiemateriaal

Deze repository bevat studiemateriaal, samenvattingen, code en opdrachten voor de bachelor Fysica.

## 📚 Repository Overzicht

Deze repository bevat **academisch studiemateriaal** dat theorie combineert met computationele methodes:
- **Markdown samenvattingen** (Nederlands of Engels) voor theoretische vakken
- **Jupyter notebooks** voor numerieke demonstraties en data-analyse
- **Python code** voor wetenschappelijke berekeningen
- **Opdrachten en examens** met modulaire code-structuren

**Dit is GEEN software project** — het is educatief materiaal georganiseerd per semester en vak.

## 📂 Structuur

### 2e Bachelor - Semester 1
- **`Biophysics/`** - Biofysica opdrachten met Python code
- **`Kwantum/`** - Kwantummechanica leerstof en samenvattingen
- **`Py4sci/`** - Python for Scientific Computing (notebooks, samenvattingen)
- **`VFR/`** - Vector en Functieruimten

### 2e Bachelor - Semester 2
- **`Groepen-en-Representaties/`** - Groepentheorie en representaties
- **`Materiaal-Fysica/`** - Materiaalkunde
- **`REM/`** - Relativiteit en Elektromagnetisme
- **`Sterrenstelsels/`** - Astrofysica: sterrenstelsels en kosmologie
- **`Thermische-Fysica/`** - Statistische mechanica en thermodynamica
- **`project-exp-2/`** - Experimenteel project semester 2

### Andere Folders
- **`1e-bachelor-semester-1/`** - Materiaal uit 1e bachelor (o.a. lineaire algebra)
- **`Cursussen:boeken/`** - Cursussen en studieboeken
- **`GIT-commands/`** - Git referentiemateriaal

## 🔧 Code Conventies

### Py4sci: Self-Contained Teaching Notebooks
Notebooks in `Py4sci/` zijn **zelfstandige leerdocumenten**:
```python
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize, integrate

def newton_method(f, fp, x0, niter):
    """Illustrative Newton's method implementation."""
    # Code prioritizes clarity over optimization
    pass
```

**Matplotlib Pattern** (voorkomt duplicate figures):
```python
plt.close("unique_name")
fig, ax = plt.subplots(num="unique_name", figsize=(8, 5))
# plotting code...
```

### Biophysics: Modular Separation Pattern
Voor opdrachten wordt code gescheiden in modules:
```
Assignment-X/
├── imports.py      # Alle library imports
├── parameters.py   # Parameters via get_parameters()
├── data.py         # Analyse functies (integratie, eigenwaarden, ...)
├── figures.py      # Alle plot functies
└── notebook.ipynb  # Schone notebook die bovenstaande importeert
```

## ⚠️ Beveiliging & Privacy

**Lees de [SECURITY.md](SECURITY.md) voor uitgebreide beveiligingsrichtlijnen.**

### Quick Reference:
- ✅ **Commit**: samenvattingen, code, notebooks, opdrachten
- ✅ **Gebruik pre-commit hooks** voor automatische controles
- ❌ **Geen persoonlijke informatie** of gevoelige data committen

### Pre-commit Hooks Installeren

```bash
pip install pre-commit
pre-commit install
```

Dit controleert automatisch op:
- Grote bestanden (>5MB)
- Private keys of credentials
- Merge conflicts
- YAML/JSON syntax

## 🚀 Setup

### Python Dependencies

Voor numerieke notebooks (vooral Py4sci):

```bash
pip install numpy scipy matplotlib jupyter ipywidgets
```

Of als er een `requirements.txt` aanwezig is:
```bash
pip install -r requirements.txt
```

### Jupyter Notebook

Start Jupyter vanuit de repository root:
```bash
jupyter notebook
```

## 📘 Handige Links

- [GITHUB_INFO.md](GITHUB_INFO.md) - Uitleg over GitHub timestamps (rode/oranje kleuren)
- [PROTECTION_GUIDE.md](PROTECTION_GUIDE.md) - Bescherming van gevoelige bestanden
- [SECURITY.md](SECURITY.md) - Volledige beveiligingsrichtlijnen

## 📝 Bijdragen

Bij het bijdragen aan deze repository:

1. **Check** wat je gaat committen: `git status`
2. **Review** je wijzigingen: `git diff`
3. **Commit** alleen bestanden die veilig gedeeld kunnen worden
4. **Gebruik** duidelijke commit messages in Nederlands of Engels

Voorbeelden van goede commit messages:
```
Add Py4sci chapter 3 summary on numerical integration
Fix typo in Biophysics Assignment 2 parameters
Update Kwantum notes with harmonic oscillator examples
```

## 🎓 Academische Integriteit

Deze repository is bedoeld voor:
- ✅ Persoonlijke studie en notities
- ✅ Samenwerking tussen studenten (waar toegestaan)
- ✅ Referentiemateriaal voor toekomstige studenten

Respecteer:
- ❌ Auteursrechten van docenten en studiemateriaal
- ❌ Academische integriteitsregels van de universiteit
- ❌ Examenreglementen (geen oplossingen van lopende examens)

## 📧 Contact

Voor vragen over deze repository, neem contact op via GitHub issues of discussies.

---

**Universiteit Antwerpen - Bachelor Fysica**
