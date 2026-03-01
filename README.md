# Bachelor Fysica - Studiemateriaal

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
- **`Materiaal-Fysica/`** - Materiaal Fysica
- **`REM/`** - Relativiteit en Elektromagnetisme
- **`Sterrenstelsels/`** - Astrofysica: sterrenstelsels en kosmologie
- **`Thermische-Fysica/`** - Statistische mechanica en thermodynamica
- **`project-exp-2/`** - Experimenteel project semester 2

### Andere Folders
- **`1e-bachelor-semester-1/`** - Vakken uit 1e bachelor (o.a. lineaire algebra)
- **`Cursussen:boeken/`** - Cursussen en (studie)boeken
- **`GIT-commands/`** - Git referentiemateriaal (samenvatting hoe Git te gebruiken)


## ⚠️ Beveiliging & Privacy

**Lees de [SECURITY.md](SECURITY.md) voor uitgebreide beveiligingsrichtlijnen.**

### Quick Reference:
- ✅ **Commit**: samenvattingen, code, notebooks, opdrachten
- ✅ **Gebruik pre-commit hooks** voor automatische controles
- 🔒 **Privé bestanden** (slides, syllabi, leerboeken): bewaar ze in de `_lokaal/` map — Git negeert die automatisch
- ❌ **Geen auteursrechtelijk materiaal** committen (professorslides, leerboeken, cursusboeken)
- ❌ **Geen persoonlijke informatie** of gevoelige data committen (mag op eigen risico)

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


## 📝 Bijdragen

Bij het bijdragen aan deze repository:

1. **Check** wat je gaat committen: `git status`
2. **Review** je wijzigingen: `git diff`
3. **Commit** alleen bestanden die veilig gedeeld kunnen worden
4. **Gebruik** duidelijke commit messages in Nederlands of Engels


## 🎓 Academische Integriteit

Deze repository is bedoeld voor:
- Persoonlijke studie en notities
- Samenwerking tussen studenten (waar toegestaan)
- Referentiemateriaal voor toekomstige studenten

## 📧 Contact

Voor vragen over deze repository, neem contact op via GitHub issues of discussies.

---

**Bachelor Fysica**
