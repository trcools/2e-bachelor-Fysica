# Bachelor Fysica - Studiemateriaal

Deze repository bevat studiemateriaal, samenvattingen, code en opdrachten voor de bachelor Fysica.

- **Markdown samenvattingen** 
- **Jupyter notebooks**
- **Python code** 
- **Opdrachten en examens**

**Dit is GEEN software project** — het is educatief materiaal georganiseerd per semester en vak.

## 📂 Structuur

### 2e Bachelor - Semester 1
- **`Biophysics/`** - Biofysica 
- **`Kwantum/`** - Kwantummechanica 
- **`Py4sci/`** - Python for Scientists
- **`VFR/`** - Vector en Functieruimten

### 2e Bachelor - Semester 2
- **`Groepen-en-Representaties/`** - Groepen en representaties
- **`Materiaal-Fysica/`** - Materiaal Fysica
- **`REM/`** - Relativiteit en Elektromagnetisme
- **`Sterrenstelsels/`** - Sterrenstelsels
- **`Thermische-Fysica/`** - Thermische Fysica
- **`exp2/`** - Experimenteren in de fysica 2

### Andere Folders
- **`1e-bachelor-semester-1/`** - Vakken uit 1e bachelor
- **`Cursussen:boeken/`** - Cursussen en (studie)boeken


## ⚠️ Beveiliging & Privacy

**Lees de [SECURITY.md](SECURITY.md) voor uitgebreide beveiligingsrichtlijnen.**

- **Commit**: samenvattingen, code, notebooks, opdrachten
- **Gebruik pre-commit hooks** voor automatische controles
- **Geen persoonlijke informatie** of gevoelige data committen (eigen verantwoordelijkheid)

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

Voor numerieke notebooks:

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

Bij het bijdragen aan deze repository is altijd welkom, maar **Check** wat je gaat committen: `git status`, **Review** je wijzigingen: `git diff`, **Commit** alleen bestanden die veilig gedeeld kunnen worden en **Gebruik** duidelijke commit messages.

Deze repository is bedoeld voor:
- Persoonlijke studie en notities
- Samenwerking tussen studenten
- Referentiemateriaal voor toekomstige studenten

