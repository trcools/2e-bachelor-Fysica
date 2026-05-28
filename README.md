# Bachelor Physics - Study Materials

This repo contains study material, summaries, code, and assignments for the bachelor of Physics. It is educational material organized per semester and course.

## Structure

### 2nd Bachelor - Semester 1
- **`2e-bachelor-semester-1/Biophysics/`** - Biophysics
- **`2e-bachelor-semester-1/Kwantum/`** - Quantum Mechanics
- **`2e-bachelor-semester-1/Py4sci/`** - Python for Scientists
- **`2e-bachelor-semester-1/VFR/`** - Vector and Function Spaces
- **`2e-bachelor-semester-1/Statistiek/`** - Statistics (currently empty)

### 2nd Bachelor - Semester 2
- **`Groepen-en-Representaties/`** - Groups and representations
- **`Materiaal-Fysica/`** - Materials Physics
- **`REM/`** - Relativity and Electromagnetism
- **`Sterrenstelsels/`** - Galaxies
- **`Thermische-Fysica/`** - Thermal Physics
- **`exp2/`** - Experimental Physics 2

### Other Folders
- **`1e-bachelor-semester-1/`** - Courses from 1st bachelor

## Security & Privacy

**Read [SECURITY.md](SECURITY.md) for detailed security guidelines.**

- **Private files** (slides, syllabi, textbooks): keep them in the `_lokaal/` folder — Git ignores it automatically
- **Do not commit copyrighted material** (professor slides, textbooks, course books)
- **Do not commit personal information** or sensitive data, including your own

### Installing pre-commit hooks

```bash
pip install pre-commit
pre-commit install
```

This automatically checks for:
- Large files (>5MB)
- Private keys or credentials
- Merge conflicts
- YAML/JSON syntax

## Setup

### Python Dependencies

```bash
pip install -r requirements.txt
```

## Getting Started

```bash
# 1. (Optional) create and activate a virtual environment
python3 -m venv .venv
source .venv/bin/activate

# 2. Install dependencies
pip install -r requirements.txt
```

To work with the notebooks, open them in VS Code or Jupyter (for example, `jupyter lab`).

## Gitignore note

Private files belong in `_lokaal/`. The `.gitignore` also blocks `*.pdf`, `*.docx`, and `*.pptx` by default.

## Tests

```bash
pytest
```

## Contributing

Contributions to this repository are always welcome. **Commit** only files that are safe to share and **Use** clear commit messages.

This repository is intended for:
- Personal study and notes
- Collaboration between students
- Reference material for future students
