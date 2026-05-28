# Security Guidelines

This repository contains study material for the 2nd-year bachelor in Physics. To keep things safe, respectful, and academically clean, please follow these guidelines. Thanks for helping keep the repo healthy for everyone.

## What NOT to commit

1. **Personal or sensitive information**
   - Student IDs
   - Private notes with sensitive info
   - API keys or passwords
2. **Copyrighted material**
   - Professor lecture slides
   - Instructor syllabi or course books
   - Commercially published textbooks
   - Current-year exam questions
3. **Large binary files**
   - Lecture videos
   - Large datasets (use external hosting)

## What IS OK to commit

1. **Your own work**
   - Self-made summaries
   - Your own code and scripts
   - Your own notes and worked solutions
2. **Open-source material**
   - Jupyter notebooks with your analyses
   - Python scripts for calculations
   - Documentation and README files

---

## Quick commit checklist

```bash
# 1. See what changed
git status

# 2. Review the diff
git diff

# 3. Add only safe files
git add <file>

# 4. Commit with a clear message
git commit -m "Short description of the change"
```

## Common commands

**Remove a file from tracking (keep it locally):**

```bash
git rm --cached file.pdf
echo "file.pdf" >> .gitignore
git commit -m "Stop tracking file.pdf"
```

**Run pre-commit on everything:**

```bash
pre-commit run --all-files
```

---

## Keep PDFs and documents private

Want to store PDFs (slides, syllabi, textbooks) without making them public? You have three options:

### Option 1 — `_lokaal/` folder (recommended for local use)

There is a special `_lokaal/` folder in this repository. Anything placed there is **never** tracked by Git.

```bash
# Copy a PDF into the local private folder:
cp ~/Downloads/cursusthermische.pdf _lokaal/Thermische-Fysica/

# Git will ignore this file:
git status  # _lokaal/*.pdf will NOT appear in the output
```

See [`_lokaal/README.md`](_lokaal/README.md) for a recommended folder structure.

### Option 2 — Make the repository private

Go to **GitHub → Settings → Danger Zone → Change visibility → Make private**.  
Then you can store everything in the repository and only you (and invited collaborators) can see it.

### Option 3 — Separate private repository

Create a second, private GitHub repository (for example `2e-bachelor-Fysica-lokaal`) and keep private files there. The public repository stays clean for your own work.

---

## .gitignore configuration

The `.gitignore` automatically blocks:

- `*.pdf`, `*.docx`, `*.pptx` — document files
- `_lokaal/**` — everything in the private folder

To commit one specific PDF anyway (for example, your own summary):

```bash
git add -f mijn_samenvatting.pdf
```

## Pre-commit hooks

This repository uses pre-commit hooks that automatically check for:

- Large files (warnings for >500KB)
- Private keys
- Merge conflicts
- Python syntax errors

Install them with:

```bash
pip install pre-commit
pre-commit install
```

## If you accidentally committed sensitive material

1. **Remove the file from tracking (keep it locally):**
   ```bash
   git rm --cached bestand.pdf
   ```
2. **Add it to .gitignore:**
   ```bash
   echo "bestand.pdf" >> .gitignore
   ```
3. **Commit the change:**
   ```bash
   git add .gitignore
   git commit -m "Remove sensitive file from tracking"
   ```
4. **If it was already pushed:** Contact the repository owner or use `git filter-branch` or BFG Repo-Cleaner.
