# Examenrooster

Scripts en data gerelateerd aan examenroosters.

- Doel: code en CSV's voor het genereren en analyseren van examenroosters.
- Voorbeelden: `examenrooster.csv`, `generate_examenrooster.py`.
# Examenrooster

This folder contains a small generator and the first CSV roster for the exam period.

## Files

- `generate_examenrooster.py` - creates the CSV file from fixed phases and special exam days.
- `examenrooster.csv` - the generated roster.

## Regenerate

```bash
python Examenrooster/generate_examenrooster.py
```

## Scope

The roster runs from 2026-05-04 through 2026-06-22 and includes:

- a daily study focus
- 90-minute study blocks
- special days for the presentation and exams