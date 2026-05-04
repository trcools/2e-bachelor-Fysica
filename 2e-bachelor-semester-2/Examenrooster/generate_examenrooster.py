from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import date, timedelta
from pathlib import Path


@dataclass(frozen=True)
class DayPlan:
    blocks: list[str]
    notes: str


START_DATE = date(2026, 5, 4)
END_DATE = date(2026, 6, 22)

OUTPUT_PATH = Path(__file__).with_name("examenrooster.csv")

SPECIAL_DAYS: dict[date, DayPlan] = {
    date(2026, 5, 18): DayPlan(
        blocks=[
            "Sterrenstelsels - presentation run-through",
            "Sterrenstelsels - presentation final checks",
            "Sterren en Planeten - light review",
        ],
        notes="Presentation day; avoid new material and focus on execution.",
    ),
    date(2026, 6, 4): DayPlan(
        blocks=[
            "Sterrenstelsels - quiet final recap",
            "Sterrenstelsels - exam day",
        ],
        notes="Exam Sterrenstelsels; only a light warm-up before the test.",
    ),
    date(2026, 6, 9): DayPlan(
        blocks=[
            "Sterren en Planeten - quiet final recap",
            "Sterren en Planeten - exam day",
        ],
        notes="Exam Sterren en Planeten; short review and then rest.",
    ),
    date(2026, 6, 10): DayPlan(
        blocks=[
            "Thermische Fysica - open-book exercise set",
            "Thermische Fysica - formula sheet and approach",
            "Thermische Fysica - exam day",
        ],
        notes="Open-book exam/practice moment Thermische Fysica; focus on structure.",
    ),
    date(2026, 6, 11): DayPlan(
        blocks=[
            "Thermische Fysica - oral core questions",
            "Thermische Fysica - oral explanation practice",
            "Thermische Fysica - exam day",
        ],
        notes="Oral Thermische Fysica; answer clearly and concisely out loud.",
    ),
    date(2026, 6, 18): DayPlan(
        blocks=[
            "Relativity & Electromagnetism - quiet final recap",
            "Relativity & Electromagnetism - exam day",
        ],
        notes="Exam Relativity & Electromagnetism; only short review.",
    ),
    date(2026, 6, 22): DayPlan(
        blocks=[
            "Groups and Representations - quiet final recap",
            "Groups and Representations - exam day",
        ],
        notes="Exam Groups and Representations; minimal warm-up.",
    ),
}

PHASES = [
    (
        date(2026, 5, 4),
        date(2026, 5, 17),
        {
            "weekday": [
                "Sterrenstelsels - theory and core concepts",
                "Sterrenstelsels - course exercises",
                "Sterren en Planeten - theory",
                "Sterrenstelsels - write summary",
                "Sterrenstelsels - presentation slide draft",
            ],
            "weekend": [
                "Sterrenstelsels - practice questions",
                "Sterren en Planeten - basic questions",
                "Sterrenstelsels - presentation structure and sources",
                "Sterrenstelsels - review block",
            ],
            "notes": "Phase 1: build the foundation for Sterrenstelsels and start the presentation.",
        },
    ),
    (
        date(2026, 5, 19),
        date(2026, 5, 31),
        {
            "weekday": [
                "Sterrenstelsels - deepening",
                "Sterren en Planeten - theory",
                "Thermische Fysica - exercises",
                "Relativity & Electromagnetism - introduction",
                "Sterrenstelsels - active recall",
            ],
            "weekend": [
                "Sterrenstelsels - practice questions",
                "Sterren en Planeten - exercises",
                "Groups and Representations - basics",
                "Thermische Fysica - summary",
            ],
            "notes": "Phase 2: after the presentation, broaden into the June exams.",
        },
    ),
    (
        date(2026, 6, 1),
        date(2026, 6, 3),
        {
            "weekday": [
                "Sterrenstelsels - timed set",
                "Sterrenstelsels - formula sheet",
                "Sterren en Planeten - maintenance",
                "Sterrenstelsels - weak spots",
                "Sterrenstelsels - short recall",
            ],
            "weekend": [
                "Sterrenstelsels - timed set",
                "Sterrenstelsels - error analysis",
                "Sterren en Planeten - light review",
                "Sterrenstelsels - short recall",
            ],
            "notes": "Phase 3A: final build-up to the Sterrenstelsels exam.",
        },
    ),
    (
        date(2026, 6, 5),
        date(2026, 6, 8),
        {
            "weekday": [
                "Sterren en Planeten - theory",
                "Sterren en Planeten - exercises",
                "Thermische Fysica - open-book set",
                "Sterren en Planeten - mock exam",
                "Sterren en Planeten - summary",
            ],
            "weekend": [
                "Sterren en Planeten - practice questions",
                "Thermische Fysica - formulas and notes",
                "Sterren en Planeten - recall",
                "Sterren en Planeten - error analysis",
            ],
            "notes": "Phase 3B: peak effort on Sterren en Planeten and parallel Thermische Fysica.",
        },
    ),
    (
        date(2026, 6, 12),
        date(2026, 6, 17),
        {
            "weekday": [
                "Relativity & Electromagnetism - theory",
                "Relativity & Electromagnetism - exercises",
                "Relativity & Electromagnetism - derivations",
                "Groups and Representations - introduction",
                "Relativity & Electromagnetism - summary",
            ],
            "weekend": [
                "Relativity & Electromagnetism - practice questions",
                "Groups and Representations - basics",
                "Relativity & Electromagnetism - error analysis",
                "Relativity & Electromagnetism - recall",
            ],
            "notes": "Phase 4: shift focus toward Relativity & Electromagnetism.",
        },
    ),
    (
        date(2026, 6, 19),
        date(2026, 6, 21),
        {
            "weekday": [
                "Groups and Representations - theory",
                "Groups and Representations - exercises",
                "Groups and Representations - proof practice",
                "Groups and Representations - mock exam",
                "Groups and Representations - recall",
            ],
            "weekend": [
                "Groups and Representations - practice questions",
                "Groups and Representations - summary",
                "Groups and Representations - error analysis",
                "Groups and Representations - final recall",
            ],
            "notes": "Phase 5: final stretch toward Groups and Representations.",
        },
    ),
]


def phase_for(current_date: date) -> dict[str, list[str] | str]:
    for phase_start, phase_end, plan in PHASES:
        if phase_start <= current_date <= phase_end:
            return plan
    raise ValueError(f"No phase configured for {current_date}")


def build_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    current = START_DATE
    while current <= END_DATE:
        special = SPECIAL_DAYS.get(current)
        if special is not None:
            blocks = special.blocks
            notes = special.notes
        else:
            plan = phase_for(current)
            template_key = "weekend" if current.weekday() >= 5 else "weekday"
            blocks = list(plan[template_key])
            notes = str(plan["notes"])

        padded_blocks = blocks + [""] * (5 - len(blocks))
        rows.append(
            {
                "date": current.isoformat(),
                "weekday": current.strftime("%A"),
                "blocks": str(len(blocks)),
                "total_hours": f"{len(blocks) * 1.5:.1f}",
                "main_focus": blocks[0].split(" - ")[0],
                "block_1": padded_blocks[0],
                "block_2": padded_blocks[1],
                "block_3": padded_blocks[2],
                "block_4": padded_blocks[3],
                "block_5": padded_blocks[4],
                "notes": notes,
            }
        )
        current += timedelta(days=1)
    return rows


def main() -> None:
    rows = build_rows()
    with OUTPUT_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "date",
                "weekday",
                "blocks",
                "total_hours",
                "main_focus",
                "block_1",
                "block_2",
                "block_3",
                "block_4",
                "block_5",
                "notes",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()