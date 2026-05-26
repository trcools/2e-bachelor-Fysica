from __future__ import annotations

import argparse
import re
import tempfile
from pathlib import Path

import fitz
from rapidocr_onnxruntime import RapidOCR


def clean_formula_text(raw_text: str) -> str:
    text = raw_text.strip()
    text = text.replace("4元E0", r"\varepsilon_0")
    text = text.replace("4元0", r"\varepsilon_0")
    text = text.replace("4E0", r"\varepsilon_0")
    text = text.replace("4元 E0", r"\varepsilon_0")
    text = text.replace("...:", r"\cdots")
    text = text.replace("…", r"\cdots")
    text = text.replace("dl';", r"d\ell'")
    text = text.replace("dq.", r"dq")
    text = re.sub(r"(?<![A-Za-z])([A-Za-z])(\d+)", r"\1_\2", text)
    text = re.sub(r"\s+", " ", text)
    return text.strip(" :;,")


def fragment_score(text: str) -> float:
    letters = sum(1 for char in text if char.isalpha())
    digits = sum(1 for char in text if char.isdigit())
    operators = sum(1 for char in text if char in "=+-*/^_()[]{}<>\\")
    specials = sum(1 for char in text if ord(char) > 127)
    score = letters * 2 + operators * 3 + len(text) * 0.1 - digits * 0.4 - specials * 2
    if re.fullmatch(r"[0-9]+", text):
        score -= 8
    return score


def ocr_formula_block(
    ocr: RapidOCR,
    page: fitz.Page,
    bbox: tuple[float, float, float, float],
    render_scale: int,
) -> str | None:
    clip = fitz.Rect(*bbox)
    try:
        pixmap = page.get_pixmap(matrix=fitz.Matrix(render_scale, render_scale), clip=clip, alpha=False)
    except Exception:
        return None

    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp_file:
        temp_path = Path(tmp_file.name)
        try:
            try:
                png_bytes = pixmap.tobytes("png")
                tmp_file.write(png_bytes)
            except Exception:
                # Some pixmap formats raise when writing PNG bytes directly.
                # Fall back to saving via pixmap.save() and read the file.
                pixmap.save(str(temp_path))
                tmp_file.flush()
                tmp_file.close()
                png_bytes = temp_path.read_bytes()
        except Exception:
            # If anything goes wrong writing the png, clean up and skip OCR for this block.
            temp_path.unlink(missing_ok=True)
            return None

    try:
        result, _elapsed = ocr(str(temp_path))
    finally:
        temp_path.unlink(missing_ok=True)

    if not result:
        return None

    fragments: list[tuple[str, float]] = []
    for _box, text, confidence in result:
        cleaned = clean_formula_text(text)
        if not cleaned:
            continue
        if re.search(r"[\u4e00-\u9fff]", cleaned):
            continue
        fragments.append((cleaned, float(confidence)))

    if not fragments:
        return None

    fragments.sort(key=lambda item: (fragment_score(item[0]) + item[1], len(item[0])), reverse=True)
    best = fragments[0][0]

    if len(best) < 4 and not any(symbol in best for symbol in ("=", "+", "-", "\\", "^", "_")):
        return None

    if not any(symbol in best for symbol in ("=", "+", "-", "\\", "^")):
        extra = next(
            (
                frag
                for frag, _conf in fragments[1:]
                if any(token in frag.lower() for token in ("varepsilon", "epsilon", "ε"))
            ),
            None,
        )
        if extra and extra not in best:
            best = f"{best} {extra}"

    if len(best) > 120:
        return None

    return best


def ocr_full_page(ocr: RapidOCR, page: fitz.Page, render_scale: int) -> str | None:
    """OCR the entire page and return cleaned text or None."""
    try:
        pixmap = page.get_pixmap(matrix=fitz.Matrix(render_scale, render_scale), alpha=False)
    except Exception:
        return None

    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp_file:
        temp_path = Path(tmp_file.name)
        try:
            try:
                png_bytes = pixmap.tobytes("png")
                tmp_file.write(png_bytes)
            except Exception:
                pixmap.save(str(temp_path))
                tmp_file.flush()
                tmp_file.close()
                png_bytes = temp_path.read_bytes()
        except Exception:
            temp_path.unlink(missing_ok=True)
            return None

    try:
        result, _elapsed = ocr(str(temp_path))
    finally:
        temp_path.unlink(missing_ok=True)

    if not result:
        return None

    # Concatenate fragments into a single cleaned paragraph-like text
    texts = []
    for _box, text, _conf in result:
        cleaned = clean_formula_text(text)
        if cleaned:
            texts.append(cleaned)

    if not texts:
        return None

    joined = " ".join(texts)
    joined = re.sub(r"\s+", " ", joined).strip()
    return joined if joined else None


def write_markdown(
    pdf_path: Path,
    output_dir: Path,
    formula_max_width: int,
    formula_max_height: int,
    render_scale: int,
    assets_suffix: str,
    ocr: RapidOCR,
    no_images: bool = False,
    page_ocr: bool = False,
) -> tuple[int, int]:
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{pdf_path.stem}.md"
    assets_dir = output_dir / f"{pdf_path.stem.replace(' ', '_')}{assets_suffix}"
    if not no_images:
        assets_dir.mkdir(parents=True, exist_ok=True)

    doc = fitz.open(pdf_path)
    chunks = [f"# {pdf_path.stem}", "", f"_Source: {pdf_path.name}_", ""]
    image_count = 0
    formula_count = 0

    for page_number, page in enumerate(doc, start=1):
        blocks = page.get_text("dict")["blocks"]
        ordered_blocks = sorted(blocks, key=lambda block: (round(block["bbox"][1], 1), round(block["bbox"][0], 1)))

        chunks.append(f"## Page {page_number}")
        chunks.append("")

        for block in ordered_blocks:
            if block["type"] == 0:
                for line in block.get("lines", []):
                    line_text = "".join(span["text"] for span in line.get("spans", [])).rstrip()
                    if line_text:
                        chunks.append(line_text)
                chunks.append("")
            elif block["type"] == 1:
                x0, y0, x1, y1 = block["bbox"]
                width = x1 - x0
                height = y1 - y0
                if width <= formula_max_width and height <= formula_max_height:
                    formula_text = ocr_formula_block(ocr, page, block["bbox"], render_scale)
                    if formula_text:
                        formula_count += 1
                        chunks.append(f"$$ {formula_text} $$")
                        chunks.append("")
                        continue

                # If the block isn't recognized as a formula, either save it
                # as an image (normal mode) or skip it (no_images mode).
                formula_text = None
                if width <= formula_max_width and height <= formula_max_height:
                    # small blocks were already attempted above; fallthrough
                    pass

                if not no_images:
                    image_count += 1
                    extension = block.get("ext", "png")
                    image_name = f"page{page_number:02d}_img{image_count:03d}.{extension}"
                    (assets_dir / image_name).write_bytes(block["image"])
                    chunks.append(f"![]({assets_dir.name}/{image_name})")
                    chunks.append("")
                else:
                    # skipping image block in no_images mode
                    continue

        # If page produced no textual content (only header), optionally OCR whole page
        # and append the OCR result when requested.
        page_body_only_header = True
        for i in range(len(chunks)-1, -1, -1):
            if chunks[i].startswith("## Page "):
                # reached the page header; if no other content found, page is empty
                break
            if chunks[i] and not chunks[i].startswith("---"):
                page_body_only_header = False
                break
        if page_body_only_header and page_ocr:
            full_text = ocr_full_page(ocr, page, render_scale)
            if full_text:
                chunks.append(full_text)
                chunks.append("")

        chunks.append("---")
        chunks.append("")

    text = "\n".join(chunks)
    text = re.sub(r"\n{3,}", "\n\n", text).strip() + "\n"
    out_path.write_text(text, encoding="utf-8")
    return formula_count, image_count


def collect_pdfs(input_dir: Path, recursive: bool) -> list[Path]:
    if recursive:
        return sorted(input_dir.rglob("*.pdf"))
    return sorted(input_dir.glob("*.pdf"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert PDFs to markdown with OCR for formula-like image blocks."
    )
    parser.add_argument("--input-dir", help="Folder containing PDF files.")
    parser.add_argument(
        "--output-dir",
        help="Output folder for .md and asset folders. Defaults to input-dir.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Also process PDFs in subfolders.",
    )
    parser.add_argument(
        "--formula-max-width",
        type=int,
        default=240,
        help="Max image width treated as candidate formula block.",
    )
    parser.add_argument(
        "--formula-max-height",
        type=int,
        default=90,
        help="Max image height treated as candidate formula block.",
    )
    parser.add_argument(
        "--render-scale",
        type=int,
        default=4,
        help="Render scale for OCR snapshots of formula blocks.",
    )
    parser.add_argument(
        "--assets-suffix",
        default="_assets",
        help="Suffix used for generated image asset folders.",
    )
    parser.add_argument(
        "--no-images",
        action="store_true",
        help="Do not save non-formula image blocks; keep only text and OCR'd formulas.",
    )
    parser.add_argument(
        "--input-file",
        help="Optional: path to a single PDF file to process instead of a directory.",
    )
    parser.add_argument(
        "--page-ocr",
        action="store_true",
        help="When a page has no extractable text, OCR the whole page and include the text.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    # Determine input(s) and output directory.
    if args.input_file:
        single = Path(args.input_file).expanduser().resolve()
        if not single.exists() or not single.is_file():
            raise SystemExit(f"Input file does not exist: {single}")
        pdf_files = [single]
        output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else single.parent
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        if not args.input_dir:
            raise SystemExit("--input-dir is required when --input-file is not provided")
        input_dir = Path(args.input_dir).expanduser().resolve()
        output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else input_dir
        output_dir.mkdir(parents=True, exist_ok=True)

        if not input_dir.exists() or not input_dir.is_dir():
            raise SystemExit(f"Input directory does not exist: {input_dir}")

        pdf_files = collect_pdfs(input_dir, args.recursive)
        if not pdf_files:
            raise SystemExit(f"No PDF files found in: {input_dir}")

    ocr = RapidOCR()

    total_formulas = 0
    total_images = 0
    for pdf_path in pdf_files:
        formulas, images = write_markdown(
            pdf_path=pdf_path,
            output_dir=output_dir,
            formula_max_width=args.formula_max_width,
            formula_max_height=args.formula_max_height,
            render_scale=args.render_scale,
            assets_suffix=args.assets_suffix,
            ocr=ocr,
            no_images=args.no_images,
            page_ocr=args.page_ocr,
        )
        total_formulas += formulas
        total_images += images
        print(f"{pdf_path.name}: {formulas} formulas, {images} images -> {pdf_path.stem}.md")

    print(
        f"Done. Processed {len(pdf_files)} PDFs, converted {total_formulas} formulas, kept {total_images} images."
    )


if __name__ == "__main__":
    main()