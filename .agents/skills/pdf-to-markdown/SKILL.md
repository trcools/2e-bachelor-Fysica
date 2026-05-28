---
name: pdf-to-markdown
description: >
  Convert a PDF file to a Markdown (.md) file, preserving document structure,
  inline math ($...$), display math ($$...$$), code blocks, and figure placeholders.
  Use this skill whenever the user wants to convert, transcribe, extract, or reformat
  a PDF into Markdown — especially for equation-heavy documents like academic papers,
  physics or math lecture notes, or textbooks. Trigger on phrases like "convert this PDF
  to markdown", "extract this paper to md", "transcribe this PDF", "reformat this PDF",
  or any time a .pdf file is mentioned alongside markdown, equations, or text extraction.
  Even if the user just says "turn this PDF into something I can edit", use this skill.
---
 
# PDF → Markdown Conversion Skill
 
## Overview
 
This skill converts PDFs to Markdown using **vision-based page transcription**.
Text extraction is deliberately avoided for math-heavy documents because PDF glyph
encoding mangles equations. Instead, each page is rasterized and read visually.
 
---
 
## Step 0 — Understand the document
 
Run a quick diagnostic before anything else:
 
```bash
pdfinfo <input.pdf>
pdffonts <input.pdf>
```
 
Note:
- **Page count** — determines chunking strategy (see Step 2)
- **Font presence** — if `pdffonts` returns nothing, the PDF is scanned/raster;
  vision is the only option regardless
- **Is it a multi-column layout?** (papers usually are) — note this; it affects
  how you read the rasterized pages
---
 
## Step 1 — Rasterize pages
 
Use `pdftoppm` to convert pages to images at 200 DPI (higher than usual to keep
equations sharp and readable):
 
```bash
mkdir -p /tmp/pdf_pages
pdftoppm -jpeg -r 200 -f 1 -l LAST_PAGE INPUT_PDF /tmp/pdf_pages/page
ls /tmp/pdf_pages/
```
 
`pdftoppm` zero-pads filenames based on total page count. Always `ls` to confirm
the exact filenames before reading them.
 
---
 
## Step 2 — Transcription strategy by document size
 
| Pages | Strategy |
|-------|----------|
| 1–10  | Transcribe all pages in one pass, sequentially |
| 11–30 | Transcribe in batches of 5–8 pages; assemble at the end |
| 31+   | Warn the user about token cost; offer to do a page range, or proceed with explicit confirmation |
 
For large documents, surface the warning clearly:
> "This PDF has N pages. Full transcription will use roughly N × 1,600 tokens
> for images. Shall I proceed with all pages, or a specific range (e.g. pages 1–20)?"
 
---
 
## Step 3 — Transcribe each page
 
For each rasterized page image, view it and produce clean Markdown. Apply the
following rules consistently:
 
### Math
- **Inline math** (expressions within a sentence): `$...$`
- **Display math** (standalone equations, derivations, numbered equations): `$$...$$`
- Preserve equation numbers when visible, e.g.:
  ```
  $$
  E = mc^2 \tag{1}
  $$
  ```
- Use standard LaTeX notation: `\frac{}{}`, `\int`, `\sum`, `\partial`,
  `\nabla`, `\mathbf{}`, `\hat{}`, `\vec{}`, Greek letters (`\alpha`, `\beta`, …), etc.
- For aligned multi-line derivations, use the `align` environment:
  ```
  $$
  \begin{align}
    \nabla \cdot \mathbf{E} &= \frac{\rho}{\varepsilon_0} \\
    \nabla \times \mathbf{B} &= \mu_0 \mathbf{J}
  \end{align}
  $$
  ```
 
### Document structure
- Map section headings to `#`, `##`, `###` etc., following the document hierarchy
- Preserve **bold** and *italic* emphasis where clearly present
- Preserve bullet lists and numbered lists
- Preserve tables using GFM table syntax
### Figures and images
- When a figure (plot, diagram, photograph, schematic) is encountered, insert:
  ```
  ![Figure N: <short description>](figure_N.png)
  ```
  where N increments globally across the whole document.
- Do not attempt to reproduce figure content in text unless it is a simple
  diagram you can describe meaningfully.
### Code
- Wrap code snippets in fenced code blocks with a language tag where identifiable:
  ````
  ```python
  # code here
  ```
  ````
 
### What to skip
- Running headers/footers (page numbers, journal names, author names repeated
  at top/bottom)
- Watermarks
- Decorative horizontal rules unrelated to section structure
---
 
## Step 4 — Assemble and save
 
Concatenate all page transcriptions in order. Separate pages with a blank line
(do **not** insert `---` page breaks unless the document itself has clear chapter
breaks).
 
Save the result to the same directory as the input PDF, with the same base name:
```
<original_name>.md
```
 
Or if the user specifies a different output path, use that.
 
---
 
## Step 5 — Quality checkpoint
 
After assembling, scan the output for common transcription errors:
 
- [ ] Inline math accidentally written as display math or vice versa
- [ ] Broken LaTeX (unmatched braces, missing `\end{...}`)
- [ ] Figure numbering is sequential and consistent
- [ ] Section heading hierarchy makes sense (no `###` appearing before any `##`)
- [ ] No garbled characters from PDF encoding leakage
If any issues are found, fix them before presenting the output to the user.
 
---
 
## Step 6 — Present the output
 
Copy the final `.md` file to `/mnt/user-data/outputs/` and present it with
`present_files`. Give the user a brief summary:
- Number of pages converted
- Number of figures placeholdered
- Any caveats (e.g. pages that were ambiguous, multi-column layouts that may
  need manual review, handwritten annotations that were skipped)
