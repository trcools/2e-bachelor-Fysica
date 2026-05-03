#!/usr/bin/env python3
"""
Extract comprehensive content from summary PDFs and expand the exam sheet.
This script systematically extracts formulas, proofs, and key definitions
from each chapter PDF and populates the exam sheet LaTeX file.
"""

import re
from pypdf import PdfReader
from pathlib import Path

# Base directory
BASE_DIR = Path(__file__).parent.parent / "Samenvatting-REM"

# Mapping of chapters to their content templates
chapter_mapping = {
    "H2": {
        "file": "H2 Electrostatics.pdf",
        "latex_section": "H2 Electrostatica",
        "keywords": ["Coulomb", "Gauss", "potential", "divergence", "curl"],
    },
    "H3": {
        "file": "H3 Potentials.pdf",
        "latex_section": "H3 Potentialen",
        "keywords": ["Poisson", "Laplace", "boundary", "image", "multipl"],
    },
    "H4": {
        "file": "H4 Electric fields in matter.pdf",
        "latex_section": "H4 E-velden in materie",
        "keywords": ["polarization", "bound charge", "displacement"],
    },
    "H5": {
        "file": "H5 Magnetostatics.pdf",
        "latex_section": "H5 Magnetostatica",
        "keywords": ["Lorentz", "Biot", "Ampere", "vector potential"],
    },
    "H6": {
        "file": "H6 Magnetic fields in matter.pdf",
        "latex_section": "H6 B-velden in materie",
        "keywords": ["magnetization", "bound current", "H field"],
    },
    "H7": {
        "file": "H7 Electrodynamics.pdf",
        "latex_section": "H7 Elektrodynamica",
        "keywords": ["Faraday", "induction", "Maxwell", "current"],
    },
    "H8": {
        "file": "H8 Conservation laws.pdf",
        "latex_section": "H8 Conservatiewetten",
        "keywords": ["continuity", "Poynting", "energy", "momentum", "stress"],
    },
    "H9": {
        "file": "H9 Electromagnetic waves.pdf",
        "latex_section": "H9 EM-golven",
        "keywords": ["wave", "plane wave", "intensity", "propagation"],
    },
    "H10": {
        "file": "H10 Potentials and fields.pdf",
        "latex_section": "H10 Potentialen en velden",
        "keywords": ["gauge", "Lorenz", "retarded", "wave equation"],
    },
    "H11": {
        "file": "H11 Radiation.pdf",
        "latex_section": "H11 Straling",
        "keywords": ["dipole", "radiation", "power", "far-field"],
    },
    "H12": {
        "file": "H12 Electrodynamics and relativity.pdf",
        "latex_section": "H12 EM en relativiteit",
        "keywords": ["Lorentz", "relativity", "transformation", "covariant"],
    },
}


def extract_text_from_pdf(filepath):
    """Extract text from PDF file"""
    try:
        with open(filepath, 'rb') as f:
            reader = PdfReader(f)
            text = ""
            for page in reader.pages:
                text += page.extract_text()
            return text
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return ""


def extract_formulas_and_key_lines(text, keywords, max_lines=50):
    """Extract lines containing formulas or keywords from text"""
    lines = text.split('\n')
    extracted = []
    
    for i, line in enumerate(lines):
        # Look for lines with math symbols or keywords
        if any(kw.lower() in line.lower() for kw in keywords):
            # Include context lines
            start = max(0, i - 1)
            end = min(len(lines), i + 3)
            for j in range(start, end):
                if lines[j].strip() and lines[j] not in extracted:
                    extracted.append(lines[j].strip())
                    if len(extracted) >= max_lines:
                        break
        
        # Also look for lines with mathematical symbols
        if any(sym in line for sym in ['$', '=', '∇', '×', '·', '∂', '∫']):
            if line.strip() and line not in extracted:
                extracted.append(line.strip())
                if len(extracted) >= max_lines:
                    break
    
    return extracted[:max_lines]


def generate_latex_content_for_chapter(chapter_key, text):
    """Generate LaTeX content for a chapter from extracted text"""
    info = chapter_mapping[chapter_key]
    extracted = extract_formulas_and_key_lines(text, info["keywords"], max_lines=40)
    
    # Format as LaTeX bullet points
    latex_content = []
    for line in extracted:
        # Clean up line
        cleaned = line.strip()
        if cleaned and len(cleaned) > 3:
            # Escape LaTeX special characters if needed
            latex_content.append(f"• {cleaned[:100]}")  # Limit line length for space
    
    return latex_content[:15]  # Max 15 bullet points per chapter


def main():
    print("=" * 70)
    print("REM Exam Sheet - Comprehensive Content Extractor")
    print("=" * 70)
    
    all_chapter_content = {}
    
    for chapter_key in sorted(chapter_mapping.keys()):
        info = chapter_mapping[chapter_key]
        pdf_path = BASE_DIR / info["file"]
        
        if pdf_path.exists():
            print(f"\n✓ Processing {chapter_key}...")
            text = extract_text_from_pdf(str(pdf_path))
            content_lines = generate_latex_content_for_chapter(chapter_key, text)
            all_chapter_content[chapter_key] = content_lines
            
            # Show sample
            print(f"  Found {len(content_lines)} key lines for {chapter_key}")
            if content_lines:
                print(f"  Sample: {content_lines[0][:70]}...")
        else:
            print(f"✗ File not found: {pdf_path}")
    
    # Now we have all content, let's generate the expanded LaTeX
    print("\n" + "=" * 70)
    print("Generating expanded LaTeX document...")
    print("=" * 70)
    
    # Store content for use in next step
    with open('chapter_content.txt', 'w') as f:
        import json
        json.dump(all_chapter_content, f, indent=2)
    
    print("\n✓ Content extraction complete. Stored in chapter_content.txt")
    print("Ready to generate expanded exam sheet...")


if __name__ == "__main__":
    main()
