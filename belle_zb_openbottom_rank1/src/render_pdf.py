#!/usr/bin/env python3
"""
Render Belle open-bottom Zb paper PDF pages to high-resolution PNGs.
Focus on pages containing the open-bottom (B B*, B* B*) spectra.
"""

import fitz  # PyMuPDF
import os

PDF_PATH = "../data/papers/belle_zb_openbottom_1512.07419.pdf"
OUTPUT_DIR = "../data/figures"
DPI = 600

def render_all_pages():
    """Render all pages to identify figure locations."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    doc = fitz.open(PDF_PATH)
    print(f"PDF has {len(doc)} pages")

    for page_num in range(len(doc)):
        page = doc[page_num]

        # Get page text to identify content
        text = page.get_text()

        # Check for open-bottom related content
        keywords = ['BB*', 'B*B*', 'B̄*', 'missing mass', 'Mmiss', 'open-bottom',
                   'Fig.', 'Figure', 'spectrum', 'Zb(10610)', 'Zb(10650)']

        has_keyword = any(kw.lower() in text.lower() for kw in keywords)

        # Render at high DPI
        zoom = DPI / 72
        matrix = fitz.Matrix(zoom, zoom)
        pix = page.get_pixmap(matrix=matrix)

        out_path = f"{OUTPUT_DIR}/page{page_num+1:02d}_{DPI}dpi.png"
        pix.save(out_path)

        print(f"Page {page_num+1}: saved to {out_path}")
        if has_keyword:
            print(f"  -> Contains relevant keywords")

        # Print first 500 chars for context
        preview = text[:500].replace('\n', ' ')
        print(f"  Text preview: {preview[:200]}...")

    doc.close()
    print(f"\nAll pages rendered to {OUTPUT_DIR}")

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    render_all_pages()
