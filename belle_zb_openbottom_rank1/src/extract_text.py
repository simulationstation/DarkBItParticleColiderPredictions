#!/usr/bin/env python3
"""
Extract text from Belle open-bottom paper, focusing on:
1. Table I (fit results)
2. Supplementary material (binned data)
"""

import fitz
import os
import re

PDF_PATH = "../data/papers/belle_zb_openbottom_1512.07419.pdf"

def extract_all_text():
    """Extract text from all pages."""
    doc = fitz.open(PDF_PATH)

    for page_num in range(len(doc)):
        page = doc[page_num]
        text = page.get_text()

        print(f"\n{'='*60}")
        print(f"PAGE {page_num + 1}")
        print('='*60)
        print(text)

    doc.close()

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    extract_all_text()
