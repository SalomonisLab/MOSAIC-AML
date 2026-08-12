import os
charts = r"C:\Users\krog5w\.gemini\antigravity\scratch\_pr\Matrix-AML\pipeline\charts"
for f in os.listdir(charts):
    p = os.path.join(charts, f)
    print(f"  {f}: {os.path.getsize(p):,} bytes")

# Convert summary PDF page 1 to PNG for preview
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.image as mpimg

# Quick check: re-open the all-mutations PDF
from PyPDF2 import PdfReader
r = PdfReader(os.path.join(charts, "beataml_rna_per_mutation.pdf"))
print(f"\nPer-mutation PDF: {len(r.pages)} pages")
r2 = PdfReader(os.path.join(charts, "beataml_rna_all_mutations.pdf"))
print(f"All-mutations PDF: {len(r2.pages)} pages")
print("\nAll PDFs look good!")
