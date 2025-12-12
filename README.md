# NGS Coverage Report Generator

## Project Overview
This project provides a Python-based tool to generate **Next-Generation Sequencing (NGS) coverage report** from Sambamba output files.  
It is designed for diagnostic labs performing targeted sequencing panels,for congenital myopathy and congenital muscular dystrophy**.  

The tool highlights **genes and exons with sub-optimal coverage that have less than 100% coverage at 30x**, which is critical for ensuring reliable variant detection in clinical genomics.


##  How to Run
### 1. Install Dependencies
Ensure you have Python 3.8+ installed.  
tkinter
pandas
pathlib 
matplotlib
re
reportlab
os
platform
subprocess
### Running the code
When you run the code , it will ask to selected the sambamba file, once you choose , you will get the NGS coverage report of the file as a PDF.

---

## Clinical Context
- The lab uses **Agilent SureSelect capture kits** to enrich DNA for target genes, followed by sequencing on an **Illumina NextSeq**.
- Performance metrics require that **every coding base is covered by ≥30 reads (30× coverage)**.
- Sambamba generates exon-level coverage statistics, but clinical scientists need **gene-level summaries** to assess test quality.
- This script:
- Takes the sambamba output and generates a report listing any genes that have less than 100% coverage at 30x
  - Aggregates exon coverage into **gene-level metrics**.
  - Flags **genes and exons with LOW coverage**.
  - Produces a **PDF report** with tables and plots for clinician review.

---

## Features
- **Gene-level coverage summary**:
  - Minimum exon coverage per gene
  - Mean % coverage
  - Median % coverage
  - Total number of exons
  - Number of exons <100% coverage
  - `LOW_COVERAGE` flagging
- **Exon-level report**: Lists specific exons failing coverage thresholds.
- **Visualization**: Bar plot of minimum exon coverage per gene.
- **PDF report generation**: Clinician-friendly output with tables and plots.

- Clinical Relevance of Metrics
Min exon coverage per gene → weakest point in sequencing; variants may be missed here.

Mean/Median % coverage → overall sequencing quality for the gene.

Total exons → context for gene size; larger genes are harder to cover uniformly.

Exons <100% coverage → counts problematic regions.

LOW_COVERAGE flagging → alerts clinicians to genes needing caution.

Exon-level report → pinpoints exact failing exons for re-sequencing or validation.

---

##  Input
- **Sambamba coverage output file** (example: `NGS148_34_139558_CB_CMCMD_S33_R1_001.sambamba_output.txt`)
- Key column: `percentage30` → indicates % of bases covered at ≥30× per exon.

---

##  Output
- **Console summary** of gene and exon coverage.
- **PDF report (`NGS_coverage_report.pdf`)** containing:
  1. Gene-level coverage summary
  2. Genes with <100% coverage at 30×
  3. Exon-level failing regions
  4. Coverage plot




