# NGS Coverage Report Generator

## 📌 Project Overview
This project provides a Python-based tool to generate **Next-Generation Sequencing (NGS) coverage reports** from Sambamba output files.  
It is designed for diagnostic labs performing targeted sequencing panels, such as the **83-gene panel for congenital myopathy and congenital muscular dystrophy**.  

The tool highlights **genes and exons with sub-optimal coverage (<100% at 30×)**, which is critical for ensuring reliable variant detection in clinical genomics.

---

## 🧬 Clinical Context
- The lab uses **Agilent SureSelect capture kits** to enrich DNA for target genes, followed by sequencing on an **Illumina NextSeq**.
- Performance metrics require that **every coding base is covered by ≥30 reads (30× coverage)**.
- Sambamba generates exon-level coverage statistics, but clinical scientists need **gene-level summaries** to assess test quality.
- This script:
  - Aggregates exon coverage into **gene-level metrics**.
  - Flags **genes and exons with insufficient coverage**.
  - Produces a **PDF report** with tables and plots for clinician review.

---

## ⚙️ Features
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

---

## 📂 Input
- **Sambamba coverage output file** (example: `NGS148_34_139558_CB_CMCMD_S33_R1_001.sambamba_output.txt`)
- Key column: `percentage30` → indicates % of bases covered at ≥30× per exon.

---

## 📊 Output
- **Console summary** of gene and exon coverage.
- **PDF report (`NGS_coverage_report.pdf`)** containing:
  1. Gene-level coverage summary
  2. Genes with <100% coverage at 30×
  3. Exon-level failing regions
  4. Coverage plot

---

## 🚀 How to Run

### 1. Install Dependencies
Ensure you have Python 3.8+ installed.  
Install required libraries:

```bash
pip install pandas matplotlib reportlab
