# coding exercise
# Libraries
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import re

from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.lib.units import inch

import os
import platform
import subprocess
from tkinter import Tk, filedialog

# function to open the report pdf in any os the user has 
def open_file_cross_platform(filepath):
    system = platform.system()

    try:
        if system == "Windows":
            os.startfile(filepath)
        elif system == "Darwin":  
            subprocess.run(["open", filepath])
        else:  
            subprocess.run(["xdg-open", filepath])
    except Exception as e:
        print(f"Could not auto-open file: {e}")

# function which gets the file from the user 
def get_sambamba_file():
    Tk().withdraw()  # hide root window
    file_path = filedialog.askopenfilename(
        title="Select Sambamba Coverage File",
        filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")]
    )
    if not file_path:
        raise SystemExit("No file selected. Exiting program.")
    return file_path

file_path = get_sambamba_file()
print(f"Selected file: {file_path}")

# file_path = r"C:\Users\adith\Downloads\NGS148_34_139558_CB_CMCMD_S33_R1_001.sambamba_output.txt"

def loading_sambamba_file(file_path):
    """
    this function loads any Sambamba-style file regardless of spacing issues, mixed delimiters,
    duplicated column names, or non-standard header formatting.
    """
    with open(file_path, "r") as f:
        raw_header = f.readline().strip()    # Reading the raw header

    fixed_header = re.sub(r" +", "\t", raw_header)    # - Replacing spaces with a single tab limited
    
    columns = fixed_header.split("\t")              

    columns = [standardize_name(c) for c in columns]

    columns = make_unique(columns)   #  Make the duplicate columns present into unique columns 

    df = pd.read_csv(
        file_path,
        sep=r"\s+",           # spliting on any space
        names=columns,        
        header=None,          
        skiprows=1
    )

    
# Spliting gene+accession into two seperate columns
    df = splitting_genesymbol_accession(df)

    return df


def splitting_genesymbol_accession(df):
    target_cols = [c for c in df.columns if "genesymbol" in c and "accession" in c]

    if not target_cols:
        return df

    col = target_cols[0]

    split = df[col].str.split(";", n=1, expand=True)

    df["genesymbol"] = split[0]
    df["accession"] = split[1]

    df = df.drop(columns=[col])

    return df


# this function Standardizes the column names

def standardize_name(name):
    name = name.strip().lower()
    name = re.sub(r"[#]", "", name)
    name = re.sub(r"[^a-z0-9]+", "_", name)    
    name = re.sub(r"_+", "_", name)            
    name = name.strip("_")                     
    return name

# this function makes duplicate column names unique

def make_unique(columns):
    seen = {}
    unique = []
    for col in columns:
        if col not in seen:
            seen[col] = 1
            unique.append(col)
        else:
            seen[col] += 1
            unique.append(f"{col}_{seen[col]}")
    return unique

df = loading_sambamba_file(file_path)
print(df.columns)
print(df.head())

# takes the sambamba file above output and generates a report listing any genes that have less than 100% coverage at 30x
def gene_coverage(df):
    """
    Creates a gene coverage report with 
      - min exon coverage per gene
      - mean % coverage
      - median % coverage
      - total number of exons
      - number of exons <100%
      - LOW_COVERAGE flagging 
      and it lists any genes that have less than 100% coverage at 30x
    """

    df["percentage30"] = pd.to_numeric(df["percentage30"], errors="coerce")
    grouped = df.groupby("genesymbol")

    gene_cov = pd.DataFrame({
    "min_percentage30": grouped["percentage30"].min(),
    "mean_percentage30": grouped["percentage30"].mean(),
    "median_percentage30": grouped["percentage30"].median(),
    "total_exons": grouped.size(),
    "exons_below_100": grouped["percentage30"].agg(lambda x: (x < 100).sum())
}).reset_index()

    gene_cov["status"] = gene_cov["min_percentage30"].apply(
        lambda x: "OK" if x >= 100 else "LOW_COVERAGE"
    )

    # genes that have less than 100% coverage at 30x
    low_cov = gene_cov[gene_cov["min_percentage30"] < 100]

    return gene_cov, low_cov

print("\nGENE COVERAGE REPORT")
full_report, low_genes = gene_coverage(df)
print(full_report)

print("\nGENES WITH LOW COVERAGE (<100%)")
print(low_genes)

# this function produces an exon-level low coverage idea .
def exon_coverage(df):

    required = ["chromosome", "startposition", "endposition", "fullposition",
                "genesymbol", "meancoverage", "percentage30"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df["meancoverage"] = pd.to_numeric(df["meancoverage"], errors="coerce")
    df["percentage30"] = pd.to_numeric(df["percentage30"], errors="coerce")

    df["exon_size"] = df["endposition"] - df["startposition"] + 1

    df["exon_status"] = df["percentage30"].apply(
        lambda x: "OK" if x >= 100 else "LOW_COVERAGE"
    )

    failing_exons = df[df["percentage30"] < 100].copy()

    minimal_cols = [
        "chromosome",
        "startposition",
        "endposition",
        "fullposition",
        "genesymbol",
        "exon_size",
        "meancoverage",
        "percentage30",
        "exon_status"
    ]

    failing_exons = failing_exons[minimal_cols].reset_index(drop=True)

    return failing_exons



print("\nEXON-LEVEL LOW COVERAGE REPORT (<100%)")
exon_low = exon_coverage(df)
print(exon_low)

# this function Displays a bar plot summarizing minimum exon %30x per gene.
def plot_gene_coverage(gene_cov):


    gene_cov = gene_cov.sort_values("min_percentage30", ascending=True)

    colors = gene_cov["min_percentage30"].apply(
        lambda x: "red" if x < 100 else "green"
    )

    plt.figure(figsize=(14, 6))
    plt.bar(gene_cov["genesymbol"], gene_cov["min_percentage30"], color=colors)
    plt.axhline(100, color="black", linestyle="--", linewidth=1)

    plt.ylabel("Minimum Percentage ≥30×")
    plt.xlabel("Gene")
    plt.title("Genes that have less than 100% coverage at 30x(red)")

    plt.xticks(rotation=90)
    plt.tight_layout()

    # plt.show()
    plt.savefig("gene_coverage_plot.png", dpi=300)
    plt.close()

plot_gene_coverage(full_report)

# this function makes the pdf with the results i am getting from the file
def create_pdf_report(
    gene_cov, low_genes, exon_low, plot_path="gene_coverage_plot.png"
):

    pdf_path = "NGS_coverage_report.pdf"

    doc = SimpleDocTemplate(pdf_path, pagesize=A4)
    styles = getSampleStyleSheet()
    story = []

    # Title
    title = Paragraph(
        "<b>NGS Coverage Report</b>",
        styles["Title"],
    )
    story.append(title)
    story.append(Spacer(1, 20))

    # Gene Coverage Summary
    story.append(Paragraph("<b>1. Gene-Level Coverage Summary</b>", styles["Heading2"]))
    story.append(Spacer(1, 10))

    gene_table_data = [
        ["Gene", "Min %30x", "Mean %30x", "Median %30x", "Total Exons", "Exons <100%", "Status"]
    ] + gene_cov.values.tolist()

    gene_table = Table(gene_table_data)
    gene_table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
                ("GRID", (0, 0), (-1, -1), 0.5, colors.black),
                ("ALIGN", (1, 1), (-1, -1), "CENTER"),
            ]
        )
    )

    story.append(gene_table)
    story.append(Spacer(1, 20))

    # Low Coverage Genes
    story.append(Paragraph("<b>2. Genes that have less than 100% coverage at 30x</b>", styles["Heading2"]))
    story.append(Spacer(1, 10))

    if not low_genes.empty:
        low_table_data = [
            ["Gene", "Min %30x", "Mean %30x", "Median %30x", "Total Exons", "Exons <100%", "Status"]
        ] + low_genes.values.tolist()

        low_table = Table(low_table_data)
        low_table.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.lightpink),
                    ("GRID", (0, 0), (-1, -1), 0.5, colors.black),
                    ("ALIGN", (1, 1), (-1, -1), "CENTER"),
                ]
            )
        )
        story.append(low_table)
    else:
        story.append(Paragraph("No low-coverage genes detected.", styles["Normal"]))

    story.append(Spacer(1, 20))

    # Exon-Level Low Coverage
    story.append(Paragraph("<b>3. Exon-Level Low Coverage - Exons that have less than 100% coverage at 30x</b>", styles["Heading2"]))
    story.append(Spacer(1, 10))

    if not exon_low.empty:
        exon_table_data = [
            ["Chr", "Start", "End", "Position", "Gene", "Size", "MeanCov", "%30x", "Status"]
        ] + exon_low.values.tolist()

        exon_table = Table(exon_table_data, repeatRows=1)
        exon_table.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.lightgoldenrodyellow),
                    ("GRID", (0, 0), (-1, -1), 0.5, colors.black),
                    ("ALIGN", (5, 1), (-1, -1), "CENTER"),
                ]
            )
        )
        story.append(exon_table)
    else:
        story.append(Paragraph("No failing exons detected.", styles["Normal"]))

    story.append(Spacer(1, 20))

    # Plot
    story.append(Paragraph("<b>4. Gene Coverage Summary Plot</b>", styles["Heading2"]))
    story.append(Spacer(1, 10))

    story.append(Image(plot_path, width=6 * inch, height=3 * inch))

    doc.build(story)
    print(f"PDF generated: {pdf_path}")

    open_file_cross_platform(pdf_path)


plot_gene_coverage(full_report)

create_pdf_report(
    full_report,
    low_genes,
    exon_low,
    plot_path="gene_coverage_plot.png"
)




