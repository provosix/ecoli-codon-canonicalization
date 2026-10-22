#!/usr/bin/env python3
"""Analyze synonymous-codon canonicalization in E. coli K-12 MG1655.

The script reads a GenBank genome, filters unsuitable CDS features, counts sense
codons, builds three canonical-codon models, validates protein preservation, and
exports reproducible tables, FASTA files, a JSON summary, and a comparison plot.
"""

from __future__ import annotations

import argparse
import json
import random
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter


# Keep the analysis settings in one place so the workflow is easy to reproduce.
GENETIC_CODE_ID = 11  # NCBI bacterial/plastid code used for E. coli.
RANDOM_SEED = 42  # Fix the random model so repeated runs produce the same mapping.
STOP_CODONS = {"TAA", "TAG", "TGA"}  # Standard stop codons in code table 11.
DNA_BASES = set("ACGT")  # Only unambiguous DNA bases are accepted in this analysis.


# Return the first annotation value, or an empty string if it is missing.
def qualifier(feature: Any, key: str) -> str:
    """Read one GenBank qualifier safely."""
    values = feature.qualifiers.get(key, [""])
    return str(values[0]) if values else ""


# Identify Sec annotations because TGA can encode selenocysteine in context.
def contains_selenocysteine_annotation(feature: Any) -> bool:
    """Return True when a CDS has an annotated Sec translation exception."""
    return any("sec" in value.lower() for value in feature.qualifiers.get("transl_except", []))


# Filter annotations that do not fit a simple, standard codon-by-codon model.
def extract_valid_cds(record: Any, table_id: int) -> tuple[list[dict[str, Any]], Counter, int]:
    """Extract valid CDS records and report why other CDS features were excluded."""
    all_cds = [feature for feature in record.features if feature.type == "CDS"]
    exclusions: Counter = Counter()
    valid_records: list[dict[str, Any]] = []

    for feature in all_cds:
        # Pseudogenes do not represent ordinary intact protein-coding sequences.
        if "pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers:
            exclusions["pseudogene"] += 1
            continue

        # Programmed frameshifts need a specialized model and are excluded here.
        if "ribosomal_slippage" in feature.qualifiers:
            exclusions["ribosomal_slippage"] += 1
            continue

        # Sec-containing CDSs use a context-dependent exception to the codon table.
        if contains_selenocysteine_annotation(feature):
            exclusions["selenocysteine annotation"] += 1
            continue

        # Extract the feature in its annotated strand orientation.
        dna = str(feature.extract(record.seq)).upper()

        # A complete CDS must consist of whole codons and unambiguous DNA bases.
        if not dna or len(dna) % 3 != 0:
            exclusions["invalid CDS length"] += 1
            continue
        if not set(dna).issubset(DNA_BASES):
            exclusions["ambiguous or non-DNA symbols"] += 1
            continue

        # cds=True checks the start, reading frame, and terminal stop codon.
        try:
            protein = str(Seq(dna).translate(table=CodonTable.unambiguous_dna_by_id[table_id], cds=True))
        except Exception:
            exclusions["translation validation failed"] += 1
            continue

        # Store sense codons only; the terminal stop is preserved separately.
        codons = [dna[index:index + 3] for index in range(0, len(dna), 3)]
        if codons[-1] not in STOP_CODONS:
            exclusions["missing terminal stop"] += 1
            continue
        sense_codons = codons[:-1]

        # Keep the original annotated protein for a later validation check.
        annotated_protein = qualifier(feature, "translation")
        if annotated_protein and annotated_protein != protein:
            exclusions["translation annotation mismatch"] += 1
            continue

        valid_records.append(
            {
                "feature": feature,
                "gene": qualifier(feature, "gene"),
                "locus_tag": qualifier(feature, "locus_tag"),
                "product": qualifier(feature, "product"),
                "dna": dna,
                "protein": protein,
                "sense_codons": sense_codons,
                "terminal_stop": codons[-1],
            }
        )

    return valid_records, exclusions, len(all_cds)


# Build amino-acid-to-codon groups using the selected NCBI genetic code.
def build_codon_groups(table: Any) -> dict[str, list[str]]:
    """Group sense codons by the amino acid they encode."""
    groups: dict[str, list[str]] = defaultdict(list)
    for codon, amino_acid in table.forward_table.items():
        groups[amino_acid].append(codon)
    return dict(groups)


# Select one representative codon per amino acid under each model.
def build_mappings(
    codon_counts: Counter, codon_groups: dict[str, list[str]], seed: int
) -> dict[str, dict[str, str]]:
    """Create most-frequent, seeded-random, and least-frequent codon mappings."""
    most_mapping = {
        amino_acid: max(codons, key=lambda codon: codon_counts[codon])
        for amino_acid, codons in sorted(codon_groups.items())
    }
    least_mapping = {
        amino_acid: min(codons, key=lambda codon: codon_counts[codon])
        for amino_acid, codons in sorted(codon_groups.items())
    }

    # Preserve codon-table order and sort amino acids for deterministic random picks.
    rng = random.Random(seed)
    random_mapping = {
        amino_acid: rng.choice(codons)
        for amino_acid, codons in sorted(codon_groups.items())
    }

    return {
        "Most frequent": most_mapping,
        "Random": random_mapping,
        "Least frequent": least_mapping,
    }


# Count nucleotide composition without relying on precomputed workbook values.
def gc_percent(sequence: str) -> float:
    """Calculate GC percentage for a DNA sequence."""
    if not sequence:
        return 0.0
    return 100.0 * (sequence.count("G") + sequence.count("C")) / len(sequence)


# Calculate codon, nucleotide, GC, and position-level changes for one mapping.
def evaluate_mapping(
    records: list[dict[str, Any]], mapping: dict[str, str], table: Any
) -> tuple[dict[str, Any], list[SeqRecord]]:
    """Apply a codon mapping and summarize sequence changes across all valid CDSs."""
    original_sense = "".join(codon for item in records for codon in item["sense_codons"])
    canonical_sense_parts: list[str] = []
    canonical_fasta_records: list[SeqRecord] = []
    changed_codons = 0
    changed_bases = 0
    position_changes = [0, 0, 0]
    protein_checks_passed = 0
    # Process each CDS separately so the original terminal stop can be retained.
    for item in records:
        original_codons = item["sense_codons"]
        canonical_codons: list[str] = []

        for codon_index, codon in enumerate(original_codons):
            # Preserve the start-methionine rule by using ATG at every CDS start.
            if codon_index == 0:
                replacement = "ATG"
            else:
                amino_acid = table.forward_table[codon]
                replacement = mapping[amino_acid]
            canonical_codons.append(replacement)
            if codon != replacement:
                changed_codons += 1
            for position, (old_base, new_base) in enumerate(zip(codon, replacement)):
                if old_base != new_base:
                    changed_bases += 1
                    position_changes[position] += 1

        canonical_dna = "".join(canonical_codons) + item["terminal_stop"]
        canonical_protein = str(Seq(canonical_dna).translate(table=CodonTable.unambiguous_dna_by_id[GENETIC_CODE_ID], cds=True))
        if canonical_protein != item["protein"]:
            raise RuntimeError(
                f"Protein-preservation check failed for {item['locus_tag']} ({item['gene']})."
            )
        protein_checks_passed += 1
        canonical_sense_parts.append("".join(canonical_codons))

        # Use stable IDs so each generated FASTA record can be traced to its CDS.
        record_id = item["locus_tag"] or item["gene"] or f"CDS_{len(canonical_fasta_records) + 1}"
        description = f"{item['gene']} | {item['product']} | canonicalized model"
        canonical_fasta_records.append(SeqRecord(Seq(canonical_dna), id=record_id, description=description))

    canonical_sense = "".join(canonical_sense_parts)
    total_codons = len(original_sense) // 3
    total_bases = len(original_sense)
    total_position_changes = sum(position_changes)
    original_gc = gc_percent(original_sense)
    canonical_gc = gc_percent(canonical_sense)

    # Keep both per-position rates and each position's share of all substitutions.
    position_metrics = {}
    for index, changed in enumerate(position_changes, start=1):
        position_metrics[f"position_{index}_changed_count"] = changed
        position_metrics[f"position_{index}_changed_percent"] = 100.0 * changed / total_codons if total_codons else 0.0
        position_metrics[f"position_{index}_share_of_substitutions_percent"] = (
            100.0 * changed / total_position_changes if total_position_changes else 0.0
        )

    metrics: dict[str, Any] = {
        "changed_codons": changed_codons,
        "total_sense_codons": total_codons,
        "changed_codons_percent": 100.0 * changed_codons / total_codons if total_codons else 0.0,
        "changed_nucleotides": changed_bases,
        "total_sense_nucleotides": total_bases,
        "changed_nucleotides_percent": 100.0 * changed_bases / total_bases if total_bases else 0.0,
        "original_gc_percent": original_gc,
        "canonical_gc_percent": canonical_gc,
        "gc_change_percentage_points": canonical_gc - original_gc,
        "protein_preservation_checks_passed": protein_checks_passed,
        **position_metrics,
    }
    return metrics, canonical_fasta_records


# Compute codon counts from the filtered dataset and prepare all report tables.
def run_analysis(input_path: Path, output_dir: Path, seed: int) -> None:
    """Run the complete analysis and export all results to the output directory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    record = SeqIO.read(str(input_path), "genbank")
    table = CodonTable.unambiguous_dna_by_id[GENETIC_CODE_ID]
    records, exclusions, raw_cds_count = extract_valid_cds(record, GENETIC_CODE_ID)
    if not records:
        raise RuntimeError("No valid CDS features remained after filtering.")

    # Count codons only in the validated sense regions, excluding terminal stops.
    codon_counts = Counter(codon for item in records for codon in item["sense_codons"])
    codon_groups = build_codon_groups(table)
    mappings = build_mappings(codon_counts, codon_groups, seed)

    # Confirm that the dataset covers all 61 sense codons in genetic code table 11.
    missing_codons = sorted(set(table.forward_table) - set(codon_counts))
    if missing_codons:
        print("Warning: no observations for these sense codons:", ", ".join(missing_codons))

    # Analyze each model and write a FASTA file containing canonicalized CDS records.
    model_results: dict[str, dict[str, Any]] = {}
    for model_name, mapping in mappings.items():
        metrics, fasta_records = evaluate_mapping(records, mapping, table)
        model_results[model_name] = metrics
        safe_name = model_name.lower().replace(" ", "_")
        SeqIO.write(fasta_records, str(output_dir / f"canonicalized_cds_{safe_name}.fasta"), "fasta")

    # Prepare the baseline data for the report and JSON summary.
    original_sense = "".join(codon for item in records for codon in item["sense_codons"])
    baseline = {
        "genome_accession": record.id,
        "genome_length_bp": len(record.seq),
        "raw_cds_features": raw_cds_count,
        "valid_cds_features": len(records),
        "excluded_cds_features": raw_cds_count - len(records),
        "sense_codon_count": len(original_sense) // 3,
        "sense_nucleotide_count": len(original_sense),
        "original_gc_percent": gc_percent(original_sense),
        "genetic_code_id": GENETIC_CODE_ID,
        "random_seed": seed,
        "excluded_by_reason": dict(exclusions),
        "protein_preservation_checks_passed_per_model": len(records),
        "sense_codon_space": len(table.forward_table),
        "amino_acid_count": len(codon_groups),
        "sense_codon_space_reduction_percent": 100.0 * (len(table.forward_table) - len(codon_groups)) / len(table.forward_table),
    }

    # Export all numeric results and model mappings in two readable Excel files.
    write_workbooks(output_dir, baseline, model_results, codon_counts, codon_groups, mappings, seed)

    # Export machine-readable results so later scripts can reuse the analysis output.
    summary = {"baseline": baseline, "models": model_results, "mappings": mappings}
    (output_dir / "results_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    # Generate a compact comparison figure from the freshly computed metrics.
    make_comparison_plot(output_dir, model_results)

    # Print a short run summary so users can check the main outputs quickly.
    print(f"Genome: {record.id} ({len(record.seq):,} bp)")
    print(f"CDS features: {raw_cds_count:,} raw; {len(records):,} retained; {raw_cds_count - len(records):,} excluded")
    print(f"Sense codons: {baseline['sense_codon_count']:,}; baseline GC: {baseline['original_gc_percent']:.2f}%")
    for model_name, metrics in model_results.items():
        print(
            f"{model_name}: {metrics['changed_codons_percent']:.2f}% codons changed; "
            f"{metrics['changed_nucleotides_percent']:.2f}% nucleotides changed; "
            f"GC {metrics['canonical_gc_percent']:.2f}%"
        )
    print(f"Outputs saved to: {output_dir.resolve()}")


# Create formatted Excel workbooks with results, QC information, and codon tables.
def write_workbooks(
    output_dir: Path,
    baseline: dict[str, Any],
    model_results: dict[str, dict[str, Any]],
    codon_counts: Counter,
    codon_groups: dict[str, list[str]],
    mappings: dict[str, dict[str, str]],
    seed: int,
) -> None:
    """Write the main-results and supplementary-table Excel workbooks."""
    main_path = output_dir / "main_results.xlsx"
    wb = Workbook()
    ws = wb.active
    if ws is None:
        raise RuntimeError("Could not create the main-results worksheet.")
    ws.title = "Main Results"
    ws.append(["Metric", "Most frequent", "Random", "Least frequent"])

    # Add the same clearly defined metrics for every model.
    result_rows = [
        ("Changed sense codons (%)", "changed_codons_percent"),
        ("Changed sense nucleotides (%)", "changed_nucleotides_percent"),
        ("Original GC (%)", "original_gc_percent"),
        ("Canonical GC (%)", "canonical_gc_percent"),
        ("GC change (percentage points)", "gc_change_percentage_points"),
        ("Position 1 changes (% of codons)", "position_1_changed_percent"),
        ("Position 2 changes (% of codons)", "position_2_changed_percent"),
        ("Position 3 changes (% of codons)", "position_3_changed_percent"),
        ("Position 1 share of substitutions (%)", "position_1_share_of_substitutions_percent"),
        ("Position 2 share of substitutions (%)", "position_2_share_of_substitutions_percent"),
        ("Position 3 share of substitutions (%)", "position_3_share_of_substitutions_percent"),
        ("Protein-preservation checks passed", "protein_preservation_checks_passed"),
    ]
    for label, key in result_rows:
        ws.append([label] + [model_results[name][key] for name in mappings])

    # Add dataset and filtering details on a separate sheet.
    design = wb.create_sheet("Baseline & Design")
    design.append(["Item", "Value"])
    for key, value in baseline.items():
        if key == "excluded_by_reason":
            for reason, count in value.items():
                design.append([f"Excluded: {reason}", count])
        else:
            design.append([key.replace("_", " "), value])

    # List all observed codon counts and the representative codon for each model.
    codon_ws = wb.create_sheet("Codon Usage")
    codon_ws.append(["Codon", "Amino acid", "Count", "Frequency among sense codons (%)", "Most frequent", "Random", "Least frequent"])
    total_codons = sum(codon_counts.values())
    table = CodonTable.unambiguous_dna_by_id[GENETIC_CODE_ID]
    for codon in sorted(table.forward_table):
        amino_acid = table.forward_table[codon]
        selected = [mappings[name][amino_acid] for name in mappings]
        codon_ws.append([
            codon,
            amino_acid,
            codon_counts[codon],
            100.0 * codon_counts[codon] / total_codons if total_codons else 0.0,
            *selected,
        ])

    # List one chosen codon for every amino acid under each model.
    map_ws = wb.create_sheet("Codon Mappings")
    map_ws.append(["Amino acid", "Most frequent", "Random", "Least frequent"])
    for amino_acid in sorted(codon_groups):
        map_ws.append([amino_acid] + [mappings[name][amino_acid] for name in mappings])
    map_ws.append([])
    map_ws.append(["Random seed", seed])

    style_workbook(wb)
    wb.save(main_path)

    # Create a supplementary workbook with codon-level and amino-acid-level tables.
    supp_path = output_dir / "supplementary_tables.xlsx"
    supp = Workbook()
    table1 = supp.active
    if table1 is None:
        raise RuntimeError("Could not create the supplementary worksheet.")
    table1.title = "Supplementary Table 1"
    table1.append(["Codon", "Amino acid", "Count", "Frequency (%)", "Most frequent", "Random", "Least frequent"])
    for codon in sorted(table.forward_table):
        amino_acid = table.forward_table[codon]
        table1.append([
            codon,
            amino_acid,
            codon_counts[codon],
            100.0 * codon_counts[codon] / total_codons if total_codons else 0.0,
            mappings["Most frequent"][amino_acid],
            mappings["Random"][amino_acid],
            mappings["Least frequent"][amino_acid],
        ])

    table2 = supp.create_sheet("Supplementary Table 2")
    table2.append(["Amino acid", "Synonymous codons", "Most frequent", "Random", "Least frequent"])
    for amino_acid in sorted(codon_groups):
        table2.append([
            amino_acid,
            ", ".join(codon_groups[amino_acid]),
            mappings["Most frequent"][amino_acid],
            mappings["Random"][amino_acid],
            mappings["Least frequent"][amino_acid],
        ])
    style_workbook(supp)
    supp.save(supp_path)


# Apply consistent formatting so the spreadsheets are readable when opened.
def style_workbook(workbook: Workbook) -> None:
    """Format worksheet headers, widths, filters, and frozen panes."""
    header_fill = PatternFill(fill_type="solid", fgColor="1F4E78")
    header_font = Font(color="FFFFFF", bold=True)

    for worksheet in workbook.worksheets:
        worksheet.freeze_panes = "A2"
        worksheet.auto_filter.ref = worksheet.dimensions
        for cell in worksheet[1]:
            cell.fill = header_fill
            cell.font = header_font
            cell.alignment = Alignment(vertical="center", wrap_text=True)
        worksheet.row_dimensions[1].height = 30

        # Size each column to fit its contents, with a reasonable maximum width.
        for column_cells in worksheet.columns:
            max_length = max((len(str(cell.value or "")) for cell in column_cells), default=0)
            column_index = column_cells[0].column
            if column_index is None:
                continue
            column_letter = get_column_letter(column_index)
            worksheet.column_dimensions[column_letter].width = min(max(max_length + 2, 12), 48)

        # Keep long text readable without making every row excessively tall.
        for row in worksheet.iter_rows(min_row=2):
            for cell in row:
                cell.alignment = Alignment(vertical="top", wrap_text=True)


# Plot model comparisons using only values calculated during this run.
def make_comparison_plot(output_dir: Path, model_results: dict[str, dict[str, Any]]) -> None:
    """Save a bar chart comparing codon changes, nucleotide changes, and GC content."""
    model_names = list(model_results)
    changed_codons = [model_results[name]["changed_codons_percent"] for name in model_names]
    changed_nt = [model_results[name]["changed_nucleotides_percent"] for name in model_names]
    gc_values = [model_results[name]["canonical_gc_percent"] for name in model_names]

    figure, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    panels = [
        (axes[0], changed_codons, "Sense codons changed (%)"),
        (axes[1], changed_nt, "Sense nucleotides changed (%)"),
        (axes[2], gc_values, "GC content after canonicalization (%)"),
    ]
    for axis, values, title in panels:
        axis.bar(model_names, values)
        axis.set_title(title)
        axis.tick_params(axis="x", rotation=25)
        axis.grid(axis="y", alpha=0.25)
        axis.set_axisbelow(True)
    figure.tight_layout()
    figure.savefig(output_dir / "model_comparison.png", dpi=200, bbox_inches="tight")
    plt.close(figure)


# Parse command-line arguments so the input and output paths are configurable.
def parse_args() -> argparse.Namespace:
    """Read input, output, and random-seed options from the command line."""
    parser = argparse.ArgumentParser(
        description="Analyze synonymous-codon canonicalization in E. coli K-12 MG1655."
    )
    parser.add_argument("--input", type=Path, default=Path("ecoli.gb"), help="Input GenBank file (default: ecoli.gb).")
    parser.add_argument("--output-dir", type=Path, default=Path("ecoli_codon_results"), help="Directory for generated results.")
    parser.add_argument("--seed", type=int, default=RANDOM_SEED, help="Seed for the random codon model.")
    return parser.parse_args()


# Run the workflow only when this file is executed directly.
if __name__ == "__main__":
    arguments = parse_args()
    run_analysis(arguments.input, arguments.output_dir, arguments.seed)
