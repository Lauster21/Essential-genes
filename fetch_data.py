"""
fetch_data.py – Download real B. subtilis gene data from SubtiWiki
===================================================================
Run this script once to (re)generate data/bacillus_subtilis_genes.csv
with genuine data straight from SubtiWiki (https://subtiwiki.uni-goettingen.de/).

Data sources
------------
* Gene names, locus tags, CDS positions → SubtiWiki REST API v4
  GET https://subtiwiki.uni-goettingen.de/v4/api/gene/
* Essentiality labels → SubtiWiki category "Essential genes"
  GET https://subtiwiki.uni-goettingen.de/v4/api/category/
* GC content per gene → NCBI E-utilities (genome sequence AL009126.3)
* Expression levels → SubtiWiki expression API
  GET https://subtiwiki.uni-goettingen.de/v4/api/expression/

Usage
-----
    pip install requests biopython pandas
    python fetch_data.py

The script writes data/bacillus_subtilis_genes.csv and prints a summary.
"""

import time
import sys
from pathlib import Path

try:
    import requests
except ImportError:
    sys.exit("Please install the 'requests' package:  pip install requests")

try:
    import pandas as pd
except ImportError:
    sys.exit("Please install the 'pandas' package:  pip install pandas")

try:
    from Bio import Entrez, SeqIO
    from Bio.Seq import Seq
except ImportError:
    sys.exit("Please install the 'biopython' package:  pip install biopython")

# ── Configuration ─────────────────────────────────────────────────────────────

SUBTIWIKI_API = "https://subtiwiki.uni-goettingen.de/v4/api"
ESSENTIAL_CATEGORY = "Essential genes"          # SubtiWiki category label
NCBI_GENOME_ACC = "AL009126.3"                  # B. subtilis 168 complete genome
Entrez.email = "user@example.com"               # ← CHANGE THIS to your own e-mail address
                                                #   NCBI requires a valid contact e-mail.

OUTPUT_CSV = Path(__file__).parent / "data" / "bacillus_subtilis_genes.csv"

# ── Helper functions ──────────────────────────────────────────────────────────

def get_all_pages(url: str) -> list:
    """Follow SubtiWiki pagination and return all result objects."""
    results = []
    while url:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        results.extend(data.get("results", []))
        url = data.get("next")          # None when on the last page
        time.sleep(0.1)                 # be polite to the server
    return results


def fetch_essential_locus_tags() -> set:
    """Return the set of locus tags (BSU numbers) classed as essential in SubtiWiki."""
    categories = get_all_pages(f"{SUBTIWIKI_API}/category/?format=json")
    essential_id = None
    for cat in categories:
        if cat.get("title", "").strip() == ESSENTIAL_CATEGORY:
            essential_id = cat["id"]
            break
    if essential_id is None:
        raise RuntimeError(f"Could not find SubtiWiki category '{ESSENTIAL_CATEGORY}'")

    genes_in_cat = get_all_pages(
        f"{SUBTIWIKI_API}/category/{essential_id}/genes/?format=json"
    )
    return {g["locus_tag"] for g in genes_in_cat}


def fetch_expression_map() -> dict:
    """
    Return a dict mapping locus_tag → mean expression value (log2).
    SubtiWiki integrates expression data from the Stülke-group microarray
    compendium; this endpoint may vary – check the SubtiWiki API docs for
    the current path.
    """
    try:
        expr_data = get_all_pages(f"{SUBTIWIKI_API}/expression/?format=json")
        return {e["locus_tag"]: float(e["mean_log2"]) for e in expr_data
                if "locus_tag" in e and "mean_log2" in e}
    except Exception as exc:
        print(f"  [warn] Could not fetch expression data: {exc}")
        return {}


def fetch_interaction_counts() -> dict:
    """
    Return a dict mapping locus_tag → number of protein–protein interactions
    curated in SubtiWiki.
    """
    try:
        interactions = get_all_pages(f"{SUBTIWIKI_API}/interaction/?format=json")
        counts: dict = {}
        for inter in interactions:
            for partner_key in ("gene_a", "gene_b"):
                lt = inter.get(partner_key, {}).get("locus_tag")
                if lt:
                    counts[lt] = counts.get(lt, 0) + 1
        return counts
    except Exception as exc:
        print(f"  [warn] Could not fetch interaction data: {exc}")
        return {}


def gc_content(sequence: str) -> float:
    """Return the fraction of G+C nucleotides in *sequence*."""
    seq = sequence.upper()
    gc = seq.count("G") + seq.count("C")
    return round(gc / len(seq), 4) if seq else 0.0


def fetch_genome_sequence() -> object:
    """Download the B. subtilis 168 genome from NCBI and return a SeqRecord."""
    print(f"  Downloading genome {NCBI_GENOME_ACC} from NCBI …")
    with Entrez.efetch(
        db="nucleotide", id=NCBI_GENOME_ACC, rettype="gb", retmode="text"
    ) as handle:
        record = SeqIO.read(handle, "genbank")
    return record


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    if Entrez.email == "user@example.com":
        sys.exit(
            "ERROR: Please set Entrez.email to your own e-mail address in fetch_data.py "
            "before running the script. NCBI requires a valid contact e-mail."
        )

    print("Fetching gene list from SubtiWiki …")
    genes = get_all_pages(f"{SUBTIWIKI_API}/gene/?format=json")
    print(f"  {len(genes)} genes retrieved.")

    print("Fetching essential-gene set …")
    essential_lts = fetch_essential_locus_tags()
    print(f"  {len(essential_lts)} essential genes found.")

    print("Fetching expression data …")
    expr_map = fetch_expression_map()

    print("Fetching interaction counts …")
    inter_map = fetch_interaction_counts()

    print("Downloading genome sequence for GC-content calculation …")
    genome = fetch_genome_sequence()
    genome_seq = str(genome.seq)

    print("Building dataset …")
    rows = []
    for g in genes:
        locus_tag = g.get("locus_tag", "")
        gene_name = g.get("gene") or locus_tag
        start = g.get("start")
        end = g.get("end")
        strand = g.get("strand", "+")

        if start is None or end is None:
            continue                    # skip genes with no position info

        length_bp = abs(int(end) - int(start)) + 1
        # Extract CDS sequence (1-based genome coordinates → 0-based slice)
        s, e = sorted([int(start) - 1, int(end)])
        cds_seq = genome_seq[s:e]
        if strand == "-":
            # Reverse complement for genes on the minus strand
            cds_seq = str(Seq(cds_seq).reverse_complement())
        gc = gc_content(cds_seq) if cds_seq else None
        expr = expr_map.get(locus_tag)
        interactions = inter_map.get(locus_tag, 0)
        is_ess = 1 if locus_tag in essential_lts else 0

        rows.append({
            "gene":            gene_name,
            "locus_tag":       locus_tag,
            "gene_length_bp":  length_bp,
            "gc_content":      gc,
            "expression_level": expr,
            "num_interactions": interactions,
            "is_essential":    is_ess,
        })

    df = pd.DataFrame(rows)
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSaved {len(df)} genes to {OUTPUT_CSV}")
    print(df["is_essential"].value_counts().rename({0: "non-essential", 1: "essential"}))


if __name__ == "__main__":
    main()
