#!/usr/bin/env python3
"""
Fetch strain and genome accession for a list of protein IDs from NCBI.

Usage:
    python fetch_protein_metadata.py protein_ids.txt output.csv

Where:
    - protein_ids.txt = one protein accession per line (e.g. WP_306294993.1)
    - output.csv      = output table with metadata
"""

import sys
import csv
import time
import re
from Bio import Entrez, SeqIO
from Bio.SeqFeature import SeqFeature

# >>>>> IMPORTANT: set your email here (NCBI requirement) <<<<<
Entrez.email = "lachlan.black.2022@uni.strath.ac.uk"
Entrez.tool = "fetch_protein_metadata_for_tree"

# polite rate limiting (NCBI suggests <= 3 requests/sec)
DELAY_BETWEEN_REQUESTS = 0.4


def parse_genome_accession_from_coded_by(coded_by_value: str) -> str | None:
    """
    Extract a genomic accession (e.g. NC_000913.3, CP000000.1, NZ_CP000000.1)
    from a 'coded_by' qualifier string.

    Examples of coded_by:
        "NC_000913.3:190..255"
        "complement(NC_000913.3:190..255)"
        "join(NZ_CP009713.1:123..456,NZ_CP009713.1:789..999)"
    """
    if not coded_by_value:
        return None

    # Look for a typical genomic accession pattern: AAA_########.#
    match = re.search(r"([A-Z]{1,3}_[0-9]+\.[0-9]+)", coded_by_value)
    if match:
        return match.group(1)

    return None


def fetch_protein_metadata(protein_id: str) -> dict:
    result = {
        "protein_id": protein_id,
        "organism": "",
        "strain": "",
        "isolate": "",
        "serovar": "",
        "biosample": "",
        "bioproject": "",
        "genome_accession": "",
        "coded_by_raw": "",
        "error": ""
    }

    try:
        handle = Entrez.efetch(db="protein", id=protein_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
    except Exception as e:
        result["error"] = f"EFETCH/parse failed: {e}"
        return result

    result["organism"] = record.annotations.get("organism", "")

    # SOURCE FEATURE (strain, isolate, BioSample)
    source_feat = None
    for feat in record.features:
        if feat.type == "source":
            source_feat = feat
            break

    if source_feat:
        quals = source_feat.qualifiers

        def q(name):
            return quals.get(name, [""])[0]

        result["strain"] = q("strain")
        result["isolate"] = q("isolate")
        result["serovar"] = q("serovar")
        result["biosample"] = q("db_xref").replace("BioSample:", "") if "BioSample:" in "".join(quals.get("db_xref", [""])) else ""
        result["bioproject"] = q("db_xref").replace("BioProject:", "") if "BioProject:" in "".join(quals.get("db_xref", [""])) else ""

    # CDS FEATURE: coded_by
    genome_acc = None
    for feat in record.features:
        if feat.type == "CDS":
            coded_by = feat.qualifiers.get("coded_by", [""])[0]
            result["coded_by_raw"] = coded_by
            genome_acc = parse_genome_accession_from_coded_by(coded_by)
            break

    # If coded_by failed — try Assembly search using BioSample
    if not genome_acc and result["biosample"]:
        try:
            asm_search = Entrez.esearch(db="assembly", term=f"BioSample:{result['biosample']}")
            asm_ids = Entrez.read(asm_search)["IdList"]
            if asm_ids:
                asm_summary = Entrez.esummary(db="assembly", id=asm_ids[0])
                doc = Entrez.read(asm_summary)
                ftp = doc['DocumentSummarySet']['DocumentSummary'][0].get("FtpPath_RefSeq", "")
                # Extract accession
                if ftp:
                    genome_acc = ftp.split("/")[-1]  # last item is the accession prefix
        except:
            pass

    if genome_acc:
        result["genome_accession"] = genome_acc
    else:
        result["error"] = result["error"] + "; No genome accession found"

    return result

def main():
    if len(sys.argv) != 3:
        print("Usage: python fetch_protein_metadata.py protein_ids.txt output.csv")
        sys.exit(1)

    ids_file = sys.argv[1]
    out_file = sys.argv[2]

    # Read protein IDs
    with open(ids_file) as f:
        protein_ids = [line.strip() for line in f if line.strip()]

    print(f"Loaded {len(protein_ids)} protein IDs from {ids_file}")

    rows = []
    for i, pid in enumerate(protein_ids, start=1):
        print(f"[{i}/{len(protein_ids)}] Fetching {pid} ...", end="", flush=True)
        meta = fetch_protein_metadata(pid)
        rows.append(meta)
        print(" done" if not meta["error"] else f" ERROR ({meta['error']})")
        time.sleep(DELAY_BETWEEN_REQUESTS)

    # Write CSV
    fieldnames = [
        "protein_id",
        "organism",
        "strain",
        "isolate",
        "serovar",
        "genome_accession",
        "coded_by_raw",
        "error",
    ]

    with open(out_file, "w", newline="") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"\nDone. Wrote metadata for {len(rows)} proteins to {out_file}")


if __name__ == "__main__":
    main()
