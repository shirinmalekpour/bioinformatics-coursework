"""Download the coding sequence (CDS) of a RefSeq transcript from NCBI and save it as FASTA.

Usage:
    python fetch_cds.py NM_000053.4

Requires Biopython (pip install biopython).
NCBI asks for a contact email with every Entrez request, so replace the placeholder below.
"""
import sys

from Bio import Entrez, SeqIO

Entrez.email = "shirin.malelpou@gmail.com"  # replace with your own address


def fetch_cds_sequence(refseq_id):
    """Fetch a GenBank record, cut out its first CDS feature and write it to <refseq_id>.fasta."""
    handle = Entrez.efetch(db="nucleotide", id=refseq_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    # Find the CDS features of the transcript
    cds_features = [feature for feature in record.features if feature.type == "CDS"]
    if not cds_features:
        print(f"No CDS found for {refseq_id}")
        return None

    # Take the first CDS if there is more than one
    cds = cds_features[0]
    cds_record = record[cds.location.start:cds.location.end]
    cds_record.id = f"{refseq_id}_CDS"
    cds_record.description = f"CDS for {record.description}"

    output_file = f"{refseq_id}.fasta"
    SeqIO.write(cds_record, output_file, "fasta")
    print(f"Saved CDS for {refseq_id} to {output_file}")
    return cds_record


if __name__ == "__main__":
    accession = sys.argv[1] if len(sys.argv) > 1 else "NM_000053.4"
    fetch_cds_sequence(accession)
