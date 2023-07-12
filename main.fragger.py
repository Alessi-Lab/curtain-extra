import pandas as pd
from uniprotparser.betaparser import UniprotSequence
from sequal.sequence import Sequence
import re

reg_positon_residue = re.compile("_(\w)(\d+)")

# a function to read fasta file into a dictionary where the key is the uniprot accession id extracted from the fasta
def read_fasta(fasta_file: str) -> dict[str, dict[str, str]]:
    fasta_dict = {}
    with open(fasta_file, 'r') as f:
        current_acc = ""
        for line in f:
            if line.startswith('>'):
                acc = UniprotSequence(line.strip(), True)

                if acc.accession:
                    fasta_dict[str(acc)] = ""
                    current_acc = str(acc)
                else:
                    fasta_dict[line.strip().replace(">", "")] = ""
                    current_acc = line.strip().replace(">", "")

            else:
                fasta_dict[current_acc] += line.strip()
    return fasta_dict


if __name__ == "__main__":
    fasta_file = r"C:\Users\Toan Phung\Downloads\2023-07-05-decoys-reviewed-contam-UP000000589.fas"
    fasta_dict = read_fasta(fasta_file)
    df = pd.read_csv(r"C:\Users\Toan Phung\Downloads\abundance_single-site_MD (1).tsv", sep='\t')
    for i, row in df.iterrows():
        match = reg_positon_residue.search(row["Index"])
        if match:
            position = int(match.group(2))
            residue = match.group(1)
            position_in_peptide = position
            if row["ProteinID"] in fasta_dict:
                print(row["ProteinID"])
                print(row["Peptide"])
                peptide_seq = row["Peptide"].split(";")[0].upper()
                print(fasta_dict[row["ProteinID"]])
                try:
                    peptide_position = fasta_dict[row["ProteinID"]].index(peptide_seq)
                except ValueError:
                    peptide_position = fasta_dict[row["ProteinID"]].replace("I", "L").index(peptide_seq.replace("I", "L"))
                    df.at[i, "Comment"] = "I replaced by L"
                if peptide_position >= -1:
                    position_in_peptide = position - peptide_position
            df.at[i, "Position"] = position
            df.at[i, "Residue"] = residue
            df.at[i, "Position.in.peptide"] = position_in_peptide

    df.to_csv(r"C:\Users\Toan Phung\Downloads\abundance_single-site_MD (1)_with_position.txt", sep='\t', index=False)
