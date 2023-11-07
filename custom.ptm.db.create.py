import pandas as pd
from uniprotparser.betaparser import UniprotSequence
import re
df = pd.read_csv(r"C:\Users\Toan Phung\Downloads\supp_RA118.000696_136158_1_supp_105958_p6rjww.txt", sep="\t")

pattern = re.compile("r(\d+)")
new_df = []
for i, r in df.iterrows():
    mods = r["Variable modifications identified by spectrum"].split(",")
    for mod in mods:
        mod = mod.strip()
        if mod.endswith(" Citrullinated (+0.98)"):
            match = pattern.search(mod)
            if match:
                for p in r["Protein groups accession numbers"].split(","):
                    header = UniprotSequence(p, True)
                    if header.accession:
                        new_df.append([str(header), int(match.group(1))+r["Peptide start index"]-1, "R"])
                    else:
                        new_df.append([p, int(match.group(1))+r["Peptide start index"]-1, "R"])
                break
                
new_df = pd.DataFrame(new_df, columns=["Protein", "Position", "Residue"])
new_df.to_csv(r"C:\Users\Toan Phung\Downloads\supp_RA118.000696_136158_1_supp_105958_p6rjww_data.txt", sep="\t", index=False)
