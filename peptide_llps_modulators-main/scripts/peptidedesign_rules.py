import random
import os

AMINO_ACIDS = list("AGILPVFWYDERHKSTCMNQ")

REPEATS = False  # AA repeats off
SOLUBLE_PEPTIDES = True  # A rule of thumb in designing soluble peptides is to ensure that at least 1/5 AA residues is charged

# Peptides with over 75% of these AAs can form intermolecular hydrogen bonds (crosslinks), leading to gel formation in aqueous solutions.
LIMITED_AA = set("DEHKNQRSTY")

# Aspartic acid can undergo hydrolysis and cause peptide cleavage under acidic conditions when paired with glycine, proline or serine.
FORBIDDEN_D_NEIGHBORS = {"G", "P", "S"}

# N-terminal Q is unstable due to cyclization under acidic conditions, while N-terminal N should be avoided as its protecting group is hard to remove during cleavage.
FORBIDDEN_NTERM = {"Q", "N"}

# A series of these AAs can form β-sheets, leading to incomplete solvation during synthesis and causing deletions.
FORBIDDEN_NEIGHBORS = set("VIYFLQT")

# Multiple of these AA residues in a sequence can cause significant deletions during synthesis
MAX_ONCE_AA = {"P", "S", "C", "M"}

# Path to txt file
desktop_path = os.path.join(os.path.expanduser("~"), "Documents", "random_peptides.txt")

def generate_peptide(length, existing_peptides):
    while True:
        peptide = []

        while len(peptide) < length:
            next_aa = random.choice(AMINO_ACIDS, 6) #in een keer maken 

            if len(peptide) == 0 and next_aa in FORBIDDEN_NTERM:
                continue

            if not REPEATS and next_aa in peptide:
                continue

            if len(peptide) >= 1 and next_aa == "D" and peptide[-1] in FORBIDDEN_D_NEIGHBORS:
                continue
            if len(peptide) >= 1 and peptide[-1] == "D" and next_aa in FORBIDDEN_D_NEIGHBORS:
                continue

            if len(peptide) >= 1 and peptide[-1] in FORBIDDEN_NEIGHBORS and next_aa in FORBIDDEN_NEIGHBORS:
                continue  

            if next_aa in MAX_ONCE_AA and peptide.count(next_aa) >= 1:
                continue  

            peptide.append(next_aa)

        
        if peptide[-1] == "W": # From our experience, peptides ending with W are often very low in (or have no) purity
            continue  

        if SOLUBLE_PEPTIDES and sum(1 for aa in peptide if aa in "DERHK") < length // 5:
            continue  

        limited_count = sum(1 for aa in peptide if aa in LIMITED_AA)
        if (limited_count / length) > 0.7:
            continue  

        peptide_str = "".join(peptide)
        if peptide_str in existing_peptides:
            continue
        
        return peptide_str


def generate_peptides(num_peptides=96, peptide_length=5, output_file=desktop_path):
    peptides = set()
    
    while len(peptides) < num_peptides:
        peptide = generate_peptide(peptide_length, peptides)
        peptides.add(peptide)
    
    with open(output_file, "w") as f:
        for seq in peptides:
            f.write(seq + "\n")

    print(f"Generated {num_peptides} unique peptides of length {peptide_length} in ({output_file})")

# change parameters
generate_peptides(num_peptides=96, peptide_length=6)








# to improve / make work
# def filter_synthesizable(sequence):
#     rule1 = False if sequence[-1] == 'W' else True
#     rule2 = ...
#     rule7 = ...
#
#     # return True if all([rule1, rule2, ..., rule7]) else False
#     return sequence if all([rule1, rule2, ..., rule7]) else None
#
#
# # list comprehension
# [filter_synthesizable(seq) for seq in sequences]
#
# # list loop
# valid_seq = []
# for seq in sequences:
#     valid_seq.append(filter_synthesizable(seq))
#
# # df apply
# df['sequences'].apply(filter_synthesizable)

### to improve / make work
# df_syn = pd.DataFrame(names=['sequence'])
# for path in os.walk('./data/seqs/'):
#     df = pd.read_csv(path, names=['sequence'])
#     df_syn = df_syn.append(
#         df['sequence'].apply(filter_synthesizable)
#     )
