codon_table = {
    'AUG': ('Met', 'M'), 'UUU': ('Phe', 'F'), 'UUC': ('Phe', 'F'),
    'UUA': ('Leu', 'L'), 'UUG': ('Leu', 'L'), 'UCU': ('Ser', 'S'), 'UCC': ('Ser', 'S'),
    'UCA': ('Ser', 'S'), 'UCG': ('Ser', 'S'), 'UAU': ('Tyr', 'Y'), 'UAC': ('Tyr', 'Y'),
    'UGU': ('Cys', 'C'), 'UGC': ('Cys', 'C'), 'UGG': ('Trp', 'W'), 'CUU': ('Leu', 'L'),
    'CUC': ('Leu', 'L'), 'CUA': ('Leu', 'L'), 'CUG': ('Leu', 'L'), 'CCU': ('Pro', 'P'),
    'CCC': ('Pro', 'P'), 'CCA': ('Pro', 'P'), 'CCG': ('Pro', 'P'), 'CAU': ('His', 'H'),
    'CAC': ('His', 'H'), 'CAA': ('Gln', 'Q'), 'CAG': ('Gln', 'Q'), 'CGU': ('Arg', 'R'),
    'CGC': ('Arg', 'R'), 'CGA': ('Arg', 'R'), 'CGG': ('Arg', 'R'), 'AUU': ('Ile', 'I'),
    'AUC': ('Ile', 'I'), 'AUA': ('Ile', 'I'), 'GUU': ('Val', 'V'), 'GUC': ('Val', 'V'),
    'GUA': ('Val', 'V'), 'GUG': ('Val', 'V'), 'GCU': ('Ala', 'A'), 'GCC': ('Ala', 'A'),
    'GCA': ('Ala', 'A'), 'GCG': ('Ala', 'A'), 'GAU': ('Asp', 'D'), 'GAC': ('Asp', 'D'),
    'GAA': ('Glu', 'E'), 'GAG': ('Glu', 'E'), 'GGU': ('Gly', 'G'), 'GGC': ('Gly', 'G'),
    'GGA': ('Gly', 'G'), 'GGG': ('Gly', 'G'), 'UAA': ('Stop', None), 'UAG': ('Stop', None),
    'UGA': ('Stop', None), 'ACU': ('Thr', 'T'), 'ACC': ('Thr', 'T'), 'ACA': ('Thr', 'T'),
    'ACG': ('Thr', 'T'), 'AAU': ('Asn', 'N'), 'AAC': ('Asn', 'N'), 'AAA': ('Lys', 'K'),
    'AAG': ('Lys', 'K'), 'AGU': ('Ser', 'S'), 'AGC': ('Ser', 'S'), 'AGA': ('Arg', 'R'),
    'AGG': ('Arg', 'R')
}

possibilities = {'Stop': ['UAA', 'UAG', 'UGA'], 'M': ['AUG'], 'F': ['UUU', 'UUC'],
                 'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'], 'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
                 'Y': ['UAU', 'UAC'], 'C': ['UGU', 'UGC'], 'W': ['UGG'], 'P': ['CCU', 'CCC', 'CCA', 'CCG'],
                 'H': ['CAU', 'CAC'], 'Q': ['CAA', 'CAG'], 'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
                 'I': ['AUU', 'AUC', 'AUA'], 'V': ['GUU', 'GUC', 'GUA', 'GUG'], 'A': ['GCU', 'GCC', 'GCA', 'GCG'],
                 'D': ['GAU', 'GAC'], 'E': ['GAA', 'GAG'], 'G': ['GGU', 'GGC', 'GGA', 'GGG'],
                 'T': ['ACU', 'ACC', 'ACA', 'ACG'], 'N': ['AAU', 'AAC'], 'K': ['AAA', 'AAG']}


def isValidDNA(dna):
    if len(dna) % 3 != 0:
        return False
    for i in dna:
        if i not in "ATCG":
            return False
    return True


def complement(dna):
    comp = ""
    for i in dna:
        if i not in "ATGC":
            return "Invalid."
        if i == "A":
            comp += "T"
        if i == "T":
            comp += "A"
        if i == "G":
            comp += "C"
        if i == "C":
            comp += "G"
    return comp


def toMrna(dna):
    comp = complement(dna)
    if comp != "Invalid.":
        return "".join(["U" if i == "T" else i for i in comp])
    else:
        return "Invalid."


def toAminoAcid(dna):
    mrna = toMrna(dna)
    aminos = []
    if mrna != "Invalid.":
        for i in range(0, len(mrna), 3):
            codon = codon_table[mrna[i: i + 3]]
            if codon[0] == "Stop":
                print("Stop Codon detected")
                break
            else:
                aminos.append(f"{codon[0]} ({codon[1]})")
        return " - ".join(aminos)
    else:
        return "Invalid."


def showStepsToAminoAcid(dna):
    print(f"Complement = {complement(dna)}")
    print(f"mRNA = {toMrna(dna)}")
    print(f"Amino acid = {toAminoAcid(dna)}")


def isValidAminoAcid(codons):
    for i in codons:
        if i not in possibilities:
            return False
    return True


def countRNACodons(rna):
    counts = {}
    for i in range(0, len(rna), 3):
        rnaCodon = rna[i:i+3]
        if rnaCodon in counts:
            counts[rnaCodon] += 1
        else:
            counts[rnaCodon] = 1
    print(f"\nmRNA = {rna}")
    for i in counts:
        print(f"{i} = {counts[i]}")


def recurse(codons, rna=""):
    if len(codons) == 0:
        countRNACodons(rna)
    else:
        for i in possibilities[codons[0]]:
            recurse(codons[1:], rna + i)


while True:
    choice = input("\nChoose mode: 1) DNA to Amino Acid, 2) Amino Acid to mRNA: ").strip()
    
    if choice == '1':
        dna = input("Input DNA = ").strip()
        while not isValidDNA(dna):
            if len(dna) % 3 != 0:
                print("Length has to be a multiple of 3")
            else:
                print("Invalid Characters Detected")
            dna = input("Input DNA = ").strip()

        print()
        showStepsToAminoAcid(dna)

    elif choice == '2':
        codons = input("Input Amino acid sequence (Max 3 characters) = ").strip().upper()
        while len(codons) > 3 or not isValidAminoAcid(codons):
            if len(codons) > 3:
                print("Maximum length is 3")
            else:
                print("Invalid codon(s) detected")
            codons = input("Input Amino acid sequence (Max 3 characters) = ").strip().upper()

        recurse(codons)
    
    else:
        print("Invalid option. Please choose 1 or 2.")
