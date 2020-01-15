import sys
import math


def create_base_counter():
    base_counter = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}

    return base_counter


def create_dinuc_counter():
    dinuc_counter = {
        "AA": 0, "AC": 0, "AG": 0, "AT": 0,
        "CA": 0, "CC": 0, "CG": 0, "CT": 0,
        "GA": 0, "GC": 0, "GG": 0, "GT": 0,
        "TA": 0, "TC": 0, "TG": 0, "TT": 0,
        "NN": 0}

    return dinuc_counter


def create_aa_dict():
    aa_dict = {"L": "leu",
               "P": "pro",
               "H": "his",
               "Q": "gln",
               "R": "arg",
               "I": "ile",
               "M": "met",
               "T": "thr",
               "N": "asn",
               "K": "lys",
               "S": "ser",
               "V": "val",
               "A": "ala",
               "D": "asp",
               "E": "glu",
               "G": "gly",
               "F": "phe",
               "Y": "tyr",
               "C": "cys",
               "W": "trp",
               "X": "stp",
               "?": "?"}

    return aa_dict


def create_aa_counter():
    aa_counter = {"L": 0,
                  "P": 0,
                  "H": 0,
                  "Q": 0,
                  "R": 0,
                  "I": 0,
                  "M": 0,
                  "T": 0,
                  "N": 0,
                  "K": 0,
                  "S": 0,
                  "V": 0,
                  "A": 0,
                  "D": 0,
                  "E": 0,
                  "G": 0,
                  "F": 0,
                  "Y": 0,
                  "C": 0,
                  "W": 0,
                  "X": 0,
                  "?": 0}

    return aa_counter


def create_codon_dict():
    codon_dict = {"AAA": "K",
                  "AAC": "N",
                  "AAG": "K",
                  "AAT": "N",
                  "ACA": "T",
                  "ACC": "T",
                  "ACG": "T",
                  "ACT": "T",
                  "AGA": "R",
                  "AGC": "S",
                  "AGG": "R",
                  "AGT": "S",
                  "ATA": "I",
                  "ATC": "I",
                  "ATG": "M",
                  "ATT": "I",
                  "CAA": "Q",
                  "CAC": "H",
                  "CAG": "Q",
                  "CAT": "H",
                  "CCA": "P",
                  "CCC": "P",
                  "CCG": "P",
                  "CCT": "P",
                  "CGA": "R",
                  "CGC": "R",
                  "CGG": "R",
                  "CGT": "R",
                  "CTA": "L",
                  "CTC": "L",
                  "CTG": "L",
                  "CTT": "L",
                  "GAA": "E",
                  "GAC": "D",
                  "GAG": "E",
                  "GAT": "D",
                  "GCA": "A",
                  "GCC": "A",
                  "GCG": "A",
                  "GCT": "A",
                  "GGA": "G",
                  "GGC": "G",
                  "GGG": "G",
                  "GGT": "G",
                  "GTA": "V",
                  "GTC": "V",
                  "GTG": "V",
                  "GTT": "V",
                  "TAA": "X",
                  "TAC": "Y",
                  "TAG": "X",
                  "TAT": "Y",
                  "TCA": "S",
                  "TCC": "S",
                  "TCG": "S",
                  "TCT": "S",
                  "TGA": "X",
                  "TGC": "C",
                  "TGG": "W",
                  "TGT": "C",
                  "TTA": "L",
                  "TTC": "F",
                  "TTG": "L",
                  "TTT": "F"}

    return codon_dict


def create_codon_counter():
    codon_counter = {"AAA": 0,
                     "AAC": 0,
                     "AAG": 0,
                     "AAT": 0,
                     "ACA": 0,
                     "ACC": 0,
                     "ACG": 0,
                     "ACT": 0,
                     "AGA": 0,
                     "AGC": 0,
                     "AGG": 0,
                     "AGT": 0,
                     "ATA": 0,
                     "ATC": 0,
                     "ATG": 0,
                     "ATT": 0,
                     "CAA": 0,
                     "CAC": 0,
                     "CAG": 0,
                     "CAT": 0,
                     "CCA": 0,
                     "CCC": 0,
                     "CCG": 0,
                     "CCT": 0,
                     "CGA": 0,
                     "CGC": 0,
                     "CGG": 0,
                     "CGT": 0,
                     "CTA": 0,
                     "CTC": 0,
                     "CTG": 0,
                     "CTT": 0,
                     "GAA": 0,
                     "GAC": 0,
                     "GAG": 0,
                     "GAT": 0,
                     "GCA": 0,
                     "GCC": 0,
                     "GCG": 0,
                     "GCT": 0,
                     "GGA": 0,
                     "GGC": 0,
                     "GGG": 0,
                     "GGT": 0,
                     "GTA": 0,
                     "GTC": 0,
                     "GTG": 0,
                     "GTT": 0,
                     "TAA": 0,
                     "TAC": 0,
                     "TAG": 0,
                     "TAT": 0,
                     "TCA": 0,
                     "TCC": 0,
                     "TCG": 0,
                     "TCT": 0,
                     "TGA": 0,
                     "TGC": 0,
                     "TGG": 0,
                     "TGT": 0,
                     "TTA": 0,
                     "TTC": 0,
                     "TTG": 0,
                     "TTT": 0,
                     "NNN": 0}

    return codon_counter


def old_create_codon_counter():
    codon_counter = {}
    bases = ["A", "C", "G", "T"]
    for base1 in bases:
        for base2 in bases:
            for base3 in bases:
                this_codon = base1 + base2 + base3
                codon_counter[this_codon] = 0

    codon_counter["NNN"] = 0

    return codon_counter


def create_codonpair_counter():
    codonpair_dict = {}
    codon_dict = create_codon_counter()

    for codon in codon_dict:
        if codon == "NNN":
            continue

        for codon2 in codon_dict:
            if codon2 == "NNN":
                continue

            codonpair = codon + codon2
            codonpair_dict[codonpair] = 0

    codonpair_dict["NNNNNN"] = 0

    return codonpair_dict


def create_aapair_counter():
    aanpair_dict = {}
    aa_dict = create_aa_counter()

    for aa in aa_dict:
        if aa == "?":
            continue

        for aa2 in aa_dict:
            if aa2 == "?":
                continue

            aapair = aa + aa2
            aanpair_dict[aapair] = 0

    aanpair_dict["??"] = 0

    return aanpair_dict


def create_aa_syn_dict():
    aa_syn_counter = create_aa_counter()
    codon_dict = create_codon_dict()

    for aa in aa_syn_counter:
        for codon in codon_dict:
            if aa == codon_dict[codon]:
                aa_syn_counter[aa] += 1

    return aa_syn_counter


def correct_seq(this_seq, this_name=""):
    this_seq = this_seq.upper().replace("U", "T").replace("-", "")
    base_dict = ["A", "C", "G", "T", "N"]
    new_seq = ""

    for base in this_seq:
        if base in base_dict:
            new_seq += base
        else:
            new_seq += 'N'
            print("Warning - non ACGTN characters found in sequence [replacing with N]: " + base + " " + this_name)

    return new_seq


def count_bases(this_seq, base_dict):

    for base in base_dict:
        base_dict[base] += this_seq.count(base)

    return base_dict


def count_bridge_bases(this_seq, base_dict):

    for i in range(2, len(this_seq)-1, 3):
        base_dict[this_seq[i]] += 1

    for i in range(3, len(this_seq), 3):
        base_dict[this_seq[i]] += 1

    return base_dict


def count_dinucleotides(this_seq, dinuc_dict):
    for i in range(len(this_seq)-1):
        this_dinuc = this_seq[i:i+2]

        if this_dinuc in dinuc_dict:
            dinuc_dict[this_dinuc] += 1
        else:
            dinuc_dict["NN"] += 1

    return dinuc_dict


def count_bridge_dinucleotides(this_seq, dinuc_dict):
    for i in range(2, len(this_seq)-1, 3):
        this_dinuc = this_seq[i:i+2]

        if this_dinuc in dinuc_dict:
            dinuc_dict[this_dinuc] += 1
        else:
            dinuc_dict["NN"] += 1

    return dinuc_dict


def count_nonbridge_dinucleotides(this_seq, dinuc_dict):
    for i in range(0, len(this_seq)-1, 3):
        this_dinuc = this_seq[i:i+2]

        if this_dinuc in dinuc_dict:
            dinuc_dict[this_dinuc] += 1
        else:
            dinuc_dict["NN"] += 1

    for i in range(1, len(this_seq)-1, 3):
        this_dinuc = this_seq[i:i+2]

        if this_dinuc in dinuc_dict:
            dinuc_dict[this_dinuc] += 1
        else:
            dinuc_dict["NN"] += 1

    return dinuc_dict


def count_codons(this_seq, codon_dict):
    for i in range(0, len(this_seq)-1, 3):
        this_codon = this_seq[i:i+3]

        if this_codon in codon_dict:
            codon_dict[this_codon] += 1
        else:
            codon_dict["NNN"] += 1

    return codon_dict


def count_aas(this_seq, aa_dict, codon_dict):
    for i in range(0, len(this_seq)-1, 3):
        this_codon = this_seq[i:i+3]

        if this_codon in codon_dict:
            this_aa = codon_dict[this_codon]
        else:
            this_aa = "?"

        aa_dict[this_aa] += 1

    return aa_dict


def count_thru_stops(this_seq, codon_dict):
    this_stops = 0

    for i in range(0, len(this_seq) - 1, 3):
        this_codon = this_seq[i:i + 3]

        if this_codon in codon_dict:
            this_aa = codon_dict[this_codon]
        else:
            this_aa = "?"

        if i < len(this_seq) - 3 and this_aa == "X":
            this_stops += 1

    return this_stops


def count_codonpairs(this_seq, codonpair_dict):
    for i in range(0, len(this_seq)-4, 3):
        this_codonpair = this_seq[i:i+6]

        if this_codonpair in codonpair_dict:
            codonpair_dict[this_codonpair] += 1
        else:
            codonpair_dict["NNNNNN"] += 1

    return codonpair_dict


def count_aapairs(this_seq, aapair_dict, codon_dict):
    for i in range(0, len(this_seq)-4, 3):
        this_codon = this_seq[i:i+3]
        this_codon2 = this_seq[i+3:i+6]

        if this_codon in codon_dict:
            this_aa = codon_dict[this_codon]
        else:
            this_aa = "?"

        if this_codon2 in codon_dict:
            this_aa2 = codon_dict[this_codon2]
        else:
            this_aa2 = "?"

        this_aapair = this_aa + this_aa2

        if this_aapair in aapair_dict:
            aapair_dict[this_aapair] += 1
        else:
            aapair_dict["??"] += 1

    return aapair_dict


def count_lines(this_filename):
    with open(this_filename) as f:
        for i, l in enumerate(f):
            pass

    return i + 1


def check_fasta(this_line):
    if this_line.rfind(">") > 0:
        print("Error - FASTA sequence start symbol '>' found mid sequence " + this_line)


def check_seq_coding(this_seq, this_name):
    seq_test = False

    if len(this_seq) % 3 == 0:
        seq_test = True
    else:
        print("Sequence not divisible by 3: " + this_name)

    return seq_test


def trim_seq(this_seq):
    new_seq = ""

    for i in range(0, len(this_seq) - 3, 3):
        new_seq += this_seq[i:i+3]

    return new_seq


def get_filename_stub(this_name, this_ext):
    this_stub = this_name

    this_pos = this_name.rfind(this_ext)
    if this_pos > 0:
        this_stub = this_name[:this_pos]

    return this_stub


def calculate_base_freqs(base_counts):
    base_sum = sum(base_counts.values())
    base_freqs = {}
    for base in base_counts:
        base_freqs[base] = (base_counts[base]/float(base_sum))

    return base_freqs


def calculate_gc_content(base_counts):
    this_total = sum(base_counts.values())
    this_content = {"GC": (base_counts["G"] + base_counts["C"])/(float(this_total)),
                    "AT": (base_counts["A"] + base_counts["T"])/(float(this_total))}

    return this_content


def calculate_aa_freqs(aa_counts):
    aa_sum = sum(aa_counts.values())
    aa_freqs = {}
    for aa in aa_counts:
        if aa == "?":
            continue

        aa_freqs[aa] = (aa_counts[aa]/float(aa_sum))

    return aa_freqs


def calculate_aapair_freqs(aapair_counts):
    aapir_sum = sum(aapair_counts.values())
    aapair_freqs = {}
    for aapair in aapair_counts:
        aapair_freqs[aapair] = (aapair_counts[aapair]/float(aapir_sum))

    # obs/exp accounting for count of AA1 and AA2 instead of freqs - or in addition to
    return aapair_freqs


def calculate_codon_freqs(codon_counts, aa_counts, codon_dict):
    codon_freqs = {}

    for codon in codon_counts:
        if codon == "NNN":
            continue

        this_aa = codon_dict[codon]
        this_aa_count = aa_counts[this_aa]

        if this_aa_count > 0:
            codon_freqs[codon] = (codon_counts[codon]/float(this_aa_count))
        else:
            codon_freqs[codon] = 0

    return codon_freqs


def calculate_dinuc_freqs(dinuc_counts, base_counts):
    dinuc_sum = sum(dinuc_counts.values())
    base_sum = sum(base_counts.values())

    dinuc_freqs = {}

    for dinuc in dinuc_counts:
        if dinuc == "NN":
            continue

        base1 = dinuc[0:1]
        base2 = dinuc[1:2]

        if base_counts[base1] == 0 or base_counts[base2] == 0:
            dinuc_freqs[dinuc] = 0
        else:
            dinuc_freqs[dinuc] = (dinuc_counts[dinuc]/float(dinuc_sum))/((base_counts[base1]/float(base_sum)) * (base_counts[base2]/float(base_sum)))

    return dinuc_freqs


def calculate_raw_dinuc_freqs(dinuc_counts):
    dinuc_sum = sum(dinuc_counts.values())

    dinuc_freqs = {}

    for dinuc in dinuc_counts:
        if dinuc == "NN":
            continue

        dinuc_freqs[dinuc] = (dinuc_counts[dinuc]/float(dinuc_sum))

    return dinuc_freqs


def calculate_cpb_scores(codonpair_counts, codon_counts, aapair_counts, aa_counts, codon_aa_dict):
    cpbs = {}

    for codonpair in codonpair_counts:
        if codonpair == "NNNNNN":
            continue

        codon1 = codonpair[0:3]
        codon2 = codonpair[3:6]

        if codonpair_counts[codonpair] > 0:
            aa1 = codon_aa_dict[codon1]
            aa2 = codon_aa_dict[codon2]
            aapair = aa1 + aa2

            cpbs[codonpair] = math.log((float(codonpair_counts[codonpair]))/float(((codon_counts[codon1] * codon_counts[codon2])/(aa_counts[aa1] * aa_counts[aa2]))*aapair_counts[aapair]))
        else:
            cpbs[codonpair] = 0

    return cpbs


def fasta_calculate_cpbs(filename, force_trim):
    output_filename = get_filename_stub(filename, ".f") + "_seqfeat.txt"
    print("Input FASTA file = " + filename)
    print("Output feature file = " + output_filename)

    last_line = count_lines(filename)

    line_count = 0
    seq_count = 0
    seq_incomplete = 0
    seq_added = 0
    seq_zero = 0
    all_seqs = []

    with open(filename) as file_handler:
        seq = ''
        name = ''

        for line in file_handler:
            line = line.rstrip()
            line_count += 1
            check_fasta(line)

            if line.find(">") != 0:
                seq += line.upper()

            if line.find(">") == 0 or line_count == last_line:
                if seq_count > 0:
                    if name.find(" ") > 0:
                        new_name = name[:name.find(" ")]
                    else:
                        new_name = name

                    seq = correct_seq(seq, name)
                    this_check = check_seq_coding(seq, name)

                    if not this_check:
                        seq_incomplete += 1

                        if force_trim:
                            seq = trim_seq(seq)

                    if len(seq) == 0:
                        seq_zero += 1
                    else:
                        all_seqs.append([seq_count, new_name, name, seq])
                        seq_added += 1

                    name = ''
                    seq = ''

                if line.find(">") == 0:
                    name = line
                    seq_count += 1

    print("Coding sequences  = " + str(seq_count))
    print("Length 0 sequences = " + str(seq_zero))
    print("Incomplete sequences [not divisible by 3]  = " + str(seq_incomplete))
    print("Usable sequences = " + str(seq_added))

    # Create the initial counters - for the header
    base_counts = create_base_counter()
    base_bridge_counts = create_base_counter()
    dinuc_counts = create_dinuc_counter()
    dinuc_bridge_counts = create_dinuc_counter()
    dinuc_nonbridge_counts = create_dinuc_counter()
    codon_counts = create_codon_counter()
    codonpair_counts = create_codonpair_counter()
    aa_counts = create_aa_counter()
    aapair_counts = create_aapair_counter()
    # Create code dict
    codon_aa_dict = create_codon_dict()

    with open(output_filename, "w") as file_output:

        # Header
        file_output.write("SeqName\tSeqs\tSeqLength\tCodons\tBadCodons\tCodonPairs\tStops\tA\tC\tG\tT\tN\tGC\tAT\tCpG\tUpA")

        for aa in aa_counts:
            if aa != "?":
                file_output.write("\t" + aa + "-Bias")

        for codon in codon_counts:
            if codon != "NNN":
                file_output.write("\t" + codon + "-Bias")

        for codonpair in codonpair_counts:
            if codonpair != "NNNNNN":
                file_output.write("\t" + codonpair[0:3] + "-" + codonpair[3:6])

        for dinuc in dinuc_counts:
            if dinuc != "NN":
                if dinuc == "CG" or dinuc == "TA":
                    continue

                base1 = dinuc[0].replace("T", "U")
                base2 = dinuc[1].replace("T", "U")
                file_output.write("\t" + base1 + "p" + base2)

        for dinuc in dinuc_bridge_counts:
            if dinuc != "NN":
                base1 = dinuc[0].replace("T", "U")
                base2 = dinuc[1].replace("T", "U")
                file_output.write("\tbr" + base1 + "p" + base2)

        for dinuc in dinuc_nonbridge_counts:
            if dinuc != "NN":
                base1 = dinuc[0].replace("T", "U")
                base2 = dinuc[1].replace("T", "U")
                file_output.write("\tNonBr" + base1 + "p" + base2)

        for base in base_counts:
            file_output.write("\tBr-" + base)

        for aapair in aapair_counts:
            if aapair != "??":
                file_output.write("\t" + aapair[0:1] + "-" + aapair[1:2])

        file_output.write("\n")
        # End of Header

        for seq in all_seqs:
            if seq[0] % 1000 == 0:
                print("Analysing sequence = " + str(seq[0]))

            base_counts = count_bases(seq[3], base_counts)
            base_bridge_counts = count_bridge_bases(seq[3], base_bridge_counts)

            dinuc_counts = count_dinucleotides(seq[3], dinuc_counts)
            dinuc_bridge_counts = count_bridge_dinucleotides(seq[3], dinuc_bridge_counts)
            dinuc_nonbridge_counts = count_nonbridge_dinucleotides(seq[3], dinuc_nonbridge_counts)

            codon_counts = count_codons(seq[3], codon_counts)
            codonpair_counts = count_codonpairs(seq[3], codonpair_counts)
            aa_counts = count_aas(seq[3], aa_counts, codon_aa_dict)
            thru_stops = count_thru_stops(seq[3], codon_aa_dict)
            aapair_counts = count_aapairs(seq[3], aapair_counts, codon_aa_dict)
            cpb_scores = calculate_cpb_scores(codonpair_counts, codon_counts, aapair_counts, aa_counts, codon_aa_dict)

            total_bases = sum(base_counts.values())
            total_codons = sum(codon_counts.values())
            total_codonpairs = sum(codonpair_counts.values())

            base_freqs = calculate_base_freqs(base_counts)
            bridge_base_freqs = calculate_base_freqs(base_bridge_counts)
            gc_content = calculate_gc_content(base_counts)
            aa_freqs = calculate_aa_freqs(aa_counts)
            aapair_freqs = calculate_aapair_freqs(aapair_counts)
            codon_freqs = calculate_codon_freqs(codon_counts, aa_counts, codon_aa_dict)
            dinuc_freqs = calculate_dinuc_freqs(dinuc_counts, base_counts)
            dinuc_bridge_freqs = calculate_dinuc_freqs(dinuc_bridge_counts, base_bridge_counts)
            dinuc_nonbridge_freqs = calculate_dinuc_freqs(dinuc_nonbridge_counts, base_counts)

            file_output.write(seq[1] + "\t" + str(seq[0]))
            file_output.write("\t" + str(total_bases) + "\t" + str(total_codons) + "\t" + str(codon_counts["NNN"]) + "\t" + str(total_codonpairs) + "\t" + str(thru_stops))

            for base in base_freqs:
                file_output.write("\t" + str(base_freqs[base]))

            file_output.write("\t" + str(gc_content["GC"]) + "\t" + str(gc_content["AT"]))
            file_output.write("\t" + str(dinuc_freqs["CG"]) + "\t" + str(dinuc_freqs["TA"]))

            for aa in aa_freqs:
                file_output.write("\t" + str(aa_freqs[aa]))

            for codon in codon_freqs:
                file_output.write("\t" + str(codon_freqs[codon]))

            # given dinuc is now explanation of codon pair - what about plain old cpb freqs
            for cpb in cpb_scores:
                file_output.write("\t" + str(cpb_scores[cpb]))

            # raw instead/as well of obs/exp?: calculate_raw_dinuc_freqs()
            for dinuc in dinuc_freqs:
                if dinuc == "CG" or dinuc == "TA":
                    continue
                file_output.write("\t" + str(dinuc_freqs[dinuc]))

            for dinuc in dinuc_bridge_freqs:
                file_output.write("\t" + str(dinuc_bridge_freqs[dinuc]))

            for dinuc in dinuc_nonbridge_freqs:
                file_output.write("\t" + str(dinuc_nonbridge_freqs[dinuc]))

            for base in bridge_base_freqs:
                file_output.write("\t" + str(bridge_base_freqs[base]))

            # obs/exp instead/as well?
            for aa in aapair_freqs:
                file_output.write("\t" + str(aapair_freqs[aa]))

            file_output.write("\n")

            base_counts = create_base_counter()
            base_bridge_counts = create_base_counter()
            dinuc_counts = create_dinuc_counter()
            dinuc_bridge_counts = create_dinuc_counter()
            dinuc_nonbridge_counts = create_dinuc_counter()
            codon_counts = create_codon_counter()
            codonpair_counts = create_codonpair_counter()
            aa_counts = create_aa_counter()
            aapair_counts = create_aapair_counter()


print("seqfeats_host_cpb.py started...")

arguments = len(sys.argv)

if arguments != 2:
    print("Error - incorrect number of arguments - example usage:")
    print("seqfeats_host_cpb.py coding_sequences.fasta")
    sys.exit(1)
else:
    fasta_calculate_cpbs(sys.argv[1], True)

print("\n...finished SeqFeat.py")
