import sys
from math import log, exp


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
               "X": "stp"}

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
                  "X": 0}
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
                     "TTT": 0}

    return codon_counter


def create_aa_syn_dict():
    aa_syn_counter = create_aa_counter()
    codon_dict = create_codon_dict()

    for aa in aa_syn_counter:
        for codon in codon_dict:
            if aa == codon_dict[codon]:
                aa_syn_counter[aa] += 1

    return aa_syn_counter


def check_fasta(line):
    if line.rfind(">") > 0:
        print("Error - FASTA sequence start symbol '>' found mid sequence " + line)


def check_seq_coding(seq, name):
    seq_test = False

    if len(seq) % 3 == 0:
        seq_test = True
    else:
        print("Sequence not divisible by 3: " + name)

    return seq_test


def correct_seq(seq):
    seq = replace_ambiguities(seq)
    seq = rna_to_dna(seq)
    seq = remove_gaps(seq)

    return seq


def replace_ambiguities(seq):
    seq = seq.upper()
    base_dict = ["A", "C", "G", "T", "U", "N", "-"]
    no_amb_seq = ''

    for base in seq:
        if base in base_dict:
            no_amb_seq += base
        else:
            no_amb_seq += 'N'

    return no_amb_seq


def rna_to_dna(seq):

    return seq.upper().replace("U", "T")


def remove_gaps(seq):

    return seq.upper().replace("-", "")


def count_lines(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass

    return i + 1


def get_filename_stub(name, ext):
    stub = name

    pos = name.rfind(ext)
    if pos > 0:
        stub = name[:pos]

    return stub


def fasta_calculate_cai(filename, allow_ns):
    print("fasta_calculate_cai")

    output_filename_table = get_filename_stub(filename, ".f") + "_cu_table.txt"
    output_filename_cai = get_filename_stub(filename, ".f") + "_cai.txt"
    print("Input FASTA file = " + filename)
    print("Output CU table file = " + output_filename_table)
    print("Output CAI data file = " + output_filename_cai)
    print("Sequences with Ns allowed [1 = True, 0 = False] = " + str(allow_ns))

    last_line = count_lines(filename)
    line_count = 0
    seq_count = 0
    proc_seq_count = 0

    ref_codons = create_codon_dict()
    ref_codons_counts = create_codon_counter()
    ref_codons_weights = create_codon_counter()
    ref_codons_freqs = create_codon_counter()
    ref_codons_thousand = create_codon_counter()
    ref_codons_rcsu = create_codon_counter()

    ref_aa_counts = create_aa_counter()
    ref_aa_syn = create_aa_syn_dict()
    ref_aa_max_syn = create_aa_counter()
    ref_aa_max_rcsu = create_aa_counter()

    with open(filename) as file_handler:
        seq = ''
        name = ''

        for line in file_handler:
            check_fasta(line)

            line = line.rstrip()
            line_count += 1

            if line.find(">") != 0:
                seq += line.upper()

            if line.find(">") == 0 or line_count == last_line:
                if seq_count > 0:
                    seq = correct_seq(seq)

                    if allow_ns == 1 or seq.count("N") == 0:
                        update_refseq_cai(seq, name, ref_codons_counts, ref_codons, ref_aa_counts)
                        proc_seq_count += 1

                    name = ''
                    seq = ''

                if line.find(">") == 0:
                    name = line
                    seq_count += 1

    print("Calculated ref codon table: processed sequences = " + str(proc_seq_count)+", total sequences = " + str(seq_count))

    calculate_ref_cai(ref_codons, ref_codons_counts, ref_codons_weights, ref_codons_thousand, ref_codons_freqs, ref_codons_rcsu, ref_aa_counts, ref_aa_syn, ref_aa_max_syn, ref_aa_max_rcsu)

    with open(output_filename_table, "w") as file_output:
        file_output.write("Codon\tCount\tWeight\tCodonBias\tNumberPer1000\tRCSU\n")

        for codon in ref_codons_counts:
            file_output.write(codon + "\t" + str(ref_codons_counts[codon]) + "\t" + str(ref_codons_weights[codon]) + "\t" + str(ref_codons_freqs[codon]) + "\t" + str(ref_codons_thousand[codon]) + "\t" + str(ref_codons_rcsu[codon]) + "\n")

        file_output.write("\n\n")
        file_output.write("Codon\tNumberPer1000\t(Count)\tCodon\tNumberPer1000\t(Count)\tCodon\tNumberPer1000\t(Count)\tCodon\tNumberPer1000\t(Count)\n")

        bases = ["T", "C", "A", "G"]

        for pos1 in bases:
            for pos3 in bases:
                counter = 0
                for pos2 in bases:
                    this_codon = pos1 + pos2 + pos3
                    uracil_codon = this_codon.replace("T", "U")

                    if counter > 0:
                        file_output.write("\t")

                    file_output.write(uracil_codon + "\t" + str(ref_codons_thousand[this_codon]) + "\t(" + str(ref_codons_counts[this_codon]) + ")")
                    counter += 1

                file_output.write("\n")

            file_output.write("\n")

    print("Outputted CU table")

    seq_count = 0
    proc_seq_count = 0
    line_count = 0

    with open(output_filename_cai, "w") as file_output:
        file_output.write("SeqName\tCAI\n")

        with open(filename) as file_handler:
            seq = ''
            name = ''

            for line in file_handler:
                line = line.rstrip()
                line_count += 1

                if line.find(">") != 0:
                    seq += line.upper()

                if line.find(">") == 0 or line_count == last_line:
                    if seq_count > 0:
                        seq = correct_seq(seq)
                        if allow_ns == 1 or seq.count("N") == 0:
                            this_cai = calculate_seq_cai(seq, name, ref_codons, ref_codons_weights)
                            file_output.write(name + "\t" + str(this_cai) + "\n")
                            proc_seq_count += 1

                        name = ''
                        seq = ''

                    if line.find(">") == 0:
                        name = line
                        seq_count += 1

                check_fasta(line)

        print("Calculated CAI for each sequence: processed sequences = " + str(proc_seq_count)+", total sequences = " + str(seq_count))


def update_refseq_cai(seq, name, ref_codons_counts, ref_codons, ref_aa_counts):
    if check_seq_coding(seq, name):
        i = 0
        while i < len(seq):
            this_codon = seq[i] + seq[i + 1] + seq[i + 2]
            ns = this_codon.count("N")

            if ns == 0:
                ref_codons_counts[this_codon] += 1
                aa = ref_codons[this_codon]
                ref_aa_counts[aa] += 1

            i += 3


def calculate_ref_cai(ref_codons, ref_codons_counts, ref_codons_weights, ref_codons_thousand, ref_codons_freqs, ref_codons_rcsu, ref_aa_counts, ref_aa_syn, ref_aa_max_syn, ref_aa_max_rcsu):
    total_codons = 0

    for codon in ref_codons_counts:
        total_codons += ref_codons_counts[codon]

    for codon in ref_codons_counts:
        ref_codons_thousand[codon] = ref_codons_counts[codon] / total_codons * 1000

    for aa in ref_aa_max_syn:
        this_max = 0

        for codon in ref_codons:
            if ref_codons[codon] == aa:
                if ref_codons_counts[codon] > this_max:
                    this_max = ref_codons_counts[codon]

        ref_aa_max_syn[aa] = this_max

    for aa in ref_aa_max_syn:
        for codon in ref_codons:
            if ref_codons[codon] == aa:
                if ref_codons_counts[codon] > 0:
                    ref_codons_weights[codon] = ref_codons_counts[codon] / ref_aa_max_syn[aa]
                else:
                    ref_codons_weights[codon] = 0

    for codon in ref_codons_freqs:
        if ref_aa_counts[ref_codons[codon]] > 0:
            ref_codons_freqs[codon] = ref_codons_counts[codon] / ref_aa_counts[ref_codons[codon]]

        elif ref_codons_counts[codon] > 0:
            print("Error - calculate_ref_cai - codon count > 0 but AA count = 0: "+codon)

    for codon in ref_codons_counts:
        aa = ref_codons[codon]
        if ref_aa_counts[aa] > 0:
            ref_codons_rcsu[codon] = ref_codons_counts[codon] / (ref_aa_counts[aa] * (1 / ref_aa_syn[aa]))

    for aa in ref_aa_max_rcsu:
        this_max = 0

        for codon in ref_codons:
            if ref_codons[codon] == aa:
                if ref_codons_rcsu[codon] > this_max:
                    this_max = ref_codons_rcsu[codon]

        ref_aa_max_rcsu = this_max


def calculate_seq_cai(seq, name, ref_codons, ref_codons_weights):
    weight = 0
    this_codon_count = 0

    # check divisible by 3
    if check_seq_coding(seq, name):
        i = 0
        while i < len(seq):
            this_codon = seq[i] + seq[i + 1] + seq[i + 2]

            ns = this_codon.count("N")

            if ns == 0:
                this_aa = ref_codons[this_codon]
                if this_aa != "M" and this_aa != "W" and this_aa != "X":
                    weight += log(ref_codons_weights[this_codon])
                    this_codon_count += 1

            i += 3

    if this_codon_count == 0:
        print("Warning - calculate_seq_cai - sequence has 0 usable codons: " + name + " " + seq)
    else:
        weight = exp(weight / this_codon_count)

    return weight


print("seqfeats_cai started...")

example_messages = ["seqfeats_cai.py fasta_file\nseqfeats_cai.py fasta_file allow_Ns\nallowNs: allow sequences to have Ns: TRUE by default"]

error_message = "Error - incorrect number of arguments - example usage:"

arguments = len(sys.argv)

if arguments == 1:
    print(error_message)

    for example in example_messages:
        print(example)

else:
    if arguments != 2 and arguments != 3:
        print(error_message + "\n" + example_messages[0])
    else:
        arg_ns = 1

        if arguments == 3:
            if sys.argv[2].lower() == "yes" or sys.argv[2].lower() == "y" or sys.argv[2].lower() == "true" or sys.argv[2].lower() == "1":
                arg_ns = 1
            elif sys.argv[2].lower() == "no" or sys.argv[2].lower() == "n" or sys.argv[2].lower() == "false" or sys.argv[2].lower() == "0":
                arg_ns = 0
            else:
                print("Error - unrecognised true/false argument can have [yes/y/true/1 or no/n/false/0]: "+sys.argv[2])
                sys.exit(1)

        fasta_calculate_cai(sys.argv[1], arg_ns)

print("...finished seqfeats_cai")
