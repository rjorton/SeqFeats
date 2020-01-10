import sys


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


def correct_seq(this_seq, this_name):
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


def count_bases(this_seq):
    base_dict = create_base_counter()

    for base in base_dict:
        base_dict[base] += this_seq.count(base)

    return base_dict


def count_dinucleotides(this_seq):
    dinuc_dict = create_dinuc_counter()

    for i in range(len(this_seq)-1):
        this_dinuc = this_seq[i:i+2]

        if this_dinuc in dinuc_dict:
            dinuc_dict[this_dinuc] += 1
        else:
            dinuc_dict["NN"] += 1

    return dinuc_dict


def count_lines(this_filename):
    with open(this_filename) as f:
        for i, l in enumerate(f):
            pass

    return i + 1


def check_fasta(this_line):
    if this_line.rfind(">") > 0:
        print("Error - FASTA sequence start symbol '>' found mid sequence " + this_line)


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


def calculate_dinuc_obsexp(dinuc_counts, base_counts):
    dinuc_sum = sum(dinuc_counts.values())
    base_sum = sum(base_counts.values())

    dinuc_obsexp = {}

    for dinuc in dinuc_counts:
        if dinuc == "NN":
            continue

        if dinuc_counts[dinuc] == 0:
            dinuc_obsexp[dinuc] = 0
        else:
            base1 = dinuc[0:1]
            base2 = dinuc[1:2]

            dinuc_obsexp[dinuc] = (dinuc_counts[dinuc]/float(dinuc_sum))/((base_counts[base1]/float(base_sum)) * (base_counts[base2]/float(base_sum)))

    return dinuc_obsexp


def calculate_dinuc_freqs(dinuc_counts):
    dinuc_sum = sum(dinuc_counts.values())

    dinuc_freqs = {}

    for dinuc in dinuc_counts:
        dinuc_freqs[dinuc] = (dinuc_counts[dinuc]/float(dinuc_sum))

    return dinuc_freqs


def fasta_calculate_dinucs(filename, remove_duplicates):
    print("fasta_calculate_dinucs")
    output_filename = get_filename_stub(filename, ".f") + "_dat.txt"
    print("Input FASTA file = " + filename)
    print("Output data file = " + output_filename)
    print("Removing duplicate genes [preserving longest length - first instance] = " + str(remove_duplicates))

    last_line = count_lines(filename)
    line_count = 0
    seq_count = 0
    seq_unavailable = 0
    seq_zero = 0
    seq_added = 0
    all_seqs = []

    with open(filename) as file_handler:
        name = ''
        sequence = ''

        for line in file_handler:
            line = line.rstrip()
            check_fasta(line)
            line_count += 1

            if line.find(">") != 0:
                sequence += line.upper()

            if line.find(">") == 0 or line_count == last_line:
                if seq_count > 0:
                    # duplicate removal is based on gene name
                    # could still be used to remove duplicates based on whole name if not in ensembl format
                    gene = name
                    transcript = ''

                    pos = name.find(" gene:")
                    if pos > 0:
                        # >ENST00000311787.5 cds chromosome:GRCh38:10:49131154:49134008:-1 gene:ENSG00000172538.6 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:FAM170B description:family with sequence similarity 170 member B [Source:HGNC Symbol;Acc:HGNC:19736]
                        transcript = name[1:name.find(" ")]
                        gene = name[name.find(" gene:") + 6:]
                        pos = gene.find(" ")
                        if pos <= 0:
                            gene = gene[0:]
                        else:
                            gene = gene[0:pos]
                        name = ">" + gene + "|" + transcript

                    if sequence.lower() == "sequence unavailable":
                        seq_unavailable += 1
                    elif len(sequence) == 0:
                        seq_zero += 1
                    else:
                        sequence = correct_seq(sequence, name)
                        all_seqs.append([seq_count, name, gene, transcript, sequence, len(sequence)])
                        seq_added += 1

                    name = ''
                    sequence = ''

                if line.find(">") == 0:
                    name = line
                    seq_count += 1

    print("Total sequences in file  = " + str(seq_count))
    print("Usable sequences  = " + str(seq_added))
    print("Zero length sequences = " + str(seq_zero))
    print("Sequence unavailable = " + str(seq_unavailable))

    if remove_duplicates:
        seq_dict = {}
        for seq in all_seqs:
            if seq[2] in seq_dict:
                this_seq = seq_dict[seq[2]]
                if seq[5] > this_seq[5]:
                    seq_dict[seq[2]] = seq
            else:
                seq_dict[seq[2]] = seq

        all_seqs = seq_dict.values()

        print("Sequences after duplicate removal = " + str(len(all_seqs)))

    with open(output_filename, "w") as file_output:

        # Header - should the T single base be a U - and ATcontent AU - AT is standard term
        file_output.write("Num\tName\tGene\tTranscript\tLength\tA\tC\tG\tT\tN")
        file_output.write("\tGCcontent\tATcontent")
        file_output.write("\tApA\tApC\tApG\tApU\tCpA\tCpC\tCpG\tCpU\tGpA\tGpC\tGpG\tGpU\tUpA\tUpC\tUpG\tUpU")
        file_output.write("\tAA%\tAC%\tAG%\tAU%\tCA%\tCC%\tCG%\tCU%\tGA%\tGC%\tGG%\tGU%\tUA%\tUC%\tUG%\tUU%")
        file_output.write("\n")

        seq_count = 0
        for seq in all_seqs:
            seq_count += 1

            if seq_count % 10000 == 0:
                print("Analysing sequence number = " + str(seq_count))

            base_counts = count_bases(seq[4])
            base_freqs = calculate_base_freqs(base_counts)
            gc_content = calculate_gc_content(base_counts)
            dinuc_counts = count_dinucleotides(seq[4])
            dinuc_obsexp = calculate_dinuc_obsexp(dinuc_counts, base_counts)
            dinuc_freqs = calculate_dinuc_freqs(dinuc_counts)

            # seq_count (output_count), name, gene, transcript, length
            file_output.write(str(seq_count) + "\t" + seq[1] + "\t" + seq[2] + "\t" + seq[3] + "\t" + str(seq[5]))

            for base in base_freqs:
                file_output.write("\t" + str(base_freqs[base]))

            file_output.write("\t" + str(gc_content["GC"]) + "\t" + str(gc_content["AT"]))

            for dinuc in dinuc_obsexp:
                if dinuc == "NN":
                    continue

                file_output.write("\t" + str(dinuc_obsexp[dinuc]))

            for dinuc in dinuc_freqs:
                if dinuc == "NN":
                    continue

                file_output.write("\t" + str(dinuc_freqs[dinuc]))

            file_output.write("\n")


print("seqfeats_isg.py started...\n")

arguments = len(sys.argv)

if arguments == 2:
    fasta_calculate_dinucs(sys.argv[1], False)
elif arguments == 3:
    if sys.argv[2].lower() == "yes" or sys.argv[2].lower() == "y" or sys.argv[2].lower() == "true" or sys.argv[2].lower() == "1":
        arg_dups = True
    elif sys.argv[2].lower() == "no" or sys.argv[2].lower() == "n" or sys.argv[2].lower() == "false" or sys.argv[2].lower() == "0":
        arg_dups = False
    else:
        print("Error - unrecognised true/false argument can have [yes/y/true/1 or no/n/false/0]: " + sys.argv[2])
        sys.exit(1)

    fasta_calculate_dinucs(sys.argv[1], arg_dups)
else:
    print("Error - incorrect number of arguments - example usage:")
    print("seqfeats_isg.py sequences.fasta")
    print("seqfeats_isg.py sequences.fasta remove_duplicates")
    print("Exiting...")
    sys.exit(1)

print("\n...finished seqfeats_isg.py ")
