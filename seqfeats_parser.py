import sys


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


def check_seq_coding(this_seq, this_name, trim_seq, truncate_seq):
    seq_test = False
    new_seq = this_seq
    seq_len = len(this_seq)

    if seq_len % 3 == 0:
        seq_test = True
    else:
        print("Warning - Sequence not divisible by 3: " + this_name)
        if trim_seq:
            new_len = seq_len - (seq_len % 3)
            new_seq = this_seq[0:new_len]
            print("Trimming sequence, old length = " + str(seq_len) + ", new length = " + str(new_len))

    trunc_seq = ""

    if len(new_seq) < 3:
        print("Warning - Sequence length < 3: " + this_name)
    else:
        if new_seq[0:3] != "ATG":
            print("Warning - First codon does not equal ATG: " + new_seq[0:3] + " : " + this_name)

        for i in range(0, len(new_seq)-3, 3):
            this_codon = new_seq[i:i+3]
            trunc_seq += this_codon

            if this_codon == "TAG" or this_codon == "TAA" or this_codon == "TGA":
                print("Warning - stop codon found mid frame in seq, pos = " + str(i+1) + ": " + this_name)

                if truncate_seq:
                    new_seq = trunc_seq
                    print("Truncating sequence at first STOP, old length = " + str(seq_len) + ", new length = " + str(len(trunc_seq)))
                    break

    return new_seq


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


def fasta_process_seqs(filename, replace_names_lcl, replace_names_space,  replace_names_convert, single_genome, trim_seq, truncate_seq, is_viral):
    print("fasta_process_seqs")
    print("Input FASTA file = " + filename)
    output_filename = get_filename_stub(filename, ".f") + "_new.fasta"
    print("Outputting sequences in single line FASTA format to = " + output_filename)

    if is_viral:
        seq_stub = "ViralGenome"
    else:
        seq_stub = "Sequence"

    if replace_names_convert:
        print("Converting '>lcl|Accession_cds_existing_header' to '>Accession existing_header'")
    elif replace_names_lcl or replace_names_space:
        if replace_names_lcl:
            print("Replacing names with '>lcl|" + seq_stub + "N_cds_existing_header', where N is genome number - see Single genome setting below")
        elif replace_names_space:
            print("Replacing names with '>" + seq_stub + "N existing_header', where N is genome number - see below")

        if single_genome:
            print("Single genome = " + str(single_genome) + ", so N always = 1")
        else:
            print("Single genome = " + str(single_genome) + ", so N increments with each sequence")

    print("Trim sequence if not multiple of 3: " + str(trim_seq))
    print("Truncate sequence if STOP codon found before end: " + str(truncate_seq))

    last_line = count_lines(filename)
    line_count = 0
    seq_count = 0
    genome_count = 0

    if single_genome:
        genome_count = 1

    with open(output_filename, "w") as file_output:
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
                        if not single_genome:
                            genome_count += 1

                        if replace_names_lcl:
                            name = ">lcl|" + seq_stub + str(genome_count) + "_cds_" + name[1:]
                        elif replace_names_space:
                            name = ">" + seq_stub + str(genome_count) + " " + name[1:]
                        elif replace_names_convert:
                            name = name.replace(">lcl|", ">").replace("_cds_", " ")

                        sequence = correct_seq(sequence, name)
                        sequence = check_seq_coding(sequence, name, trim_seq, truncate_seq)
                        file_output.write(name + "\n" + sequence + "\n")

                        name = ''
                        sequence = ''

                    if line.find(">") == 0:
                        name = line
                        seq_count += 1

    print("Total sequences in file  = " + str(seq_count))

    if replace_names_lcl or replace_names_space:
        print("Total genomes  = " + str(genome_count))


print("seqfeats_parser.py started...")

arguments = len(sys.argv)

arg_replace_lcl = False
arg_replace_space = False
arg_replace_convert = False
arg_single_genome = False
arg_trim_seq = False
arg_truncate_seq = False
arg_virus = True

if arguments < 2:
    print("Error - incorrect number of arguments [" + str(arguments-1) + "] - example usage:\n")
    print("seqfeats_parser.py sequences.fasta")
    print("seqfeats_parser.py sequences.fasta [optional]replace-names-lcl [optional]replace-names-space [optional]replace-names-convert [optional]single-genome [optional]trim-seq [otpional]truncate-seq")
    print("\nreplace-names-lcl: will replace the '>' with '>lcl|ViralGenomeN_cds_ExistingHeader' where N is genome number (see below).")
    print("replace-names-space: will replace '>' with '>ViralGenomeN ExistingHeader' - where N is genome number (see below)")
    print("replace-names-convert: convert '>lcl|Accession_cds_ExistingHeader' to '>Accession ExistingHeader'")
    print("single-genome: all the coding sequences are from the same genome, so N in ViralGenomeN (see above) is always 1. If not set, N increments with each sequence")
    print("trim-seq: if coding sequence not a multiple of three then trim off last incomplete codon")
    print("truncate_seq: if coding sequence has a STOP before the last codon, then truncate the sequence at the first STOP codon")
    # not-virus - keep it hidden just now - for host transcripts - add later
    # fix_slippage option
    print("\nOutput will be in single line FASTA format - all sequence on single line")
    print("Exiting...")
    sys.exit(1)

if arguments > 2:
    for argu in sys.argv:
        if sys.argv.index(argu) < 2:
            continue

        if argu.lower() == "replace-names-lcl":
            arg_replace_lcl = True
        elif argu.lower() == "replace-names-space":
            arg_replace_space = True
        elif argu.lower() == "replace-names-convert":
            arg_replace_convert = True
        elif argu.lower() == "single-genome":
            arg_single_genome = True
        elif argu.lower() == "trim-seq":
            arg_trim_seq = True
        elif argu.lower() == "truncate-seq":
            arg_truncate_seq = True
        elif argu.lower() == "not-virus":
            arg_virus = False
        else:
            print("Error - unrecognised argument: " + argu)
            sys.exit(1)

# can only have one of these set
if arg_replace_convert:
    arg_replace_lcl = False
    arg_replace_space = False
elif arg_replace_space:
    arg_replace_lcl = False

if arg_single_genome:
    if not arg_replace_lcl and not arg_replace_space:
        print("Warning - setting single-genome without replace-names-lcl or replace-names-space will have no effect")

fasta_process_seqs(sys.argv[1], arg_replace_lcl, arg_replace_space, arg_replace_convert, arg_single_genome, arg_trim_seq, arg_truncate_seq, arg_virus)

print("\n...finished seqfeats_parser.py")
