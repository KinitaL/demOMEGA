import os
import sys
import random
import re
import shutil
import seaborn as sns
import yaml
import subprocess
import glob
from Bio.Seq import Seq
from matplotlib import pyplot as plt

DEFAULT_INPUT_DIR="./input"
DEFAULT_OCCUPANCY_THRESHOLD=0.85
DEFAULT_SEQUENCE_TYPE="NT"
PATTERN_FOR_TAXA_NAMES_IN_TREE = r'([A-Za-z0-9_]+)(?=:)'

# Function to create a new directory, after checking for it's existence
def create_new_directory(dir_name):
    if os.path.isdir(dir_name):
        subprocess.run(f"rm -rf {dir_name}",
                    shell = True)
        subprocess.run(f"mkdir -p {dir_name}",
                        shell = True)
    else:
        subprocess.run(f"mkdir -p {dir_name}",
                        shell = True)

def exec_busco(input_dir, output_dir, busco_db):
    """
    Execute BUSCO to obtain complete single-copy genes
    :param input_dir: the directory with fasta files
    :param output_dir: the output directory
    :param busco_db: the DB of the BUSCO we are goint to use
    :return:
    """
    subprocess.run(
        f"busco -i {input_dir} -o {output_dir}/busco_output -m genome -l {busco_db}",
        shell=True,
    )

# TODO: simplify it
def extract_single_copy_genes(output_dir, busco_dir, occupancy_threshold, sequence_type):
    """
    Aggregate BUSCO complete single-copy genes into separate fasta files
    :param output_dir: the output directory
    :param busco_dir: the busco output directory
    :param occupancy_threshold: genes with the occupancy lower than this threshold will be excluded
    :param sequence_type: the type of sequence
    :return:
    """
    # Get the list of unique identifiers per each species
    speciesID_list = [os.path.basename(i).split("_")[0] for i in os.listdir(busco_dir)]
    print(f"{len(speciesID_list)} species identifier found in {busco_dir}/:")
    print("    " + ", ".join(speciesID_list) + ".")
    print(
        f"    {int(round(len(speciesID_list) * occupancy_threshold))} is the minimum number of species to keep a gene for downstream analysis\n"
        f"    (occupancy threshold: {int(occupancy_threshold * 100)}%).")

    # Get the list of all the complete BUSCO genes found in the input species
    gene_file = output_dir + "/complete_busco_genes.ls"

    print("Retrieving the list of all the complete BUSCO genes found in the input species...")
    try:
        awk_process = subprocess.run(f"cat {busco_dir}/*/run*/full_table.tsv | "
                                     "awk -F \"\t\" '{if ($2 == \"Complete\") print $1}' | "
                                     "sort -u | "
                                     f"shuf > {gene_file}",
                                     shell=True)

        awk_process.check_returncode()

    except subprocess.CalledProcessError as err:
        print("An error occured:", err)

    # Store the obtained gene IDs in a list
    gene_list = []
    with open(gene_file) as input_list:
        for line in input_list.readlines():
            gene_list.append(line.strip())

    print(f"    {len(gene_list)} unique BUSCO gene identifiers found.")

    # Create a directory to store fasta files
    fasta_dir = f"{output_dir}/fasta_files"

    # Remove the dir in order to avoid doubles
    create_new_directory(fasta_dir)

    # Counters for discarded and kept genes
    discarded_genes = 0
    kept_genes = 0

    if sequence_type == "NT":
        if not glob.glob(f"{busco_dir}/*/run*/busco_sequences/single_copy_busco_sequences/*.fna"):
            print()
            print(
                f"### Mmh, you selected the {sequence_type} fasta type (*.fna), but no such file has been found. ###\n"
                f"### You may want to check your input directory ({busco_dir}/). ###\n"
                "Quitting the analysis...")
            print()

            quit()
        else:
            print(f"    The analyses will be performed on nucleotide sequences!")
            print()

    elif sequence_type == "AA":
        if not f"{busco_dir}/*/run*/busco_sequences/single_copy_busco_sequences/*.faa":
            print()
            print(
                f"### Mmh, you selected the {sequence_type} fasta type (*.faa), but no such file has been found. ###\n"
                f"### You may want to check your input directory ({busco_dir}/). ###\n"
                "Quitting the analysis...")
            print()

            quit()

        else:
            print(f"    The analyses will be performed on amino acid sequences!")
            print()

    # Generate separate fasta files with genes from species meeting inclusion threshold
    print(f"Generating fasta files for each gene with all the annotated species...")
    for gene in gene_list:

        if sequence_type == "NT":
            # Get the list of each genes in each species (nucleotides)
            gene_file_list = glob.glob(f"{busco_dir}/*/run*/busco_sequences/single_copy_busco_sequences/{gene}.fna")

        elif sequence_type == "AA":
            # Get the list of each genes in each species (amino acids)
            gene_file_list = glob.glob(f"{busco_dir}/*/run*/busco_sequences/single_copy_busco_sequences/{gene}.faa")

        # Get the occupancy thershold for each gene
        occupancy = len(gene_file_list) / len(os.listdir(busco_dir))

        # print(f"{gene=} {occupancy=} {args.occupancy_threshold=}")

        # Check if occupancy requirements are met, in which case do not continue with the analysis
        if float(occupancy) < float(occupancy_threshold):
            discarded_genes += 1
            continue
        else:
            kept_genes += 1

        # For each gene in each species, change the header to include just the unique species identifier
        for id in speciesID_list:
            if sequence_type == "NT":
                input_fasta = glob.glob(
                    f"{busco_dir}/{id}*/run*/busco_sequences/single_copy_busco_sequences/{gene}.fna")
                output_fasta = f"{output_dir}/fasta_files/{gene}.fna"
            elif sequence_type == "AA":
                input_fasta = glob.glob(
                    f"{busco_dir}/{id}*/run*/busco_sequences/single_copy_busco_sequences/{gene}.faa")
                output_fasta = f"{output_dir}/fasta_files/{gene}.faa"

            if input_fasta:
                # Progressively add sequences to each fasta file
                subprocess.run(f"touch {output_fasta}",
                               shell=True,
                               text=True)

                # Do not write doubles
                species_in_file_count = subprocess.run(f"tr -cd \">\" < {output_fasta} | wc -c",
                               shell=True,
                               capture_output=True)

                if int(species_in_file_count.stdout) >= len(speciesID_list):
                    continue

                subprocess.run(f"cat {input_fasta[0]} | sed -E 's/^>.+$/>{id}/' >> {output_fasta}",
                               shell=True,
                               text=True)
            else:
                continue

    if discarded_genes == 1:
        print(
            f"    {discarded_genes} gene was discarded because didn't meet the occupancy threshold ({occupancy_threshold}).")
    else:
        print(
            f"    {discarded_genes} genes were discarded because didn't meet the occupancy threshold ({occupancy_threshold}).")

    print(f"    {kept_genes} genes were kept for downstream phylogenomic analysis.")
    print("Genes extraction is done")


def translate_nucleotides_to_protein(input_dir, output_dir):
    """
    Translates nucleotide sequences of the gene to amino acid sequences
    :param input_dir: the input directory with genes represented as nucleotide sequences
    :param output_dir: the output directory with amino acid sequences files
    :return:
    """
    os.makedirs(output_dir, exist_ok=True)
    for subdir, dirs, files in os.walk(input_dir):
        for file in files:
            with open(f"{input_dir}/{file}", 'r') as input, open(f"{output_dir}/{file.split(".")[0]}.fa", 'w') as output:
                to_translate = ""
                for line in input:
                    if line.startswith('>'):
                        if to_translate != "":
                            output.write(str(Seq(to_translate).translate()))
                        output.write(f"\n{line.strip()}\n")
                        to_translate = ""
                    else:
                        to_translate += line.strip()
                if to_translate != "":
                    output.write(str(Seq(to_translate).translate()))


def align_with_mafft(input_dir, output_dir):
    """
    Align protein sequences using MAFFT
    :param input_dir: the directory with genes represented as nucleotide sequences
    :param output_dir: the output directory with alignments
    :return:
    """
    os.makedirs(output_dir, exist_ok=True)
    for subdir, dirs, files in os.walk(input_dir):
        for file in files:
            subprocess.run(
                f"mafft --auto {input_dir}/{file} > {output_dir}/{file}",
                shell=True,
            )


def translate_to_codon_alignments(input_dir, output_dir, sequence_type, nucleotides_dir):
    """
    Translates protein alignments to codon alignments
    :param input_dir: the input directory with protein alignments
    :param output_dir: the output directory with codon alignments
    :param sequence_type: the type of sequence
    :param nucleotides_dir: the directory with nucleotide sequences
    :return:
    """
    if sequence_type == "NT":
        file_format = ".fna"
        nt_dir = f"{input_dir}/fasta_files"
    else:
        file_format = ".faa"
        nt_dir = nucleotides_dir
    os.makedirs(output_dir, exist_ok=True)
    for subdir, dirs, files in os.walk(f"{input_dir}/protein_alignments"):
        for file in files:
            # # convert the sequence to upper case
            # with open(f"{input_dir}/fasta_files/{file.split(".")[0]}{file_format}", 'r') as fasta_file:
            #     lines = fasta_file.readlines()
            #
            # with open(f"{input_dir}/fasta_files/{file.split(".")[0]}{file_format}", 'w') as out:
            #     for fasta_line in lines:
            #         if fasta_line.startswith('>'):
            #             out.write(fasta_line.rstrip('\n') + '\n')
            #         else:
            #             out.write(fasta_line.upper())

            # execute pal2nal
            subprocess.run(
                f"pal2nal.pl {input_dir}/protein_alignments/{file} {nt_dir}/{file.split(".")[0]}{file_format} -output fasta > {output_dir}/{file.split(".")[0]}{file_format}",
                shell=True,
            )

def concat(input_dir, output_dir, sequence_type):
    """
    Concatenates genes
    :param input_dir: the input directory
    :param output_dir: the output directory
    :return:
    """
    if sequence_type == "NT":
        file_format = ".fna"
    else:
        file_format = ".faa"
    os.makedirs(output_dir, exist_ok=True)
    for subdir, dirs, files in os.walk(input_dir):
        # input = f"{input_dir}/"+f" {input_dir}/".join(files)
        subprocess.run(
            f"AMAS.py concat -i {input_dir}/*{file_format} -f fasta -d dna -t {output_dir}/concatenated.out -u fasta -p {output_dir}/partitions.txt --part-format nexus",
            shell=True,
        )

def calculate_dnds_for_species_tree(input_dir, output_dir, tree_dir):
    """
    Calculate dN/dS using codeml
    :param input_dir: the input directory
    :param output_dir: the output directory
    :param tree_dir: the tree directory
    :return:
    """
    os.makedirs(output_dir, exist_ok=True)
    build_config(input_dir, output_dir, tree_dir)

    subprocess.run(
        f"codeml {output_dir}/codeml.ctl",
        shell=True,
        input="\n",
        text=True,
    )

def construct_tree(input_dir, output_dir):
    """
    Construct a tree
    :param input_dir: the input directory
    :param output_dir: the output directory
    :return:
    """
    os.makedirs(output_dir, exist_ok=True)
    subprocess.run(
        f"iqtree -s {input_dir}/concatenated.out -pre {output_dir}/tree",
        shell=True,
    )

def label_tree(tree_dir, attach_labels = False):
    """
    Assign labels to the tree
    :param tree_dir: a directory where the tree is stored
    :return:
    """
    if attach_labels:
        with open(f"{tree_dir}/tree.treefile", "r") as file, open(f"{tree_dir}/labeled_tree.treefile", "w") as new_file:
            tree = file.read()
            modified_tree = re.sub(PATTERN_FOR_TAXA_NAMES_IN_TREE, add_suffix, tree)
            new_file.write(modified_tree)
    else:
        shutil.copyfile(f"{tree_dir}/tree.treefile", f"{tree_dir}/labeled_tree.treefile")


def add_suffix(match):
    return f"{match.group(1)}#1"

def build_config(input_dir, output_dir, tree_dir):
    """
    Build codeml config using the template
    :param input_dir: a directory when the concatenation file is stored
    :param output_dir: a directory when the output is going to be saved
    :param tree_dir: a directory when the tree is going to be saved
    :return:
    """
    with open("./codeml.ctl.template", "r") as template, open(f"{output_dir}/codeml.ctl", "w") as conf:
        for line in template:
            if "$seqfile" in line:
                line = line.replace("$seqfile", f"{input_dir}/concat/concatenated.out", -1)
            if "$outfile" in line:
                line = line.replace("$outfile", f"{output_dir}/codeml.out")
            if "$treefile" in line:
                line = line.replace("$treefile", f"{tree_dir}/labeled_tree.treefile")
            conf.write(line)

def sample_alignment_files(dir, n):
    """
    Randomly sample n file paths from all_files without replacement.
    """
    for subdir, dirs, files in os.walk(dir):
        if n > len(files):
            raise ValueError(f"Requested {n} files, but only {len(files)} available.")
        return random.sample(files, n)

def extract_dnds_distribution_data(input_dir):
    """
    Extract dN/dS values from codeml output
    :param input_dir: codeml output dir
    :return: List of dN/dS values
    """
    # Take the tree from codeml output
    tree = subprocess.run(
        f"tac {input_dir}/codeml.out | grep -m1 -v '^$' | tr -d '[:space:]'",
        shell=True,
        capture_output=True,
    )
    # Decode bytes to string if needed
    if isinstance(tree.stdout, bytes):
        newick_str = tree.stdout.decode('utf-8')
    else:
        newick_str = tree.stdout

    branch_lengths = re.findall(r'#(\d+\.\d+)', newick_str)

    return [float(length) for length in branch_lengths]

def draw_distributions(data, dir):
    """
    Draw and save to a file a distribution of dN/dS for a branch
    :param data: the data for the distribution
    :param dir: the path to the dir where the distribution will be saved
    :return:
    """
    plt.figure()
    fig = sns.histplot(data, kde=True, stat="probability").get_figure()
    fig.savefig(f"{dir}/distribution.png")

def label_final_tree(tree_dir, results_dir, means):
    """
    Draw labels for the final tree
    :param tree_dir: the path to the file with the final tree
    :param results_dir: the path to the directory with the results
    :param means: the array with means of dN/dS for each branch
    :return:
    """
    with open(f"{tree_dir}/tree.treefile", "r") as file, open(f"{results_dir}/tree.treefile", "w") as new_file:
        tree = file.read()
        parts = tree.split(':')
        n_colons = len(parts) - 1
        if len(means) != n_colons:
            raise ValueError(f"{n_colons} colons but {len(means)} values provided")

        out = parts[0]
        for val, suffix in zip(means, parts[1:]):
            out += f"#{val}:{suffix}"

        new_file.write(out)


if __name__ == '__main__':
    # Read config
    with open('config.yaml', 'r') as file:
        config = yaml.safe_load(file)

    # Create the output directory
    if "output_dir" in config and config["output_dir"] is not None:
        print(f"Using {config['output_dir']} as a root output directory.")
        os.makedirs(config["output_dir"], exist_ok=True)
    else:
        sys.exit("No output directory specified (see config)")

    # extract sequence type from config for simpler usage
    sequence_type = config["sequence_type"] if "sequence_type" in config else DEFAULT_SEQUENCE_TYPE

    # check the presence of nucleotides sequences for pal2nal
    if sequence_type == "AA" and not ("nt_dir" in config and config["nt_dir"] is not None):
        sys.exit("It is required to specify the directory with nucleotide sequences.")
    else:
        # Execute busco for NT sequences since we have to use it for pal2nal.
        exec_busco(
            config["nt_dir"],
            f"{config["output_dir"]}/nt",
            config["busco_db"],
        )
        extract_single_copy_genes(
            f"{config["output_dir"]}/nt",
            busco_dir=f"{config["output_dir"]}/nt/busco_output",
            occupancy_threshold=config[
                "occupancy_threshold"] if "occupancy_threshold" in config else DEFAULT_OCCUPANCY_THRESHOLD,
            sequence_type="NT",
        )

    # Prepare busco output
    if "input_type" in config and config["input_type"] is not None:
        print(f"Using {config['input_type']} as an input type.")
    else:
        sys.exit("No input type is specified (see config)")

    if config["input_type"] == "B":
        # extract genes from busco output
        extract_single_copy_genes(
            config["output_dir"],
            busco_dir=config["input_dir"],
            occupancy_threshold=config[
                "occupancy_threshold"] if "occupancy_threshold" in config else DEFAULT_OCCUPANCY_THRESHOLD,
            sequence_type=sequence_type,
        )

    elif config["input_type"] == "F":
        # execute busco
        exec_busco(
            config["input_dir"],
            config["output_dir"],
            config["busco_db"],
        )

        # extract genes from busco output
        extract_single_copy_genes(
            config["output_dir"],
            busco_dir=f"{config["output_dir"]}/busco_output",
            occupancy_threshold=config[
                "occupancy_threshold"] if "occupancy_threshold" in config else DEFAULT_OCCUPANCY_THRESHOLD,
            sequence_type=sequence_type,
        )
    else:
        sys.exit(f"Input type {config["input_type"]} not recognized (see config)")

    # TODO: remove paralogs

    if sequence_type == "NT":
        translate_nucleotides_to_protein(
            f"{config["output_dir"]}/fasta_files",
            f"{config["output_dir"]}/aa_seq_files",
        )
    else:
        os.makedirs(f"{config["output_dir"]}/aa_seq_files", exist_ok=True)
        subprocess.run(
            f"cp -R {config["output_dir"]}/fasta_files/. {config["output_dir"]}/aa_seq_files",
            shell=True,
            capture_output=True, )

    align_with_mafft(
        f"{config["output_dir"]}/aa_seq_files",
        f"{config["output_dir"]}/protein_alignments",
    )

    translate_to_codon_alignments(
        f"{config["output_dir"]}",
        f"{config["output_dir"]}/codon_alignments",
        sequence_type,
        f"{config["output_dir"]}/nt/fasta_files",
    )

    concat(
        f"{config["output_dir"]}/codon_alignments",
        f"{config["output_dir"]}/concat",
        sequence_type
    )

    # Use a provided tree or construct a tree using concatenation
    if "tree_template" in config and config["tree_template"] is not None and os.path.isfile(config["tree_template"]):
        shutil.move(config["tree_template"], f"{config["output_dir"]}/tree/tree.treefile")
    else:
        construct_tree(
            f"{config["output_dir"]}/concat",
            f"{config["output_dir"]}/tree",
        )

    # Label species tree
    tree_dir = f"{config["output_dir"]}/tree"

    need_labeling = False
    with open("./codeml.ctl.template", "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("model"):
                parts = line.split('=')
                if len(parts) == 2:
                    if int(parts[1].strip()) == 2:
                        need_labeling = True
    label_tree(tree_dir, need_labeling)


    # Here goes the following:
    # - If no bootstrapping:
    #   - calculate dn/ds using species tree
    # - If bootstrapping:
    #   - in each iteration:
    #       - sample from codon_alignments
    #       - concat the samples taken in this iteration
    #       - calculate dn/ds using concatenation and species tree (calculated priorly or provided by the user)
    #   - build a distribution of dn/ds for each branch
    if "bootstrapping" in config and "do" in config["bootstrapping"] and config["bootstrapping"]["do"]:
        os.makedirs(f"{config["output_dir"]}/bootstrapping", exist_ok=True)
        sample_number = config["bootstrapping"]["sample_number"] if "sample_number" in config["bootstrapping"] else 0
        sample_size = config["bootstrapping"]["sample_size"] if "sample_size" in config["bootstrapping"] else 0
        if sample_number <= 0:
            raise ValueError("You use bootstrapping option but do not provide sample_number: sample number must be greater than zero")
        if sample_size <= 0:
            raise ValueError("You use bootstrapping option but do not provide sample_size: sample size must be greater than zero")
        counter = 0
        distribution_data = []
        while counter < sample_number:
            create_new_directory(f"{config["output_dir"]}/bootstrapping/{counter + 1}")
            files = sample_alignment_files(f"{config["output_dir"]}/codon_alignments", sample_size)
            for file in files:
                create_new_directory(f"{config["output_dir"]}/bootstrapping/{counter + 1}/codon_alignments")
                shutil.copyfile(f"{config["output_dir"]}/codon_alignments/{file}", f"{config["output_dir"]}/bootstrapping/{counter + 1}/codon_alignments/{file}")
            concat(
                f"{config["output_dir"]}/bootstrapping/{counter + 1}/codon_alignments",
                f"{config['output_dir']}/bootstrapping/{counter + 1}/concat",
            )
            calculate_dnds_for_species_tree(
                f"{config["output_dir"]}/bootstrapping/{counter + 1}",
                f"{config["output_dir"]}/bootstrapping/{counter + 1}/dnds",
                tree_dir,
            )
            ratios = extract_dnds_distribution_data(f"{config["output_dir"]}/bootstrapping/{counter + 1}/dnds")
            for i in range(len(ratios)):
                if counter == 0:
                    distribution_data.append([])
                distribution_data[i].append(ratios[i])
            counter += 1
        create_new_directory(f"{config["output_dir"]}/results/distributions")
        means = []
        for i in range(len(distribution_data)):
            create_new_directory(f"{config["output_dir"]}/results/distributions/{i + 1}")
            draw_distributions(distribution_data[i],f"{config["output_dir"]}/results/distributions/{i + 1}")
            means.append(sum(distribution_data[i]) / len(distribution_data[i]))
        label_final_tree(tree_dir, f"{config["output_dir"]}/results", means)
    else:
        calculate_dnds_for_species_tree(
            f"{config["output_dir"]}",
            f"{config["output_dir"]}/dnds",
            tree_dir,
        )
        means = extract_dnds_distribution_data(f"{config["output_dir"]}/dnds")
        label_final_tree(tree_dir, f"{config["output_dir"]}/results", means)

    # Draw a tree with dN/dS attached
    subprocess.run(
        f"Rscript plot.R --file={config["output_dir"]}/results",
        shell=True,
        capture_output=True,)