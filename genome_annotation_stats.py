import argparse
from Bio import SeqIO

def calculate_genome_statistics(genome_fasta_file):
    total_bases = 0
    gc_count = 0

    for record in SeqIO.parse(genome_fasta_file, 'fasta'):
        total_bases += len(record.seq)
        gc_count += record.seq.count('G') + record.seq.count('C')

    genome_size = total_bases
    gc_content = (gc_count / total_bases) * 100

    return genome_size, gc_content

def parse_gff3_file(gff3_file):
    gene_count = 0
    gene_lengths = []
    exon_count = 0
    gene_types = {}
    cds_per_gene = []

    with open(gff3_file, 'r') as gff_file:
        for line in gff_file:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if len(fields) == 9:
                    feature_type = fields[2]
                    if feature_type == 'gene':
                        gene_count += 1
                        gene_length = int(fields[4]) - int(fields[3]) + 1
                        gene_lengths.append(gene_length)

                    elif feature_type == 'exon':
                        exon_count += 1
                    elif feature_type == 'CDS':
                        cds_length = int(fields[4]) - int(fields[3]) + 1
                        gene_id = fields[8].split(';')[0].split('=')[1]
                        gene_types[gene_id] = gene_types.get(gene_id, []) + [cds_length]
                        cds_per_gene.append(gene_id)

    if gene_count == 0:
        raise ValueError("No gene features found in the GFF3 file.")

    average_gene_length = sum(gene_lengths) / len(gene_lengths)
    average_exons_per_gene = exon_count / gene_count
    average_cds_length = sum([sum(cds_lengths) for cds_lengths in gene_types.values()]) / sum([len(cds_lengths) for cds_lengths in gene_types.values()])
    cds_per_gene_count = len(set(cds_per_gene))

    return gene_count, average_gene_length, exon_count, average_exons_per_gene, average_cds_length, cds_per_gene_count

def parse_amino_acid_sequence_file(amino_acid_sequence_file):
    protein_count = 0
    total_protein_length = 0

    for record in SeqIO.parse(amino_acid_sequence_file, 'fasta'):
        protein_count += 1
        total_protein_length += len(record.seq)

    if protein_count == 0:
        raise ValueError("No protein sequences found in the amino acid sequence file.")

    average_protein_length = total_protein_length / protein_count

    return protein_count, average_protein_length

def parse_cds_sequence_file(cds_sequence_file):
    gene_lengths = []

    protein_coding_genes = set()
    predicted_protein_sequences = []

    with open(cds_sequence_file, 'r') as cds_file:
        for record in SeqIO.parse(cds_file, 'fasta'):
            gene_id = record.id.split('.')[0]  # Assuming gene IDs are before the first dot
            protein_coding_genes.add(gene_id)
            predicted_protein_sequences.append(record)
            gene_lengths.append(len(record.seq))

    if not protein_coding_genes:
        raise ValueError("No protein-coding genes found in the CDS sequence file.")

    number_of_protein_coding_genes = len(protein_coding_genes)
    number_of_predicted_protein_sequences = len(predicted_protein_sequences)
    mean_gene_length = sum(gene_lengths) / len(gene_lengths)

    return number_of_protein_coding_genes, number_of_predicted_protein_sequences, mean_gene_length

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate genome annotation statistics.", allow_abbrev=False)
    parser.add_argument("--genome", required=True, help="Path to the genome fasta file.")
    parser.add_argument("--gff3", required=True, help="Path to the annotation gff3 file.")
    parser.add_argument("--amino_acids", required=True, help="Path to the amino acid sequence file.")
    parser.add_argument("--cds", required=True, dest="cds_file", help="Path to the coding sequence (CDS) sequence file.")
    args = parser.parse_args()

    genome_size, gc_content = calculate_genome_statistics(args.genome)
    gene_count, average_gene_length, exon_count, average_exons_per_gene, average_cds_length, cds_per_gene_count = parse_gff3_file(args.gff3)
    protein_count, average_protein_length = parse_amino_acid_sequence_file(args.amino_acids)
    number_of_protein_coding_genes, number_of_predicted_protein_sequences, mean_gene_length = parse_cds_sequence_file(args.cds_file)

    gene_ratio = gene_count / (genome_size / 1000000)
    exon_ratio = exon_count / (genome_size / 1000000)
    cds_ratio = cds_per_gene_count / gene_count

    print(f"Genome Size: {genome_size} bases")
    print(f"GC Content: {gc_content:.2f}%")
    print(f"Gene Count: {gene_count}")
    print(f"Gene Ratio (genes per Mb): {gene_ratio:.2f}")
    print(f"Average Gene Length: {average_gene_length:.2f} bases")
    print(f"Exon Count: {exon_count}")
    print(f"Average Exons per Gene: {average_exons_per_gene:.2f}")
    print(f"Exon Ratio (exons per Mb): {exon_ratio:.2f}")
    print(f"Number of CDSs per Gene: {cds_per_gene_count:.2f}")
    print(f"Average CDS Length: {average_cds_length:.2f} bases")
    print(f"CDS Ratio (genes with CDS per total genes): {cds_ratio:.2f}")
    print(f"Protein Count: {protein_count}")
    print(f"Average Protein Length: {average_protein_length:.2f} amino acids")
    print(f"Number of Protein-Coding Genes: {number_of_protein_coding_genes}")
    print(f"Number of Predicted Protein Sequences: {number_of_predicted_protein_sequences}")
    print(f"Mean Gene Length: {mean_gene_length:.2f} bases")

