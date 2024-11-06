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

def parse_braker_annotation(annotation_file):
    gene_count = 0
    gene_lengths = []
    exon_count = 0
    exon_lengths = []
    cds_lengths = []
    intron_lengths = []

    with open(annotation_file, 'r') as ann_file:
        for line in ann_file:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if len(fields) == 9:
                    feature_type = fields[2]
                    start, end = int(fields[3]), int(fields[4])
                    feature_length = end - start + 1

                    if feature_type == 'gene':
                        gene_count += 1
                        gene_lengths.append(feature_length)

                    elif feature_type == 'exon':
                        exon_count += 1
                        exon_lengths.append(feature_length)

                    elif feature_type == 'CDS':
                        cds_lengths.append(feature_length)

                    elif feature_type == 'intron':
                        intron_lengths.append(feature_length)

    # Calculate means and ratios
    average_gene_length = sum(gene_lengths) / gene_count if gene_count else 0
    average_exons_per_gene = exon_count / gene_count if gene_count else 0
    average_exon_length = sum(exon_lengths) / exon_count if exon_count else 0
    cds_per_gene = len(cds_lengths) / gene_count if gene_count else 0
    average_cds_length = sum(cds_lengths) / len(cds_lengths) if cds_lengths else 0
    average_intron_length = sum(intron_lengths) / len(intron_lengths) if intron_lengths else 0
    total_cds_length = sum(cds_lengths)
    total_intron_length = sum(intron_lengths)

    return {
        "gene_count": gene_count,
        "average_gene_length": average_gene_length,
        "exon_count": exon_count,
        "average_exons_per_gene": average_exons_per_gene,
        "average_exon_length": average_exon_length,
        "cds_per_gene": cds_per_gene,
        "average_cds_length": average_cds_length,
        "intron_count": len(intron_lengths),
        "average_intron_length": average_intron_length,
        "total_cds_length": total_cds_length,
        "total_intron_length": total_intron_length
    }

def main():
    parser = argparse.ArgumentParser(description="Calculate genome and annotation statistics from genome and Braker files.")
    parser.add_argument("--genome", required=True, help="Path to the genome fasta file.")
    parser.add_argument("--annotation", required=True, help="Path to the Braker annotation file.")
    args = parser.parse_args()

    genome_size, gc_content = calculate_genome_statistics(args.genome)
    annotation_stats = parse_braker_annotation(args.annotation)

    gene_ratio = annotation_stats["gene_count"] / (genome_size / 1_000_000)
    exon_ratio = annotation_stats["exon_count"] / (genome_size / 1_000_000)
    cds_ratio = (annotation_stats["total_cds_length"] / genome_size) * 100
    intron_ratio = (annotation_stats["total_intron_length"] / genome_size) * 100

    print(f"Genome Size: {genome_size} bases")
    print(f"GC Content: {gc_content:.2f}%")
    print(f"Gene Count: {annotation_stats['gene_count']}")
    print(f"Gene Ratio (genes per Mb): {gene_ratio:.2f}")
    print(f"Average Gene Length: {annotation_stats['average_gene_length']:.2f} bases")
    print(f"Exon Count: {annotation_stats['exon_count']}")
    print(f"Average Exons per Gene: {annotation_stats['average_exons_per_gene']:.2f}")
    print(f"Exon Ratio (exons per Mb): {exon_ratio:.2f}")
    print(f"Average Exon Length: {annotation_stats['average_exon_length']:.2f} bases")
    print(f"Number of CDSs per Gene: {annotation_stats['cds_per_gene']:.2f}")
    print(f"Average CDS Length: {annotation_stats['average_cds_length']:.2f} bases")
    print(f"CDS Ratio (percentage of genome): {cds_ratio:.2f}%")
    print(f"Intron Count: {annotation_stats['intron_count']}")
    print(f"Average Intron Length: {annotation_stats['average_intron_length']:.2f} bases")
    print(f"Intron Ratio (percentage of genome): {intron_ratio:.2f}%")

if __name__ == "__main__":
    main()
