# genome_annotation_stats
It can use for genome_annotation_report file

Yon can use command line as follows.

python genome_annotation_stats.py --genome genome.fasta --gff3 genome.gff3 --amino_acids aa.fasta --cds cds.fasta

Before use it, you have to install biopython

then output report come out as follows.

Genome Size: .. bases
GC Content: ..
Gene Count: ..
Gene Ratio (genes per Mb): ..
Average Gene Length: .. bases
Exon Count: ..
Average Exons per Gene: ..
Exon Ratio (exons per Mb): ..
Number of CDSs per Gene: .
Average CDS Length: .. bases
CDS Ratio (genes with CDS per total genes): ..
Protein Count: ..
Average Protein Length: ..
Number of Protein-Coding Genes: ..
Number of Predicted Protein Sequences: ..
Mean Gene Length: .. bases

Thank you
