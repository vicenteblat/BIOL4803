#Import necessary libraries you will need
import argparse, pysam
import matplotlib.pyplot as plt

# 1. Set up the options for the program. Three positional arguments are required in the order <Genome file> <LSJ2 bam> <N2 bam>. You will lose points for requiring flagged arguments or utilizing a different order
parser = argparse.ArgumentParser()
parser.add_argument('Genome_file', help='Name of a Fasta file containing the DNA sequence of an organism')
parser.add_argument('Bam_file_LSJ2', help='Name of a bam file containing the reads of the LSJ2 strain')
parser.add_argument('Bam_file_N2', help='Name of a bam file containing the reads of the N2 strain')
args = parser.parse_args()

# 2. Open the reference and two bam files using appropriate pysam methods (Unzip the reference first using gunzip!)
genome = pysam.FastaFile(args.Genome_file)
LSJ2 = pysam.AlignmentFile(args.Bam_file_LSJ2)
N2 = pysam.AlignmentFile(args.Bam_file_N2)

# 3. Loop through all the reads of the both bam file and count the total number of reads, the number of reads with mapping quality zero, and the number of reads with one or more mismatches from the reference. For your convience I have created variables to store this info and a print statement to report it. Hint (analyze the cigartuples attribute to identify reads with mismatches without explicitly comparing their sequence to the reference)
total_LSJ2 = 0
mq0_LSJ2 = 0
mismatch_LSJ2 = 0

for read in LSJ2:
    total_LSJ2+=1
    if read.mapping_quality == 0:
        mq0_LSJ2+=1
    if len(read.cigartuples) > 1:
        mismatch_LSJ2+=1

LSJ2.reset()

print('LSJ2: ' + str(total_LSJ2) + ' total reads, ' + str(mq0_LSJ2) + ' have a mapping quality of zero, and ' + str(mismatch_LSJ2) + ' contain one or more mismatches')
total_N2 = 0
mq0_N2 = 0
mismatch_N2 = 0

for read in N2:
    total_N2+=1
    if read.mapping_quality == 0:
        mq0_N2+=1
    if len(read.cigartuples) > 1:
        mismatch_N2+=1

N2.reset()

print('N2: ' + str(total_N2) + ' total reads, ' + str(mq0_N2) + ' have a mapping quality of zero, and ' + str(mismatch_N2) + ' contain one or more mismatches')

# 4. Create code to calculate the # of reads that mapped every 50,000 bp from the LSJ2 and N2 file (i.e x1 reads map between 0-50000, x2 reads map between 50000-100000). Save this information as lists (you don't need to keep track of the chrom and location yet) Hint: Use the count method provided by the AlignmentFile object
bin_size = 50000
coverages_LSJ2 = []
coverages_N2 = []
chrom = 'CHROMOSOME_II'
len_LSJ2_chrom = LSJ2.get_reference_length(chrom)
len_N2_chrom = N2.get_reference_length(chrom)

for i in range(0, len_LSJ2_chrom, bin_size):
    coverages_LSJ2.append(LSJ2.count(chrom, i, i + bin_size))

for i in range(0, len_N2_chrom, bin_size):
    coverages_N2.append(N2.count(chrom, i, i + bin_size))


# 5. Create a scatter plot plotting these two lists. Make sure you add appropriate labels.
plt.scatter(coverages_LSJ2,coverages_N2,alpha=0.4)
plt.plot([0,3*total_LSJ2*bin_size/15279345],[0,3*total_N2*bin_size/15279345]) # This line adds a reference line that normalizes for the coverages
plt.title('Scatterplot of reads binned from LSJ2 and N2 bam files')
plt.xlabel('# LSJ2 reads (50,000 bp bins)')
plt.ylabel('# N2 reads (50,000 bp bins)')
plt.savefig('figure1_Blat.png')
plt.show()

# 6. Once you have ran this successfully to create the above plot, add logic to print out the bins that are falling off the line. Note that they should be next to each other!

# 7. Once you have the minimum and maximum regions of the bins that are unusual, recalculate the coverage 10kb downstream and upstream of these bins with a bin size of 5000 bp. Plot this using matplotlib with appropriate labels

locations = [] # Keep track of the bins
unusual_region_LSJ2 = []
unusual_region_N2 = []
chrom = "CHROMOSOME_II"
bin_size = 5000

plt.plot(locations,unusual_region_LSJ2)
plt.plot(locations,unusual_region_N2)
plt.show()

# 8. Finally, use the pileup command to loop through all the positions of the genome and all the reads from LSJ2. Identify all the positions where there are two or more reads with mismatches to the reference. Calculate mismatch frequency in LSJ2 (num_mismatches/total_reads) and also calculate the mismatch frequency in N2. Store this data in lists.
LSJ2_mismatch_freq = []
N2_mismatch_freq = []
chrom = "CHROMOSOME_II"

# 9. Create a scatter plot of the two frequencies with appropriate labels

plt.scatterplot(LSJ2_mismatch_freq,N2_mismatch_freq)
plot.show()
