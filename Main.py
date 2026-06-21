#Load modules
import regex as re
import sys
import os
import time
import csv
import multiprocessing as mp
import random
from math import *
import os.path
import subprocess
import scipy.stats as stats
start_time = time.time()
seed_number=135
random.seed(seed_number)

# Launched per-element by run_exinator_el.sh. 

# $PYTHON "$MAIN_SCRIPT" -i "$MUTATION_FILE" -o "$OUTPUT_DIR" -f "$FASTA_FILE" -g "$GTF_FILE" -e "$BED_FILE" -r "$BLACKLIST_FILE" -s "$CHROM_FILE" -k "$KMERS_FILE" -w "$WHOLE_GTF" -c $CORES -n $N -b $BOOTSTRAPS -ss $STRAND
#
# Equivalent direct call for:
# BLCA cohort, strand-specific, "rbp_combined" element - Path to be changed accordingly
#
# python2.7 Main.py \
#   -i  /path/to/base_dir/bed_files/unique_mutations_cancerwise/BLCA_passvariants_exinator2_sorted_uniq.bed \
#   -o  /path/to/pipeline_outputs/Outputs/strand_specific/BLCA/rbp_combined \
#   -f  /path/to/base_dir/Inputs/hg38_nochr_allcaps.fa \
#   -g  /path/to/cohort.gtf \
#   -e  /path/to/base_dir/element_bed_files/rbp_combined.bed \
#   -r  /path/to/base_dir/blacklist/ExInAtor_final_blacklist_hg38_sorted_merged_boyle_GIAB_segdup_merged.bed \
#   -s  /path/to/base_dir/Inputs/chromosomes_38.bed \
#   -k  /path/to/base_dir/Inputs/3mers.txt \
#   -w  /path/to/base_dir/Inputs/gencode.v39.long_noncoding_RNAs_noPARY_nochrM_cpatfilter_noPCG_CLC3.gtf \
#   -c  6  -n 119  -b 10000  -ss true

if len(sys.argv)==1 or "-h" in sys.argv: 
    print "Mandatory arguments: -i <--input_file> -o <--output_folder> -f <--fasta_file> -g <--gtf_file> -e <--element_file> -r <--exc_regions> -s <--chr_sizes> -k <--kmers_file> -w <--whole_genome> -n <--number_of_genomes>"
    sys.exit()
n=1
while n<len(sys.argv):
    arg = sys.argv[n]
    if arg == "-i" or arg == "--input_file":
        input_file = sys.argv[n+1]
    elif arg == "-o" or arg == "--output_folder":
        output_folder = sys.argv[n+1]
    elif arg == "-f" or arg == "--fasta_file":
        fasta_file = sys.argv[n+1]
    elif arg == "-g" or arg == "--gtf_file":
        gtf_file = sys.argv[n+1]
    elif arg == "-e" or arg == "--element_file":
        element_file = sys.argv[n+1]
    elif arg == "-r" or arg == "--exc_regions": #BED file with regions from the genome to ignore(such as those with low mappability, high repetitive sequences,etc).
        exc_regions = sys.argv[n+1]
    elif arg == "-s" or arg == "--chr_sizes":
        chr_sizes = sys.argv[n+1]
    elif arg == "-k" or arg == "--kmers_file":
        kmers_file = sys.argv[n+1]
    elif arg == "-w" or arg == "--whole_genome":
        whole_genome = sys.argv[n+1]
    elif arg == "-b" or arg == "--background_size":
        background_size = sys.argv[n+1]
    elif arg == "-n" or arg == "--number_of_genomes":
        number_of_genomes = sys.argv[n+1]
    elif arg == "-c" or arg == "--cores":
        cores = int(sys.argv[n+1])
    elif arg == "-y" or arg == "--exonic_filter":
        exonic_filter = sys.argv[n+1]
    elif arg == "-x" or arg == "--background_filter":
        background_filter = sys.argv[n+1]    
    elif arg == "-ss" or arg == "--strand_specific":
        strand_specific = sys.argv[n+1].lower() in ['true', 'yes', '1']
    n+=2
if 'input_file' not in globals():
    print "Please indicate \"-i\" or \"--input_file\""
    print "type \"Main.py -h\" for help"
    sys.exit()
if 'output_folder' not in globals():
    print "Please indicate \"-o\" or \"--output_folder\""
    print "type \"Main.py -h\" for help"
    sys.exit()
if 'fasta_file' not in globals():
    print "Please indicate \"-f\" or \"--fasta_file\""
    print "type \"Main.py -h\" for help"
    sys.exit()
if 'gtf_file' not in globals():
    print "Please indicate \"-g\" or \"--gtf_file\""
    print "type \"v1.py -h\" for help"
    sys.exit()
if 'element_file' not in globals():
    print "Please indicate \"-e\" or \"--element_file\""
    print "type \"v1.py -h\" for help"
    sys.exit()
if 'exc_regions' not in globals():
    exc_regions="NO"
if 'chr_sizes' not in globals():
    print "Please indicate \"-s\" or \"--chr_sizes\""
    print "type \"Main.py -h\" for help"
    sys.exit()
if 'kmers_file' not in globals():
    print "Please indicate \"-k\" or \"--kmers_file\""
    print "type \"Main.py -h\" for help"
    sys.exit()
if 'whole_genome' not in globals():
    print "Please indicate \"-w\" or \"--whole_genome\""
    print "type \"Main.py -h\" for help"
    sys.exit()
if 'number_of_genomes' not in globals():
    print "Please indicate \"-n\" or \"--number_of_genomes\""
    print "type \"Main.py -h\" for help"
    sys.exit()
if 'background_size' not in globals():
    background_size=str(10000)
if 'cores' not in globals():
    cores=1
if 'exonic_filter' not in globals():
    exonic_filter=str(1)
if 'background_filter' not in globals():
    background_filter=str(1)
if 'strand_specific' not in globals():
    strand_specific=False  # Default to False if not specified
    

print("Arguments correctly inserted: %.0f seconds " % (time.time() - start_time))
del arg

print("Starting the analysis of "+input_file.split(".")[0]+": %.0f seconds " % (time.time() - start_time))

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
input_folder=os.getcwd()+"/"
os.chdir(output_folder)
try:
    # Get genes

    genes=[]
    file = open ("genes.bed","r")
    for line in file:
        line=line.rstrip()
        line=line.split("\t")
        genes.append(line[3])

    print("\tGenes gathered: %.0f seconds " % (time.time() - start_time))

    # Get kmers

    file = open (kmers_file,"r")

    kmers=[]
    for line in file:
        line=line.rstrip()
        kmers.append(line)
    file.close()

    print("\tKmers gathered: %.0f seconds " % (time.time() - start_time))

except: adsada=12

if not os.path.isfile("table_kmer_counts.txt"):
        
    # Create a BED file of the merged exons, introns and background for each lncRNA gene.
    subprocess.call("cat "+whole_genome+' | tr --delete \; | tr --delete \\" | grep -v \\# | awk \'$3==\"exon\"\' | awk \'{split($10,array,\".\")}{print $1\"\t\"$4-1\"\t\"$5\"\t\"$10\"\t.\t\"$7}\' | sortBed -i - > exons_n2.bed', shell=True)
    subprocess.call("mergeBed -i exons_n2.bed -s -c 4,6 -o distinct,distinct > temp_merged_exons_n2.bed", shell=True)
    subprocess.call("awk '{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}' temp_merged_exons_n2.bed > merged_exons_fn2.bed", shell=True)

    # merge exons of lncrna
    subprocess.call("cat "+gtf_file+' | tr --delete \; | tr --delete \\" | grep -v \\# | awk \'$3==\"exon\"\' | awk \'{split($10,array,\".\")}{print $1\"\t\"$4-1\"\t\"$5\"\t\"$10\"\t.\t\"$7}\' |sortBed -i - > exons_n.bed', shell=True)
    subprocess.call("mergeBed -i exons_n.bed -s -c 4,6 -o distinct,distinct | awk \'{split($1,array,\",\")} {print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}\'> merged_exons_n.bed", shell=True)
    subprocess.call("cat "+gtf_file+' | tr --delete \; | tr --delete \\" | grep -v \\# | awk \'$3==\"gene\"\' | awk \'{print $1\"\t\"$4-1\"\t\"$5\"\t\"$10\"\t.\t\"$7}\'  > genes_n.bed', shell=True)
    subprocess.call("bedtools slop -i genes_n.bed -g "+chr_sizes+" -b "+background_size+" > genes_extended_n.bed", shell=True)
    

    subprocess.call("awk \'{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}\' merged_exons_n.bed | sortBed -i - > merged_exonsn.bed",shell=True)#removed $4
    subprocess.call("awk \'{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}\' merged_exons_fn2.bed | sortBed -i - > merged_exonsn2.bed",shell=True)#addednewly

    if exc_regions!="NO":
        subprocess.call("subtractBed -a merged_exonsn2.bed -b "+exc_regions+" > merged_exons_n2_clean.bed", shell=True)
        subprocess.call("subtractBed -a merged_exonsn.bed -b "+exc_regions+" > merged_exons_n_clean.bed", shell=True)
        subprocess.call("mv merged_exons_n2_clean.bed merged_exons_n2.bed", shell=True)
        subprocess.call("mv merged_exons_n_clean.bed merged_exons_n.bed", shell=True)
    subprocess.call("awk \'{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5}\' merged_exons_n.bed | sortBed -i - > merged_exons.bed",shell=True)

    # Collapse all merged exon regions under a single pooled gene label so the
    # element type yields one contingency table. The specific ID is arbitrary,
    # but it MUST exist in genes.bed (same GTF) or table_kmer_counts ends up all
    # zeros. 
    with open("merged_exons.bed") as _fh:
        pooled_gene = _fh.readline().split("\t")[3].split(",")[0].strip()
    subprocess.call(
        "awk -v OFS='\\t' -v g='%s' '{$4 = g} 1' merged_exons.bed > output.bed" % pooled_gene,
        shell=True,
    )
    if strand_specific:
        # For strand-specific analysis
        print("Performing strand-specific analysis")

        subprocess.call("mergeBed -i "+element_file+" -c 6 -o distinct -s > u_merged_element.bed", shell=True)
        subprocess.call("sortBed -i u_merged_element.bed > merged_element.bed",shell=True)

        #Merge the regions based on strands and split the regions into plus and minus stranded merged element
        # Command to split merged_element.bed into + strand
        subprocess.call("awk '$4 == \"+\"' merged_element.bed > merged_element_plus.bed", shell=True)
        # Command to split merged_element.bed into - strand
        subprocess.call("awk '$4 == \"-\"' merged_element.bed > merged_element_minus.bed", shell=True)
        print("Splitting complete: + strand -> merged_element_plus.bed, - strand -> merged_element_minus.bed")
    

        # Command to extract + stranded genes
        subprocess.call("awk '$5 == \"\\+\"' output.bed > output_plus.bed", shell=True)
        # Command to extract - stranded genes
        subprocess.call("awk '$5 == \"-\"' output.bed > output_minus.bed", shell=True)
        print("Successfully split output.bed into 'output_plus.bed' and 'output_minus.bed'")
    
        # Command to intersect element and exons for + strand
        subprocess.call("intersectBed -a output_plus.bed -b merged_element_plus.bed -wb > intersect_plus.bed", shell=True)
        # Command to intersect element and exons for - strand
        subprocess.call("intersectBed -a output_minus.bed -b merged_element_minus.bed -wb > intersect_minus.bed", shell=True)
        print("Intersection complete: 'intersect_plus.bed' and 'intersect_minus.bed' created.")


        # Combine intersect_plus.bed and intersect_minus.bed
        subprocess.call("cat intersect_plus.bed intersect_minus.bed > combined_temp.bed", shell=True)
        # Sort the combined file
        subprocess.call("sort -k1,1 -k2,2n combined_temp.bed > merged_exons_element.bed", shell=True)
        # Remove the temporary file
        subprocess.call("rm combined_temp.bed", shell=True)
        print("Combined and sorted 'intersect_plus.bed' and 'intersect_minus.bed' into 'merged_exons_element.bed'")


        subprocess.call("subtractBed -a output_plus.bed -b merged_element_plus.bed > merged_exons_nonelement_plus.bed", shell=True)
        subprocess.call("subtractBed -a output_minus.bed -b merged_element_minus.bed > merged_exons_nonelement_minus.bed", shell=True)
        subprocess.call("cat merged_exons_nonelement_plus.bed merged_exons_nonelement_minus.bed > combined_temp_nonelement.bed", shell=True)
        # Sort the combined file
        subprocess.call("sort -k1,1 -k2,2n combined_temp_nonelement.bed > merged_exons_nonelement.bed", shell=True)
        # Remove the temporary file
        subprocess.call("rm combined_temp_nonelement.bed", shell=True)
        print("Combined and sorted 'intersect_plus.bed' and 'intersect_minus.bed' into 'merged_exons_nonelement.bed'")

    else:
        # For non-strand-specific analysis
        print("Performing non-strand-specific analysis") 
        subprocess.call("mergeBed -i "+element_file+"  > merged_element.bed", shell=True)

        subprocess.call("subtractBed -a output.bed -b merged_element.bed > merged_exons_nonelement.bed", shell=True)
        subprocess.call('intersectBed -a output.bed -b merged_element.bed -sorted -wb | awk \'{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}\' > merged_exons_element.bed', shell=True)

    subprocess.call('awk \'{split($1,array,\",\")} {print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$6}\' exons_n.bed | sortBed -i - > exons.bed', shell=True)
    subprocess.call('awk \'{print $1\"\t\"$2\"\t\"$3\"\t\"$4}\' genes_extended_n.bed | sortBed -i - > genes_extended.bed', shell=True)
    subprocess.call('awk \'{print $1\"\t\"$2\"\t\"$3\"\t\"$4}\' genes_n.bed | sortBed -i - > genes.bed', shell=True)
    subprocess.call('awk \'{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5}\' merged_exons_n.bed | sortBed -i - > merged_exons.bed', shell=True)#check again
    subprocess.call('awk \'{print $1\"\t\"$2\"\t\"$3\"\t\"$4}\' merged_exons_nonelement.bed | sortBed -i - > exons_nonelement.bed', shell=True)
    subprocess.call('awk \'{print $1\"\t\"$2\"\t\"$3\"\t\"$4}\' merged_exons_element.bed | sortBed -i - > exons_element.bed', shell=True)
    print("\tBED files created for merged exonic elements and non-elements: %.0f seconds " % (time.time() - start_time))

    # Count mutation in merged exons and background

    subprocess.call('intersectBed -a exons_element.bed -b '+input_file+' -sorted -wb | awk \'{print $1"\t"$2-1"\t"$3+1"\t"$4"\t"$5"\t"$6"\t"$8}\' | sort | uniq >exons_element_mutations.bed', shell=True)
    subprocess.call('intersectBed -a exons_nonelement.bed -b '+input_file+' -sorted -wb | awk \'{print $1"\t"$2-1"\t"$3+1"\t"$4"\t"$5"\t"$6"\t"$8}\' | sort | uniq > exons_nonelement_mutations.bed', shell=True)
    subprocess.call("sed 's/chr//g' exons_element_mutations.bed > exons_element_mutations_for_fasta.bed", shell=True)
    subprocess.call("sed 's/chr//g' exons_nonelement_mutations.bed > exons_nonelement_mutations_for_fasta.bed", shell=True)
    subprocess.call('fastaFromBed -fi '+fasta_file+' -bed exons_element_mutations_for_fasta.bed -fo exons_element_mutations.fa -nameOnly', shell=True)
    subprocess.call('fastaFromBed -fi '+fasta_file+' -bed exons_nonelement_mutations_for_fasta.bed -fo exons_nonelement_mutations.fa -nameOnly', shell=True)

    print("\tMutations counted in merged exons elements and non-elements: %.0f seconds " % (time.time() - start_time))

    # Get genes

    genes=[]
    file = open ("genes.bed","r")
    for line in file:
        line=line.rstrip()
        line=line.split("\t")
        genes.append(line[3])

    print("\tGenes gathered: %.0f seconds " % (time.time() - start_time))

    # Get kmers

    file = open (kmers_file,"r")

    kmers=[]
    for line in file:
        line=line.rstrip()
        kmers.append(line)
    file.close()

    print("\tKmers gathered: %.0f seconds " % (time.time() - start_time))
    
    # Function to count kmers
    
    def count_kmers(list):
        file = open (list[0],"r")
        kmers=list[1]
        order=list[2]
        exons={}
        header = None
        for line in file:
            line=line.rstrip()
            if line.startswith('>'):
                header = line[1:]
            else:
                sequence=line
                if header not in exons:
                    exons[header]={}
                for motif in kmers:
                    if motif not in exons[header]:
                        exons[header][motif]=0
                    if motif in exons[header]:
                        exons[header][motif]+=len(re.findall(motif, sequence, overlapped=True))
        return((order,exons))
    
    # Function to print the counted kmers
    
    def print_kmers(list):
        file=list[1]
        list=list[0]
        file=open(file,"w")
        for header in list:
            for motif in list[header]:
                line_to_print=header+"\t"+motif+"\t"+str(list[header][motif])+"\n"
                file.write(line_to_print)
        file.close()
        return()
    
    # Function to read the counted kmers
            
    def read_kmers(file_name):
        dir={}
        file=open(file_name,"r")
        for line in file:
            line=line.rstrip()
            line=line.split("\t")
            header=line[0]
            motif=line[1]
            count=line[2]
            if header not in dir:
                dir[header]={}
            if motif not in dir[header]:
                dir[header][motif]=count
        file.close()
        return(dir)
    
    # Calculate or read kmers depending if they have been calculated before or not    
        
    if not os.path.isfile("element_kmers.txt"):
        if not os.path.isfile("nonelement_kmers.txt"):
            subprocess.call("sed 's/chr//g' exons_element.bed > exons_element_for_fasta.bed", shell=True)
            subprocess.call("sed 's/chr//g' exons_nonelement.bed > exons_nonelement_for_fasta.bed", shell=True)
            subprocess.call('awk \'{print $1"\t"$2-1"\t"$3+1"\t"$4"\t"$5"\t"$6}\' exons_element_for_fasta.bed > exons_element_for_fasta2.bed', shell=True)
            subprocess.call('awk \'{print $1"\t"$2-1"\t"$3+1"\t"$4"\t"$5"\t"$6}\' exons_nonelement_for_fasta.bed > exons_nonelement_for_fasta2.bed', shell=True)    
            subprocess.call('fastaFromBed -fi '+fasta_file+' -bed exons_element_for_fasta2.bed -fo exons_element.fa -nameOnly', shell=True)
            subprocess.call('fastaFromBed -fi '+fasta_file+' -bed exons_nonelement_for_fasta2.bed -fo exons_nonelement.fa -nameOnly', shell=True)
            pool = mp.Pool(processes=cores)
            results = pool.map(count_kmers, [["exons_nonelement.fa",kmers,2],["exons_element.fa",kmers,1],["exons_element_mutations.fa",kmers,3],["exons_nonelement_mutations.fa",kmers,4]])
        else:
            subprocess.call("sed 's/chr//g' exons_element.bed > exons_element_for_fasta.bed", shell=True)
            subprocess.call('awk \'{print $1"\t"$2-1"\t"$3+1"\t"$4"\t"$5"\t"$6}\' exons_element_for_fasta.bed > exons_element_for_fasta2.bed', shell=True)
            subprocess.call('fastaFromBed -fi '+fasta_file+' -bed exons_element_for_fasta2.bed -fo exons_element.fa -nameOnly', shell=True)
            pool = mp.Pool(processes=cores)
            results = pool.map(count_kmers, [["exons_element.fa",kmers,1],["exons_element_mutations.fa",kmers,3],["exons_nonelement_mutations.fa",kmers,4]])
            introns=read_kmers("nonelement_kmers.txt")
    else:
        if not os.path.isfile("nonelement_kmers.txt"):
                subprocess.call("sed 's/chr//g' exons_nonelement.bed > exons_nonelement_for_fasta.bed", shell=True)
                subprocess.call('awk \'{print $1"\t"$2-1"\t"$3+1"\t"$4"\t"$5"\t"$6}\' exons_nonelement_for_fasta.bed > exons_nonelement_for_fasta2.bed', shell=True)
                subprocess.call('fastaFromBed -fi '+fasta_file+' -bed exons_nonelement_for_fasta2.bed -fo exons_nonelement.fa -nameOnly', shell=True)
                pool = mp.Pool(processes=cores)
                results = pool.map(count_kmers, [["exons_nonelement.fa",kmers,2],["exons_element_mutations.fa",kmers,3],["exons_nonelement_mutations.fa",kmers,4]])
                exons=read_kmers("element_kmers.txt")
        else:
                pool = mp.Pool(processes=cores)
                results = pool.map(count_kmers, [["exons_element_mutations.fa",kmers,3],["exons_nonelement_mutations.fa",kmers,4]])
                exons=read_kmers("element_kmers.txt")
                introns=read_kmers("nonelement_kmers.txt")
                
    for element in results:
        if element[0]==1: elements=element[1]
        if element[0]==2: nonelements=element[1]
        if element[0]==3: element_mutations=element[1]
        if element[0]==4: nonelement_mutations=element[1]
                
    if not os.path.isfile("element_kmers.txt"):
        if not os.path.isfile("nonelement_kmers.txt"):
            pool = mp.Pool(processes=cores)
            variable = pool.map(print_kmers, [[elements,"element_kmers.txt"],[nonelements,"nonelement_kmers.txt"]])
        else:
            print_kmers(elements,"element_kmers.txt")
    else:
        if not os.path.isfile("nonelement_kmers.txt"):
            print_kmers(nonelements,"nonelement_kmers.txt")
            

    print("\tKmers counted: %.0f seconds " % (time.time() - start_time))

    # Combine all the information before subsampling

    file=open("table_kmer_counts.txt","w")

    for gene in genes:
        for motif in kmers:
            if gene not in elements:
                a="0"
                b="0"
            if gene in elements:
                if motif not in elements[gene]:
                    a="0"
                    b="0"
                if motif in elements[gene]:
                    b=str(elements[gene][motif])
                    if gene not in element_mutations:
                        a="0"
                    if gene in element_mutations:
                        if motif in element_mutations[gene]:
                            a=str(element_mutations[gene][motif])
                        if motif not in element_mutations[gene]:
                            a="0"
            if gene not in nonelements:
                c="0"
                d="0"
            if gene in nonelements:
                if motif not in nonelements[gene]:
                    c="0"
                    d="0"
                if motif in nonelements[gene]:
                    d=str(nonelements[gene][motif])
                    if gene not in nonelement_mutations:
                        c="0"
                    if gene in nonelement_mutations:
                        if motif in nonelement_mutations[gene]:
                            c=str(nonelement_mutations[gene][motif])
                        if motif not in nonelement_mutations[gene]:
                            c="0"
            line_to_print=gene+"\t"+motif+"\t"+a+"\t"+b+"\t"+c+"\t"+d+"\n"
            file.write(line_to_print)
                
    file.close()
        
    print("\tAll information combined: %.0f seconds " % (time.time() - start_time))
    subprocess.call('cp table_kmer_counts.txt table_kmer_counts_backup.txt', shell=True)
    subprocess.call('grep -v -P "\.\t" table_kmer_counts_backup.txt > table_kmer_counts.txt', shell=True)
if not os.path.isfile("counts.txt"):
    # Get all the information before subsampling

    file=open("table_kmer_counts.txt","r")

    table_kmer_counts={}
    for line in file:
        line=line.rstrip()
        line=line.split("\t")
        gene=line[0]
        motif=line[1]
        el_mut=int(line[2])
        el_len=int(line[3])
        ne_mut=int(line[4])
        ne_len=int(line[5])
        tot_mut=el_mut+ne_mut
        tot_len=el_len+ne_len
        if gene not in table_kmer_counts:
            table_kmer_counts[gene]={}
        table_kmer_counts[gene][motif]=[el_mut,el_len,ne_mut,ne_len]
    file.close()

    for gene in table_kmer_counts:
        total=0.0
        for motif in table_kmer_counts[gene]:
            total+=table_kmer_counts[gene][motif][1]
        table_kmer_counts[gene]["Total"]=total
    print("\tTable_kmer_count.txt parsed: %.0f seconds " % (time.time() - start_time))

    # Function to subsample the background region

    def define(gene):
        random.seed(seed_number)
        try:
            if "ENSGR" not in gene:
                if gene in table_kmer_counts:
                    gene=str(gene)
                    elements_counts={}
                    nonelements_counts={}
                    eyes=0
                    eno=table_kmer_counts[gene]["Total"]
                    neyes=0
                    neno=0
                    for nuc in kmers:
                        im=int(table_kmer_counts[gene][nuc][2])
                        it=int(table_kmer_counts[gene][nuc][3]-im)
                        if nuc not in nonelements_counts: nonelements_counts[nuc]=[]
                        nonelements_counts[nuc].extend([nuc+"+"]*im)
                        nonelements_counts[nuc].extend([nuc+"-"]*it)
                        eyes+=table_kmer_counts[gene][nuc][0]
                        if nuc not in elements_counts: elements_counts[nuc]=0
                        elements_counts[nuc]=int(round(float(table_kmer_counts[gene][nuc][1])/float(table_kmer_counts[gene]["Total"])*1000.0))

                    for nuc in kmers:
                        random.shuffle(nonelements_counts[nuc])

                    llllllist=[]
                    for nuc in kmers:
                        if elements_counts[nuc]>0:
                            a_number=(len(nonelements_counts[nuc])-(len(nonelements_counts[nuc])%elements_counts[nuc]))/elements_counts[nuc]
                            llllllist.append(a_number)
                    min_number=min(llllllist)
                    
                    for nuc in kmers:
                        if elements_counts[nuc]>0: 
                            
                            tmpNum = elements_counts[nuc]*min_number
                            flist = list(filter(lambda x:'+' in x, nonelements_counts[nuc][:int(tmpNum)]))                     
                            found=len(flist)#no of intr mutations
                            neyes+=found
                            neno+=int(elements_counts[nuc]*min_number)
                    return((gene,eyes,eno,neyes,neno))

        except: dadadsr=0

    # Label the per-element files by the element type (from the element BED
    # filename) rather than the pooled gene id, so these files are
    # self-identifying and cannot be mistaken for single-gene results.
    element_name = os.path.splitext(os.path.basename(element_file))[0]
    file2=open("counts.txt","w")
    pool = mp.Pool(processes=cores)
    results = pool.map(define, genes)    
    for element in results:
        if element is not None:
            line=element_name+"\t"+str(element[1])+"\t"+str(int(element[2]))+"\t"+str(element[3])+"\t"+str(element[4])+"\n"
            file2.write(line)
    print("\tCounts calculated: %.0f seconds " % (time.time() - start_time))               
    file2.close()    
    
    #To make counts_output.tsv and kmer_count files
    file3=open("counts.txt","r")
    al, bl, cla, dl = [], [], [], []
    for line in file3:
        w=line.split('\t')
        al.append(int(w[1]))
        bl.append(int(w[2]))
        cla.append(int(w[3]))
        dl.append(int(w[4].strip('\n')))
    ta, tb, tc,td = sum(al), sum(bl), sum(cla),sum(dl)

    # creating data
    dat = [[ta, tc], [tb-ta, td-tc]]
    odd_ratio, p_value = stats.fisher_exact(dat,alternative='greater')# CHecks the probability of obtaining a table as extreme or more extreme than the observed 
    chistat, chip, dof, expected = stats.chi2_contingency(dat,correction = False)
    
    # Open the TSV file for writing
    with open('counts_output.tsv', 'w') as tsv_file:
            # Create a writer object
            writer = csv.writer(tsv_file, delimiter='\t')

            # Write the values to the TSV file
            writer.writerow(['Element', 'Element Mutations(M)', 'Element Not Mutated(N-M)', 'Non Element Mutations:(m)', 'Non Element not mutated:(n-m)', 'OR-FET','p-value FET', 'chisq-stat','chi p-value'])
            writer.writerow([element_name, ta, tb - ta, tc, td - tc, odd_ratio, p_value, chistat, chip])


    file4=open("table_kmer_counts.txt","r")
    kmera, kmerb, kmerc, kmerd = [], [], [], []
    for line in file4:
        w=line.split('\t')
        kmera.append(int(w[2]))
        kmerb.append(int(w[3]))
        kmerc.append(int(w[4]))
        kmerd.append(int(w[5].strip('\n')))
    kmra, kmrb,kmrc,kmrd = sum(kmera),sum(kmerb),sum(kmerc),sum(kmerd)
    print("Element Mutations(M):",kmra,"Element Not Mutated(N-M):",kmrb-kmra,"Non Element Mutations(m):",kmrc,"Non Element not mutated(n-m):",kmrd-kmrc)
        # creating data
    data = [[kmra, kmrc], [kmrb-kmra, kmrd-kmrc]]
    odd_ratiok, p_valuek = stats.fisher_exact(data,alternative='greater')# checks OR>1
    chistatk, chipk, dofk, expectedk = stats.chi2_contingency(data,correction=False)
    with open('kmer_counts_output.tsv', 'w') as tsv_file:
        # Create a writer object
        writer1 = csv.writer(tsv_file, delimiter='\t')

        writer1.writerow(['Element', 'Element Mutations(M)', 'Element Not Mutated(N-M)', 'Non Element Mutations:(m)', 'Non Element not mutated:(n-m)', 'OR-FET','p-value FET', 'chisq-stat','chi p-value'])
        writer1.writerow([element_name, kmra, kmrb - kmra, kmrc, kmrd - kmrc, odd_ratiok, p_valuek, chistatk, chipk])

# Execute statistics with R            

subprocess.call('Rscript '+input_folder+'Stats.r '+exonic_filter+" "+background_filter+" "+number_of_genomes, shell=True)        
print("\tQvalues calculated: %.0f seconds " % (time.time() - start_time))
