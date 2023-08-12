import sys
import os
import gzip
from collections import defaultdict

#====================================================================#

# now read in the input fastq and split it up:     

def read_fastq_file_and_split(input_fastq_file, seqs_per_output_file, output_file_prefix):

    # open an output file:
    output_file_cnt = 1
    output_file = "%s_%d.fastq.gz" % (output_file_prefix, output_file_cnt)
    outputfileObj = gzip.open(output_file, "wb") # write out the output file in gzipped format
    print("Opening",output_file,"...")

    # read in the input file:
    fileObj = gzip.open(input_fastq_file, "rt") # this opens a gzipped file in text mode
    seqcnt = 0
    linecnt = 0
    for line in fileObj:
        line = line.rstrip()
        if line.startswith('@') and (linecnt == 0 or linecnt == 4): 
            linecnt = 0
            seqcnt += 1
            if seqcnt == (seqs_per_output_file + 1):
                outputfileObj.close()
                output_file_cnt += 1
                output_file = "%s_%d.fastq.gz" % (output_file_prefix, output_file_cnt)
                outputfileObj = gzip.open(output_file, "wb") # write out the output file in gzipped format
                print("Opening",output_file,"...")
                seqcnt = 1
        outputline = "%s\n" % line
        outputfileObj.write(outputline.encode()) # need to write bytes (not a string) to the output file, see https://stackoverflow.com/questions/49286135/write-strings-to-gzip-file
        linecnt += 1 
    fileObj.close()
    outputfileObj.close()

    # the fastq file looks like this:
    # @M03558:259:000000000-BH588:1:1101:15455:1333 2:N:0:NTTGTA
    # NCAAGCATCTCATTTTGTGCATATACCTGGTCTTTCGTCTTCTGGCGTGAAGTCGCCGACTGAATGCCAGCAATCTCTTTTTGAGTCTCATTTTGCATCTCGGCAATCTCTTTCTGATTGTCCAGTTGCATTTTAGTAAGCTCTTTTTGATTCTCAAATCCGGCGTCAACCATACCAGCAGAGGAAGCATCAGCACCAGCACGCTCCCAAGCATTAAGCTCAGGAAATGCAGCAGCAAGA
TAATCACGAGT
    # +
    # #>>ABBFFFFFFGGGGGGGGGGHHHHHHHHHHHGHHGH2FFHHHHHGGGGGHHHGGGGGGGGHHHGHHHHGHHHHHHHHHHHGGGHHHHHHHHHHHHHHHHGGGGGHGFHHHHHHHHHHHHHGHHHHGHFHHHHHEFFFHGFFHHHHHGGHHHHHFGHHHHGGGCGGGFHHGGHGHHHHHHHHAGFFHHHHHGGGAEFHHHFDGDDGGGGEEBFGGGGGFGGGGEFGFFFGGGGGGFFFF
FFBBFFFDFAA

    return 

#====================================================================#

def main():
    
    # check the command-line arguments:
    if len(sys.argv) != 4 or os.path.exists(sys.argv[1]) == False:
        print("Usage: %s input_fastq_file seqs_per_output_file output_file_prefix" % sys.argv[0]) 
        sys.exit(1)
    input_fastq_file = sys.argv[1] # input fastq file                  
    seqs_per_output_file = int(sys.argv[2]) # number of sequences to put into each output file
    output_file_prefix = sys.argv[3] # prefix to use for the output file names

    # now read in the input fastq and split it up:     
    read_fastq_file_and_split(input_fastq_file, seqs_per_output_file, output_file_prefix)
    
    #====================================================================#

if __name__=="__main__":
    main()

#====================================================================#
