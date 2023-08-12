import sys
import os
from collections import defaultdict

#====================================================================#

# define a function to read in the genes that are in families, for our species of interest:

def find_genes_in_families(families_file, our_species_name, locus_tag):
    """read in the genes that are in families, for our species of interest """

    # define a set to store the genes in families:
    genes_in_families = set()

    # read in the input file of genes in families:
    fileObj = open(families_file, "r")
    for line in fileObj:
        line = line.rstrip()
        if line.startswith('family'): # eg. family 5168 : SPAL_0001096400.1:mRNA (strongyloides_papillosus) Bm11175 (brugia_malayi) HCOI_1164600.1 (haemonchus_contortus)...
            temp = line.split()
            family = temp[1] # eg. 4376
            temp = line.split(': ')
            line = temp[1] # SPAL_0000608500.1:mRNA (strongyloides_papillosus) Bm10870f (brugia_malayi) F23B12.1 (caenorhabditis_elegans)...
            temp = line.split()
            # make a list of tuples, with gene names and species names
            # taken from http://stackoverflow.com/questions/2631189/python-every-other-element-idiom
            # uses the zip function: http://docs.python.org/2/library/functions.html#zip
            # see also on zip: http://nbviewer.ipython.org/github/pycam/python-intro/blob/release/Introduction_to_python_session_2.ipynb
            # note the array 'pairs' seems to be empty after we have looked at it once, this must be a feature of 'zip'.
            set_of_species = set() # define the set of species we have seen
            family_genes = set() # define the set of genes in the family
            pairs = zip(temp[0::2], temp[1::2])
            for pair in pairs:
                (gene, species) = pair # unpack the tuple, eg. gives gene= SPAL_0000608500.1:mRNA species= (strongyloides_papillosus)
                species = species[1:-1] # remove the parentheses
                # get the gene name:
                # eg. SRAE_X000112100.t1:mRNA/SSTP_0000225500.1:mRNA/SPAL_0001294400.1:mRNA/SVEN_0373900.1/PTRK_0001428200.1:mRNA/RSKR_0000682500.1:mRNA
                if ':mRNA' in gene:
                    temp2 = gene.split(':mRNA') # eg. SRAE_X000112100.t1:mRNA
                    gene = temp2[0] # eg. SRAE_X000112100.t1
                if '-mRNA-' in gene: # eg. ASIM_0000655901-mRNA-1
                    temp2 = gene.split('-mRNA-')
                    gene = temp2[0] # eg. ASIM_0000655901
                if '_' in gene and '.' in gene: # eg. SRAE_X000112100.t1 or MhA1_Contig1285.frz3.fgene2
                    temp2 = gene.split('.')
                    afterdot = temp2[1] # eg. t1 or fgene2
                    if (afterdot.startswith('t') or afterdot.isdigit() == True) and species != 'caenorhabditis_elegans':
                        gene = temp2[0] # eg. SRAE_X000112100                 
                if '.t' in gene: # eg. HCOI02162100.t1
                    temp2 = gene.split('.t')
                    afterdot = temp2[1] # eg. 1
                    if afterdot.isdigit() == True:
                        gene = temp2[0] # eg. HCOI02162100 
                if species == our_species_name:
                    if locus_tag == 'NOO' or locus_tag == 'NLS' or locus_tag == 'NAV':
                                              gene = get_new_gene_name(gene, locus_tag)
                    genes_in_families.add(gene)
    fileObj.close()

    return genes_in_families

#====================================================================#

# get the gene name from a line of the embl file:

def get_gene_name2(line):
    """ get the name of a gene from line of an embl file
    >>> get_gene_name2('FT                   /note="ID:gene:NAV_0000008934"')
    'NAV_0000008934'
    >>> get_gene_name2('FT                   /note="ID:cds:SMUV_0000034201-mRNA-1"')
    'SMUV_0000034201'
    """

    temp = line.split()
    gene = temp[1] # /note="ID:gene:NAV_0000008934"
    temp = gene.split("\"")
    gene = temp[1] # ID:gene:NAV_0000008934 or ID:cds:SMUV_0000034201-mRNA-1
    temp = gene.split(":")
    gene = temp[2] # NAV_0000008934 or SMUV_0000034201-mRNA-1
    temp = gene.split("-mRNA-")
    gene = temp[0] # NAV_0000008934 or SMUV_0000034201

    return gene

#====================================================================#

# read in the input embl file, and write out a new embl file:

def read_input_embl_and_write_output_embl(input_embl, output_embl, dodgy_genes, genes_in_families):

    # open an output file:
    outputfileObj = open(output_embl,"w")

    # read in the input file:
    fileObj = open(input_embl, "r")
    for line in fileObj:
        line = line.rstrip()
        if 'source:WormBase_imported' in line:
            pass # don't print this line to the output
        else:
            if 'AC *' in line: # AC * _ALUE_contig0000001 length=91883
                output_line = make_embl_output_line1(line)
                output_line = "%s\n" % output_line
                outputfileObj.write(output_line)
            elif '/note="ID:gene:' in line or 'note="ID:cds:' in line or 'note="ID:exon:' in line : # FT                   /note="ID:gene:NAV_0000008934"
                                                                                                  # FT                   /note="ID:cds:SMUV_0000034201-mRNA-1"
                                                                                                  # FT                   /note="ID:exon:SMUV_0000587401-mRNA-1.1"
                gene = get_gene_name2(line) # eg. NAV_0000008934
                output_line = "%s\n" % line
                outputfileObj.write(output_line) # write out the current line
                                if gene in dodgy_genes: 
                    output_line = "FT                   /pseudo\n"
                    outputfileObj.write(output_line)
                    output_line = "FT                   /note=\"odd gene structure probably caused by assembly error\"\n" 
                    outputfileObj.write(output_line)
                    if gene in genes_in_families:
                        output_line = "FT                   /note=\"predicted protein belongs to a gene family in the Compara database of the International Helminth Genomes Consortium\"\n"
                        outputfileObj.write(output_line)
            else: 
                output_line = "%s\n" % line
                outputfileObj.write(output_line)
    fileObj.close() 
    outputfileObj.close()

    return

#====================================================================#

# make a new output line with the accession:

def make_embl_output_line1(line):
    """make a new output line with the accession
    >>> make_embl_output_line1('AC * _ALUE_contig0000001 length=91883')
    'AC * _ALUE_contig0000001'
    """

    temp = line.split()
    scaffold = temp[2] # _ALUE_contig0000001 
    scaffold = scaffold[1:] # ALUE_contig0000001 
    output_line = "AC * _%s" % scaffold

    return output_line

#====================================================================#

# define a function to trim off '0's at the start of a scaffold name: 

def trim_off_zeroes(scaffold_name):
    """ trim off zeroes at the start of a scaffold name:
    >>> trim_off_zeroes("00001")
    '1'
    >>> trim_off_zeroes("00010")
    '10'
    """

    # trim off zeroes at the start:
    finished = False
    while finished == False:
        if scaffold_name.startswith('0'):
            scaffold_name = scaffold_name[1:]
        else:
            finished = True

    return scaffold_name

#====================================================================#

# get a gene name from a transcript name:

def get_gene_name_from_transcript_name(transcript):
    """get gene name from a transcript name
    >>> get_gene_name_from_transcript_name('nAv.1.0.1.t01019-RA')
    'nAv.1.0.1.g01019'
    """

    if '-RA' in transcript:
        temp = transcript.split('-RA')
        transcript = temp[0] # eg. nAv.1.0.1.t01019
    if '.t' in transcript: # nAv.1.0.1.t01019 
        temp = transcript.split('.t')
        transcript = temp[0] + ".g" + temp[1] # nAv.1.0.1.g01019
    gene = transcript

    return gene

#====================================================================#

# read set of genes/transcripts that look dodgy:

def read_set_of_dodgy_genes(dodgy_gene_list, locus_tag):

    dodgy_genes = set()
  
    fileObj = open(dodgy_gene_list, "r")
    for line in fileObj:
        line = line.rstrip()
        temp = line.split()
        transcript = temp[0] # eg. nAv.1.0.1.t01019-RA
        gene = get_gene_name_from_transcript_name(transcript) # eg. nAv.1.0.1.g01019
        # get the new name for the gene:
        if locus_tag == 'NOO' or locus_tag == 'NLS' or locus_tag == 'NAV':
            gene = get_new_gene_name(gene, locus_tag)
        dodgy_genes.add(gene)
    fileObj.close()

    return dodgy_genes

#====================================================================#

# find the new name for a gene:

def get_new_gene_name(old_gene_name, locus_tag):
    """find the new name for a gene
    >>> get_new_gene_name('nAv.1.0.1.g01019', 'NOO')
    'NOO_0000001019'
    """

    temp = old_gene_name.split('.g')
    gene_name = temp[1] # 08934
    gene_name = trim_off_zeroes(gene_name) # 8934 

    #====================================================================#

# find the new name for a gene:

def get_new_gene_name(old_gene_name, locus_tag):
    """find the new name for a gene
    >>> get_new_gene_name('nAv.1.0.1.g01019', 'NOO')
    'NOO_0000001019'
    """

    temp = old_gene_name.split('.g')
    gene_name = temp[1] # 08934
    gene_name = trim_off_zeroes(gene_name) # 8934 
    gene_name_len = len(gene_name)
    num_zeroes_to_add = 10 - gene_name_len
    zeroes_to_add = "0" * num_zeroes_to_add
    new_gene_name = "%s_%s%s" % (locus_tag, zeroes_to_add, gene_name)

    return new_gene_name

#====================================================================#

def main():
    
    # check the command-line arguments:
    if len(sys.argv) != 7 or os.path.exists(sys.argv[1]) == False or os.path.exists(sys.argv[3]) == False or os.path.exists(sys.argv[6]) == False:
        print("Usage: %s input_embl output_embl dodgy_gene_list locus_tag our_species families_file" % sys.argv[0]) 
        sys.exit(1)
    input_embl = sys.argv[1] # e.g. onchocerca_ochengi.embl
    output_embl = sys.argv[2] # e.g. onchocerca_ochengi_new.embl 
    dodgy_gene_list = sys.argv[3] # list of genes/transcripts that look dodgy eg. onchocerca_ochengi.proteins_with_stops2 
    locus_tag = sys.argv[4] # eg. NOO 
    our_species = sys.argv[5] # eg. onchocerca_ochengi
    families_file = sys.argv[6] # /nfs/helminths02/analysis/50HGP/00ANALYSES/final_families/complete_families.txt_27mar2017 

    # read in the set of genes/transcripts that look dodgy:
    dodgy_genes = read_set_of_dodgy_genes(dodgy_gene_list, locus_tag)

    # read in the genes that are in families, for our species of interest:
    genes_in_families = find_genes_in_families(families_file, our_species, locus_tag)

    # read in the input embl file, and write out a new embl file:
    read_input_embl_and_write_output_embl(input_embl, output_embl, dodgy_genes, genes_in_families)

    print("FINISHED\n")

#====================================================================#

if __name__=="__main__":
    main()

#====================================================================#
