FASTA file format
-----------------
In bioinformatics, FASTA format is a text-based format for representing DNA sequences, in which base pairs are represented using a single-letter code [A,C,G,T,N] where A=Adenosine, C=Cytosine, G=Guanine, T=Thymidine and N= any of A,C,G,T. The format also allows for sequence names and comments to precede the sequences.
A sequence in FASTA format begins with a single-line identifier description, followed by lines of DNA sequence data. The identifier description line is distinguished from the sequence data by a greater-than ('>') symbol in the first column. The word following the ">" symbol is the identifier of the sequence, and the rest of the line is a description (optional) separated from the identifier by a white space or tab. The sequence data starts on the next line following the text line and ends if another line starting with a ">" appears; this indicates the start of another sequence.

For DNA sequences from one species :
The sequence identifier typically refers to the target gene regulated by the hidden motif in the dataset:

        >target_gene1     description(optional)
        lines of {ACGTNacgtn}
        >target_gene2     description(optional)
        ....
        >target_geneEnd     description(optional)
        lines of {ACGTNacgtn}

For orthologous DNA sequences from multiple species (input for PHMS and NOMS):
Sets of orthologous sequences (regulating the same target genes in different species) are separated by a line consisting of a double greater-than ('>>') symbol in the first column. The word following the '>>' symbol is further called the group identifier, optionally followed by a tab-spaced more detailed description. The group('>>') identifier typically refers to the target gene regulated by the hidden motif in the dataset. The following sequence('>') identifiers now refer to the species genome of the respective orthologous sequences:

        >>target_gene1     description(optional)
        >species1     description(optional)
        lines of {ACGTNacgtn}
        >species2     description(optional)
        ....
        >>target_gene2     description(optional)
        >species1     description(optional)
        lines of {ACGTNacgtn}
        >species2     description(optional)
        ....
        >>target_geneEnd     description(optional)
        ....
        >speciesEnd     description(optional)
        lines of {ACGTNacgtn}
