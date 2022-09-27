# Ortholog_Finder

## Usage: 

1.	Change input features in `ortholog_finder.py` (Example input and its corresponding output has been given)
2.	Run `python ortholog_finder.py` for the main function(Runtime is approx.15-20min per organism. It might be faster in non-busy hours and significantly slower in busier hours for Blast). The output will be stored in `orthologs.csv`
3.	Run `python ortholog_finder_tests.py`, runtime is 25sec per ortholog, and your output will be in the terminal. 


## Approach:

Ortholog_Finder uses NCBIWWW to connect to the web-based server, run blastp for the protein, and cut off results with a user-defined threshold. It then picks the one with the highest percent identity to run tblastn, similarly with a user-defined threshold. 

After tblastn, it locates the partial CDS and obtains the protein accession number of the whole CDS. Looking up that accession number would give us the start and end location(and whether it’s reverse) of the whole CDS in the genome sequence. This allows it to download the whole CDS in the genome sequence as well as the 200 sequences up and down stream in the correct order, which can be used for further sequence analysis/manipulation.

The final output is a pandas dataframe formatted to include name of organism, 200bp upstream, cds, 200bp downstream, and genome accession number.

There is also a short test for basic features, such as CDS starting with ATG, ending with a stop codon, and that the input organism corresponds to the output organism. The test file can be modified to include more customized tests. Note that only two organism input formats are supported, Genus+Species and only Genus. 


