codonOverview is a class that provides codon-centric information in regards to a given gene list.
    
The input for this class is is a .fasta file containing only
CODING SEQUENCES (and the different outputs will be incorrect
if any sequence with e.g. UTR or intronic regions are included).

Import requirements:
	operator
	pandas
	biopython
	itertools
	openpyxl

Functions:

	HighestOccurrenceCodon(self, query, n = 10, proportionInsteadOfValue = False):
    	    - query, codon of interest (must be a defined legal codon and must be uppercase letters).
        	- n, top n genes to return.
        	- proportionInsteadOfValue, if True, the top n genes will be for proportions relative to
        	  gene length and not raw count.
        
        	Outputs a list of tuples of gene names and number/proportion of query.

	AverageCodons(self, proportionInsteadOfValue = False):
	        Returns a dictionary with either the number or the proportion of every codon relative to genome.

	PrintSequence(self, query, queryFormat):
	        Takes gene name and outputs either DNA, RNA or protein sequence

	ExportList(self, listOfTuples, outName, outFormat):
	        Takes a list of tuples with gene name and number/proportion
	        and outputs a .csv file containing ID, number/proportion, gene description and gene length

Example workflow:
	# Call class on C. crescentus cds fasta
	codon_overview = codonOverview('../data/Caulobacter_crescentus_na1000.ASM2200v1.cds.all.fa')

	# Get list of all codons with average proportions across entire cds
	allProportions = codon_overview.AverageCodons(proportionInsteadOfValue = True)

	# Get top 10 genes with highest number of ACT codon
	top10genesACT = codon_overview.HighestOccurrenceCodon('ACT', 10, proportionInsteadOfValue = False)

	# Print protein sequence from ACL94628
	proteinSeq_ACL94060 = codon_overview.PrintSequence('ACL94060', 'protein')
	print('Protein sequence for ACL94060:\n'+proteinSeq_ACL94060)
	print('Number of threonin in ACL94060:\n'+str(proteinSeq_ACL94060.count('T'))) # <- Simple way to count all T's

	# Export a list made from the HighestOccurrenceCodon-function into a csv-/excel-file
	codon_overview.ExportList(top10genesACT, '../output/top10ACT', 'excel')