import operator
import pandas as pd
from Bio import SeqIO
from itertools import islice

# legalCodons can be customized to only analyse on desired codons 
legalCodons = [
        'TTT', 'TTC', 'TTA', 'TTG', 'CTT', 
        'CTC', 'CTA', 'CTG', 'ATT', 'ATC', 
        'ATA', 'ATG', 'GTT', 'GTC', 'GTA', 
        'GTG', 'TAT', 'TAC', 'TAA', 'TAG', 
        'CAT', 'CAC', 'CAA', 'CAG', 'AAT', 
        'AAC', 'AAA', 'AAG', 'GAT', 'GAC', 
        'GAA', 'GAG', 'TCT', 'TCC', 'TCA', 
        'TCG', 'CCT', 'CCC', 'CCA', 'CCG', 
        'ACT', 'ACC', 'ACA', 'ACG', 'GCT', 
        'GCC', 'GCA', 'GCG', 'TGT', 'TGC', 
        'TGA', 'TGG', 'CGT', 'CGC', 'CGA', 
        'CGG', 'AGT', 'AGC', 'AGA', 'AGG', 
        'GGT', 'GGC', 'GGA', 'GGG']

emptyCodonDict = dict(zip(legalCodons, [0] * len(legalCodons)))

class codonOverview:
    '''
    A class that provides codon-centric information in regards
    to a given gene list.
    
    The input for this class is is a .fasta file containing only
    CODING SEQUENCES (and the different outputs will be incorrect
    if any sequence with e.g. UTR or intronic regions are included).
    '''
    
    def __init__(self, fastaFile):
        '''The class initializes by parsing the given .fasta file with biopython,
        the parsed fasta file is made into a codon object and a dictionary with
        additional information'''
        self.parsedFasta = SeqIO.parse(fastaFile, 'fasta')
        self.codonObject,self.informationDict = self._countCodons()
    
    def _countCodons(self):
        '''Creates and returns a codon object and informationDict from a fasta
        file parsed by Bioppython's SeqIO module.
        
        Codon object is a dictionary of tuples of dictionaries. The outer dictionary
        having gene names as keys and the two inner dictionaries having codons as keys
        and number of codons and proportion of codons in respect to the length of the
        gene respectivly'''
        
        codonObject = {}
        informationDict = {}
    
        # Useing Biopython to parse .fasta
        for gene in self.parsedFasta:
        
            # For each gene, initialize by saving:
            # Codon dict with all counts = 0 and the sequence,
            # and initializing the total number of codons in a sequence to 0
            tmpNumberCodonDict = emptyCodonDict.copy()
            tmpSequence = gene.seq
            tmpNCodonsInGene = 0

            informationDict[gene.id] = (gene.description,len(gene.seq),tmpSequence)
            
            # Make sure that the sequence length divides by 3
            if len(tmpSequence) % 3 != 0:
                continue
        
            # Iterate over sequence by for each 3 letters in the sequence,
            # adding 1 to the value of the corrosponding key in tmpCodonDict,
            # while also adding 1 to tmpNCodons (total number of codons in sequence),
            # then removing the codon.
            # Repeat untill sequence length = 0
            while len(tmpSequence) > 0:
                tmpCodon = tmpSequence[0:3]
                tmpNumberCodonDict[str(tmpCodon)] += 1
                tmpSequence = tmpSequence[3:]
                tmpNCodonsInGene += 1
            
            tmpProportionCodonDict = tmpNumberCodonDict.copy()
            
            for individualCodon in tmpProportionCodonDict.keys():
                tmpProportionCodonDict[individualCodon] /= tmpNCodonsInGene
            
            # Save as self defined CodonObject:
            # A dictionary with gene id as key, and a tuple of two codon dictionaries,
            # (one with numbers and one with proportions) as values.
            
            codonObject[gene.id] = (tmpNumberCodonDict, tmpProportionCodonDict)
        return(codonObject, informationDict)
    
    def HighestOccurrenceCodon(self, query, n = 10, proportionInsteadOfValue = False):
        '''Takes:
        - query, codon of interest (must be a defined legal codon and must be uppercase letters)
        - n, top n genes to return
        - proportionInsteadOfValue, if True, the top n genes will be for proportions relative to
        gene length and not raw count.
        
        Outputs a list of tuples of gene names and number/proportion of query'''
        
        tmpDict = {}
    
        # Make sure that query is a real codon
        if query in legalCodons:
        
            # For each gene in the codonObject,
            # add the codon count for the query and the gene name
            # to tmpDict as key and value respectivly.
            for geneName, codonTuple in self.codonObject.items():
                
                if proportionInsteadOfValue == True:
                    queryCount = codonTuple[1][query]
                else:
                    queryCount = codonTuple[0][query]
                
                tmpDict[geneName] = queryCount
            
            # Reverse sort tmpDict by codon count
            #tmpDict = dict(sorted(tmpDict.items(), reverse = True))
            resultList = sorted(tmpDict.items(), key=operator.itemgetter(1), reverse = True)
        
            # Return top n values from tmpDict
            #return(list(islice(tmpDict.items(),n)))
            return(resultList[:n])
    
        else:
            print('query: \'' + query + '\' is not a legal codon')
            
    def AverageCodons(self, proportionInsteadOfValue = False):
        '''Returns a dictionary with either the number or the proportion of every codon relative to genome'''
        avgDict = emptyCodonDict.copy()
        count = 0
        
        for geneName, codonTuple in self.codonObject.items():
            
            if proportionInsteadOfValue == True:
                proportionDict = codonTuple[1].copy()
            else:
                proportionDict = codonTuple[0].copy()
            
            for codonName, codonProportion in proportionDict.items():
                avgDict[codonName] += codonProportion
            count += 1
            
        for codonName, codonProportion in avgDict.items():
            
            if proportionInsteadOfValue == True:
                avgDict[codonName] = round(codonProportion/count,4)
            else:
                avgDict[codonName] = round(codonProportion/count,1)
        
        resultList = sorted(avgDict.items(), key=operator.itemgetter(1), reverse = True)
        return(resultList)

    def PrintSequence(self, query, queryFormat):
        '''Takes gene name and outputs either DNA, RNA or protein sequence'''
        
        legalPro = ['PRO','PROTEIN','AMINO']
        
        sequence = self.informationDict[query][2]
        
        if queryFormat.upper() == 'DNA':
            return(sequence)
        elif queryFormat.upper() == 'RNA':
            return(sequence.transcribe())
        elif queryFormat.upper() in legalPro:
            return(sequence.translate())
        else:
            print('ERROR: '+queryFormat.upper()+' is not a valid format')
        
    def ExportList(self, listOfTuples, outName, outFormat):
        '''Takes a list of tuples with gene name and number/proportion
        and outputs a .csv file containing ID, number/proportion, gene description and gene length'''
        nameList = []
        valueList = []
        descriptionList = []
        lengthList = []
                
        for i in range(0,len(listOfTuples)):
            
            name = listOfTuples[i][0]
            value = listOfTuples[i][1]
            
            nameList.append(name)
            valueList.append(value)
            descriptionList.append(self.informationDict[name][0])
            lengthList.append(self.informationDict[name][1])
            
            if isinstance(value, float) == True:
                valueFormat = 'Proportion of codon'
            else:
                valueFormat = 'Number of codon'
            
        tmpDict = {'ID': nameList,
                   valueFormat: valueList,
                   'Description': descriptionList,
                   'Gene length': lengthList}

        df = pd.DataFrame(tmpDict)
        
        if outFormat == 'csv':
            df.to_csv(outName+'.csv', index = False)
        elif outFormat == 'excel':
            df.to_excel(outName+'.xlsx', index = False)