import unittest
import pandas as pd
from Bio import Entrez, SeqIO
Entrez.email = "hanqiaoz@berkeley.edu"

class TestMethods(unittest.TestCase):

    def setUp(self):
        '''Loads csv file. Enter your csv file name below.'''
        test_file_name = 'orthologs.csv'
        self.csv = pd.read_csv(test_file_name, index_col=[0])

    def test_start(self):
        '''Tests that each CDS starts with a start codon'''
        for index, row in self.csv.iterrows():
            self.assertEqual('ATG', row.cds[:3].upper())

    def test_stopcodon(self):
        '''Tests that each CDS ends with a stop codon'''
        for index, row in self.csv.iterrows():
            self.assertIn(row.cds[-3:].upper(),['TAA','TAG','TGA'])

    def test_organism(self):
        '''Tests that the organism name of ortholog matches out initial input requirement.
        Accepted input organism format : 
            Only genus(eg. 'Citrobacter') or genus + spieces (eg. 'Citrobacter koseri') '''
        for index, row in self.csv.iterrows():
            accession = row.accession
            ortholog_organism = row.organism
            open("test_organism.gb", "w").close()
            #fetch the genbank file given the ortholog accession number
            genbank_handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
            genbank_result = open("test_organism.gb", "w")
            genbank_result.write(genbank_handle.read())
            genbank_result.close()
            genbank_handle.close()
            #record the name of organism from NCBI entry
            for rec in SeqIO.parse("test_organism.gb", "genbank"):
                organism = rec.annotations['organism']
            #compare NCBI organism name with input organism
            if ' '  in ortholog_organism: #Format: Genus + species
                self.assertEqual(organism.lower(), ortholog_organism.lower())
            else: #Format: Genus
                self.assertEqual(organism.split(' ')[0].lower(), ortholog_organism.lower())

if __name__ == '__main__':
    unittest.main()