#!/usr/bin/env python
# coding: utf-8

from Bio.Seq import Seq
from Bio.Blast import NCBIWWW,NCBIXML
from Bio import Entrez, SeqIO
import pandas as pd
import re
#change your email for NCBI Entrez requests
Entrez.email = "hanqiaoz@berkeley.edu"

#Enter your cds, list of organisms to search, and customized blastp and tblastn cutoff threshold
cds = "atggactttccgcagcaactcgaagcctgcgttaagcaggccaaccaggcgctgagccgttttatcgccccactgccctttcagaacactcccgtggtcgaaaccatgcagtatggcgcattattaggtggtaagcgcctgcgacctttcctggtttatgccaccggtcatatgttcggcgttagcacaaacacgctggacgcacccgctgccgccgttgagtgtatccacgcttactcattaattcatgatgatttaccggcaatggatgatgacgatctgcgtcgcggtttgccaacctgccatgtgaagtttggcgaagcaaacgcgattctcgctggcgacgctttacaaacgctggcgttctcgattttaagcgatgccgatatgccggaagtgtcggaccgcgacagaatttcgatgatttctgaactggcgagcgccagtggtattgccggaatgtgcggtggtcaggcattagatttagacgcggaaggcaaacacgtacctctggacgcgcttgagcgtattcatcgtcataaaaccggcgcattgattcgcgccgccgttcgccttggtgcattaagcgccggagataaaggacgtcgtgctctgccggtactcgacaagtatgcagagagcatcggccttgccttccaggttcaggatgacatcctggatgtggtgggagatactgcaacgttgggaaaacgccagggtgccgaccagcaacttggtaaaagtacctaccctgcacttctgggtcttgagcaagcccggaagaaagcccgggatctgatcgacgatgcccgtcagtcgctgaaacaactggctgaacagtcactcgatacctcggcactggaagcgctagcggactacatcatccagcgtaataaataa"
organisms = ['Pseudomonas aeruginosa','Pseudomonas putida','Bacillus megaterium','Sporosarcina']
blastp_perc_id_cutoff=0.4
tblastn_perc_id_cutoff=0.95

def main():
    #Create the dataframe for final output 
    df_final = pd.DataFrame(columns = ['organism','pre_200','cds','post_200','accession'])
    for organism in organisms:
        # Add [ORGN] to fit the NCBIWWW input format
        organism_search = organism + '[ORGN]' 
        try:
            #Run blastp, output a XML file with all hits, and parse it    
            blastp_xml_path = generate_blastp_xml(cds,organism_search,blastp_perc_id_cutoff)
            df_blastp = parse_blastp(blastp_xml_path)
        except ValueError:
            #If this organism doesn't work, update the csv and move on. 
            print('No protein sequences above percent identity cutoff! Skipping to new organism...')
            df_final = df_final.append({'organism':organism,'pre_200':'N/A','cds':'N/A','post_200':'N/A','accession' :'N/A'},ignore_index=True)
            continue
        '''A for loop through every blastp hit. In case the first blastp hit doesn't yield a result
        (which is very rare), it automatically uses the next hit. '''
        for index,row in df_blastp.iterrows():
            try:
                # Obtain protein sequence of the hit and run tblastn, parse result. 
                protein_seq = (get_protein_sequences(row.accession))
                tblastn_xml_path = generate_tblastn_xml(protein_seq,organism_search,tblastn_perc_id_cutoff)
                df_tblastn = parse_tblastn(tblastn_xml_path)
                '''A for loop through every tblastn hit. In case a tblastn hit ends up not working
                (for example, if it does not correspond to a cds sequence), use next hit'''
                for index,row in df_tblastn.iterrows():
                    #A bit of a hard code here. I want to make sure the hit is long enough to not be a gene sequence
                    if row.start > 1000:
                        try:
                            #Obtain the accession number for the ortholog's protein sequence.
                            ortholog_accession = get_ortholog_accession(row.accession,row.start,row.end)
                            #Accession number for the genome's sequence
                            genome_accession = row.accession
                            print('ortholog accession number collected')
                            break
                        except ValueError:
                            print("No CDS.. trying next hit")
                #Obtain the start and end bases of cds in the genome sequence, also mark whether it is reversed.
                start,end,reverse = get_complete_cds_ortholog(ortholog_accession)
                #Clean data and prepare them for the dataframe entry 
                list_seq = finish_dataframe(genome_accession, start, end,reverse)
                df_final = df_final.append({'organism':organism,'pre_200':list_seq[0],'cds':list_seq[1],'post_200':list_seq[2],'accession' :genome_accession},ignore_index=True)
                break
            except ValueError:
                print('For this protein, no sequences above cutoff. Skipping to new protein...')
            except NameError:
                print('No genome sequence corresponds to out protein ortholog, skipping to new organism...')
                break        
    df_final.to_csv('orthologs.csv')

def generate_blastp_xml(cds,organism,blastp_perc_id_cutoff=0.4):
    '''Run blastp and output a XML file.'''
    #Erases the file before each run.
    open("my_blastp_query.xml", "w").close()
    protein = Seq(cds).translate()
    result_handle = NCBIWWW.qblast("blastp", "nr",protein,entrez_query = organism,perc_ident =blastp_perc_id_cutoff)
    blast_result = open("my_blastp_query.xml", "w")
    blast_result.write(result_handle.read())
    blast_result.close()
    result_handle.close()
    print('blastp completed')
    return "my_blastp_query.xml"

def parse_blastp(blastp_xml_path):
    '''Parse the resulting XML file into a Dataframe with blastp results'''
    blastp_result = open(blastp_xml_path, 'r')
    blast_records = NCBIXML.read(blastp_result)
    df = pd.DataFrame(columns = ['accession','length','perc_identity'])
    for alignment in blast_records.alignments:
        for hsp in alignment.hsps:
            #Use regex to extract accession number from the hsp in blast records
            m = re.search('(?<=\|)(.*?)(?=\|)', alignment.hit_id).group(1)
            df = df.append({'accession':m,
                       'length':alignment.length,
                       'perc_identity':hsp.identities/ hsp.align_length},ignore_index=True)
    if df.size >= 1:
        # Only returns dataframe if it's not empty
        print('Parsed blastp result')
        return df.sort_values(by='perc_identity',ascending=False)
    else:
        raise ValueError('No protein sequences above percent identity cutoff!')

def get_protein_sequences(accession):
    '''Fetch protein sequence from an accession number'''
    handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="fasta")
    record = handle.read()
    print('Fetched protein sequence')
    return record.split('\n')[1]

def generate_tblastn_xml(protein_seq,organism,tblastn_perc_id_cutoff=0.95):
    '''Run tblastn and outputs a XML file'''
    open("my_tblastn_query.xml", "w").close()
    t_blastn_result_handle = NCBIWWW.qblast("tblastn", "nt",protein_seq,entrez_query = organism,perc_ident=tblastn_perc_id_cutoff)
    t_blastn_blast_result = open("my_tblastn_query.xml", "w")
    t_blastn_blast_result.write(t_blastn_result_handle.read())
    t_blastn_blast_result.close()
    t_blastn_result_handle.close()
    print('tblastn completed')
    return "my_tblastn_query.xml"

def parse_tblastn(tblastn_xml_path):
    '''Parse the resulting XML to a dataframe'''
    tblastn_result = open(tblastn_xml_path, 'r')
    tblastn_records = NCBIXML.read(tblastn_result)
    df_tblastn = pd.DataFrame(columns = ['accession','length','perc_identity','start','end'])
    for alignment in tblastn_records.alignments:
        for hsp in alignment.hsps:
            #Extract the accession number
            m = alignment.hit_id.rsplit('|',1)[0].rsplit('|',1)[1]
            df_tblastn = df_tblastn.append({'accession':m,
                       'length':hsp.align_length,
                       'perc_identity':hsp.identities/ hsp.align_length,
                        'start':hsp.sbjct_start,
                        'end':hsp.sbjct_end},ignore_index=True
                    )
    if df_tblastn.size >=1:
        print('parsed tblastn result')
        return df_tblastn
    else:
        raise ValueError('No cds sequence above percent identity cutoff!')

def get_ortholog_accession(seq_id,start,end):
    '''Given a genome sequence and start/end location, fetch the accession number
    of the CDS's protein translation in the part'''
    open("my_genbank.gb", "w").close()
    genbank_handle = Entrez.efetch(db="nuccore", id=seq_id, rettype="gb", retmode="text",seq_start=start,seq_stop=end)
    genbank_result = open("my_genbank.gb", "w")
    genbank_result.write(genbank_handle.read())
    genbank_result.close()
    genbank_handle.close()
    for rec in SeqIO.parse("my_genbank.gb", "genbank"):
        if rec.features:
            for feature in rec.features:
                if feature.type == "CDS":
                    #The first CDS is the one we need
                    return feature.qualifiers["protein_id"][0]
                    
    raise ValueError('No CDS found in respective location of genome')

def get_complete_cds_ortholog(genbank_acc):
    '''Given accession number of protein, return its location in the genome sequence'''
    open("my_ortholog.gb", "w").close()
    genbank_handle = Entrez.efetch(db="protein", id=genbank_acc, rettype="gb", retmode="text")
    genbank_result = open("my_ortholog.gb", "w")
    genbank_result.write(genbank_handle.read())
    genbank_result.close()
    genbank_handle.close()
    for rec in SeqIO.parse("my_ortholog.gb", "genbank"):
        if rec.features:
            for feature in rec.features:
                if feature.type == "CDS":
                    #"complement" means the CDS is reversed
                    if(feature.qualifiers['coded_by'][0]).startswith('complement'):
                        reverse = True
                    else:
                        reverse = False
                    # Text manipulations to extract the start and end base location
                    temp = feature.qualifiers['coded_by'][0].replace("(","").replace(")","")
                    start = temp.split(':')[1].split('..')[0]
                    end = temp.split(':')[1].split('..')[1]
                    print('start end location of genome cds reached')
                    return int(start),int(end),reverse

def finish_dataframe(genome_accession, start,end,reverse):
    '''Output the cds and pre/post 200 based on whether it's reversed'''
    open("final_cds", "w").close()
    if reverse:
        final_post_200 = Entrez.efetch(db="nuccore", id=genome_accession, rettype="fasta", retmode="text",seq_start=start-200,seq_stop=start-1)
        final_cds = Entrez.efetch(db="nuccore", id=genome_accession, rettype="fasta", retmode="text",seq_start=start,seq_stop=end)
        final_prev200 = Entrez.efetch(db="nuccore", id=genome_accession, rettype="fasta", retmode="text",seq_start=end+1,seq_stop=end+200)
    else:
        final_prev200 = Entrez.efetch(db="nuccore", id=genome_accession, rettype="fasta", retmode="text",seq_start=start-200,seq_stop=start-1)
        final_cds = Entrez.efetch(db="nuccore", id=genome_accession, rettype="fasta", retmode="text",seq_start=start,seq_stop=end)
        final_post_200 = Entrez.efetch(db="nuccore", id=genome_accession, rettype="fasta", retmode="text",seq_start=end+1,seq_stop=end+200)
    final_result = open("final_cds", "w")
    final_result.write(final_prev200.read())
    final_result.write(final_cds.read())
    final_result.write(final_post_200.read())
    final_result.close()
    final_cds.close()
    list_seq = []
    for seq in SeqIO.parse("final_cds", "fasta"):
        if reverse:
            list_seq.append(str(seq.seq.reverse_complement()))
        else:
            list_seq.append(str(seq.seq))
    print('finished seq list')
    return list_seq

if __name__=="__main__":
   main()



