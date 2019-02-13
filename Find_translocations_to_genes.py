import argparse
import glob
import os
import subprocess as sub
import sys
from collections import Counter
import itertools
from os import path

import pandas as pd


class Tlx():
    ''' Processes .tlx files ''' 
    def __init__(self, file_name):

        self.file_name = file_name
        self.name = os.path.splitext(os.path.basename(self.file_name))[0]
        self.bed = None
        
        self.bed_file = f'Bed_files/{self.name}.bed'
        self.fasta_file = f'Fasta_files/{self.name}.fasta'
        self.rss_file = f'RSS_files/{self.name}.rss.txt'

        try:
            os.mkdir('Bed_files')
            os.mkdir('Fasta_files')
            os.mkdir('RSS_files')
            os.remove('Excluded_Regions.bed')
        except:
            pass
   
    def make_bed(self, exclusions = None, window = 0, full = False):
        
        ''' Generates a .bed file from Tlx files. 
            The full argument can be used if the entire mapped translocation is required''' 

        def form_window(l, window, full= False):
            if not full:       
                if l[4] == '1':
                    return [l[2], str(int(l[5]) - window) , str(int(l[5]) + window), l[0], '+']
                else:
                    return [l[2], str(int(l[6]) - window) , str(int(l[6]) + window), l[0], '-']
            else:
                return [l[2], l[5], l[6], l[0], l[4]]               
       
        def form_exclusion_bed(exclusions):
            ''' makes .bed file with excluded regions '''
            
            bed = 'Excluded_Regions.bed'
            with open(bed, 'w') as w:
                for ex in exclusions:
                    c = [ex.split(':')[0], 
                         ex.split(':')[1].split('-')[0],
                         ex.split(':')[1].split('-')[1],]           
                    w.write('\t'.join(c) + '\n')
            return bed
        
        def write_bed():
            w = open(self.bed_file, 'w')
            for l in self.bed:
                w.write(l + '\n')

        with open(self.file_name) as r:
            tlx_coords = [x.split('\t') for x in r][1:]
            tlx_windows = ['\t'.join(form_window(x, window, full = full)) for x in tlx_coords]
            self.bed = '\n'.join(tlx_windows)
        
        if exclusions:
            ex = form_exclusion_bed(exclusions)
            cmd = f'bedtools subtract -a - -b {ex}'
            bedtools_out = sub.check_output(cmd.split(), input = str.encode(self.bed))
            self.bed = bytes.decode(bedtools_out).split('\n')[:-1]

        write_bed()
        
        return self.bed
          
    def make_fasta(self, genome):
        if self.bed:
            fasta = self.name + '.fasta'
            cmd = f'bedtools getfasta -fi {genome} -bed {self.bed_file} -fo {self.fasta_file}'
            sub.run(cmd.split())  
        else:
            raise ValueError('No bed file specified, make one first')

class gene_parser():
    ''' Extracts required genes from a .gtf file ''' 
    def __init__(self, genes, gtf):

        self.gene_file = genes
        self.gene_names = []
        self.gtf = gtf

        self.genes_required = {x.strip() for x in open(self.gene_file)}

    def get_gene_coords(self):
        ''' Takes the longest protein coding transcript for each gene in the list'''

        def parse_gtf_record(r):
            gene = {}
            gene['name'] = r.split('gene_name')[1].split('\"')[1]
            
            r = r.split('\t')
            gene['chr'] = 'chr' + r[0]
            gene['start'] = int(r[3])
            gene['end'] = int(r[4])
            return gene  
                
        if self.gene_file:
            grep_cmd  = ['grep','-f', self.gene_file, self.gtf]
            rec =  bytes.decode(sub.check_output(grep_cmd)).split('\n')[:-1]
            rec = (x for x in rec if 'transcript' in x.split('\t')[2] and 'protein_coding' in x)
            
            get_gene_name = lambda x: x.split('gene_name')[1].split('\"')[1]
            rec = [x for x in rec if get_gene_name(x) in self.genes_required]
        else:
            raise ValueError('No genes specified')

        rec.sort(key= lambda x: int(x.split('\t')[4]) - int(x.split('\t')[3]), reverse=True)
             
        seen = set()
        genes = []
        for l in rec:
            n = l.split('gene_id')[1].split('\"')[1]
            if not n in seen:
                genes.append(l)
                seen.add(n)
        
        return [parse_gtf_record(x) for x in genes]

def count_overlaps(genes, bed):
    ''' Counts the number of translocations overlapping the extracted genes '''

    gb = '\n'.join(['\t'.join([g['chr'], str(g['start']), str(g['end']), g['name']]) for g in genes])
    bedtools_cmd = f'bedtools intersect -a - -b {bed} -c'

    bed_out = bytes.decode(sub.check_output(bedtools_cmd.split(), input = str.encode(gb))).split('\n')
    s = [x.split('\t') for x in bed_out][:-1]
    hits = [(x[-2], int(x[-1])) for x in s]
    return hits

def display_data(data, totals,  normalise = True):

    ''' Groups data by experiment, summarises counts and will normalise 
        the counts obtained if normalise is true'''

    def group_data(data):
        ''' Groups the supplied data by locus and experiment''' 

        #key_func = lambda x: (x.split('_')[0], x.split('_')[1])
        key_func =  lambda x: x.split('_')[1]
        data = sorted(iter(data), key = key_func)
        for k, g in itertools.groupby(data, key = key_func):
            yield list(g) 

    data_df = pd.DataFrame()
    for group in group_data(data):
        experiments = [data[x] for x in group]
        total = sum([totals[x] for x in group])

        series =  pd.DataFrame(experiments).sum(axis=0)

        if normalise:
            series = series.apply(lambda x: (x / total) * (1*10**5))
        
        data_df = pd.concat([data_df, series], axis=1, sort = True)
    
    get_name = lambda x: x[0].split('_')[0] + '_' + x[0].split('_')[1]
    col = [get_name(x) for x in group_data(data)]
    data_df.columns = col
    data_df.to_csv('Translocations_to_genes.csv')

def main():

    def get_args():
        parser = argparse.ArgumentParser()
        parser.add_argument('genes', help= 'File with gene names')
        parser.add_argument('gtf', help = 'GTF file to parse')
        parser.add_argument('-n', '--normalise', help = 'Flag to normalise counts based on number of translocations',
                            action = 'store_true')
        return parser.parse_args()

    
    args = get_args()

    ''' Parse a .gtf file for the required gene names and returns the longest protein coding transcript'''
    gp = gene_parser(args.genes, args.gtf)
    genes = gp.get_gene_coords()
    
    ''' Generate a .bed file from the .tlx file and overlap with the required genes obtained from the supplied .gtf file
        the number of translocations in each .bed file is also stored for normalisation of the overlaps obtained. ''' 
    data = {}
    total = {}
    for tlx_file in sorted(glob.glob('Tlx_files/*.tlx')):
        
        tlx = Tlx(tlx_file)
        bed = tlx.make_bed(exclusions=('chr6:74225755-74235755',), full= True)
        genes_hit = count_overlaps(genes, tlx.bed_file)
        data[tlx.name] = dict(genes_hit)
        total[tlx.name] = len(bed)
    
    display_data(data, total, normalise = args.normalise)


    





            




            


  
            






    
    
    
    





    

if __name__ == '__main__':
    main()
