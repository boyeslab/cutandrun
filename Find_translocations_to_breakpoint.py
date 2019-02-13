import argparse
import glob
import itertools
import os
import subprocess as sub

import pandas as pd

from pybedtools import BedTool


class Tlx():

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

def check_overlaps(a, b, window = 0, count = False):

    a = BedTool(a)
    b = BedTool(b)

    return a.window(b, w=window, c=count)

def make_grouped_dataframe(data):

    ''' Makes a DataFrame from translocations overlapping patient breakpoints,
        grouped by locus and experiment'''

    def group_data(data):
        
        ''' Groups the supplied data by locus and experiment''' 

        key_func = lambda x: (x.split('_')[0], x.split('_')[1])
        data = sorted(iter(data), key = key_func)
        for k, g in itertools.groupby(data, key = key_func):
            yield list(g) 

    for group in group_data(data):
        values = (data.get(x) for x in group)
        grouped_data = filter(None, itertools.chain.from_iterable(values))
        df = pd.DataFrame([list(x) for x in grouped_data])
        if not df.empty:
            df.columns = ['chr_trans','start_trans', 'end_trans', 'trans_ID', 'strand', 'chr_brk', 'start_brk', 'end_brk', 'brk_ID']
            df.name = f'{group[0].split("_")[0]}_{group[0].split("_")[1]}.csv'
            yield df

def main():
    def get_args():
        p = argparse.ArgumentParser()
        p.add_argument('breakpoints', help = 'Breakpoints file in .bed format')
        p.add_argument('-ow', '--overlap_window', type = int, help = 'Sets the window between translocations and breakpoints', default = 50)
        return p.parse_args()
    

    if not os.path.exists('Output_tables'):
        os.mkdir('Output_tables')

    ''' Remove previous run files '''
    for f in glob.glob('Output_tables/*.csv'):
        os.remove(f)

    args = get_args()   
    print(f'\nFinding patient breakpoints within {args.overlap_window} bp of the translocation breakpoint')
    print(f'Finding motifs within {args.motif_window} bp of the translocation\n')
     
    ''' Make a bed file from the Tlx and check for overlaps with patient breakpoints'''
    data = {}
    for f in glob.glob('Tlx_files/*.tlx'):
        tlx = Tlx(f)
        tlx.make_bed(exclusions=('chr6:74225755-74235755',))
        data[tlx.name] = check_overlaps(tlx.bed_file, args.breakpoints, args.overlap_window)
           
    ''' Generate .bed from the translocations overlapping patient breakpoints '''
    for df in make_grouped_dataframe(data):
        bed = df.iloc[:,:4]
        bed.to_csv('tmp.bed', sep='\t', header=False, index=False)
      
    ''' Remove temporary files '''
    for f in glob.glob('*.bed'):
        os.remove(f)

if __name__ == '__main__':
    main()
