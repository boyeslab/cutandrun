import os
import glob
import subprocess as sub
import pandas as pd
import itertools
import matplotlib.pyplot as plt

class Tlx():

    def __init__(self, file_name):

        self.file_name = file_name
        self.name = os.path.splitext(os.path.basename(self.file_name))[0]
        self.bed = None
        self.bed_file = self.name + '.bed'
        self.fasta_file = self.name + '.fasta'
        self.rss_file = self.name + '.rss.txt'

        try:
            os.remove(self.name + '.bed')
            os.remove(self.name + '.fasta')
            os.remove(self.name + '.rss.txt')
        except:
            pass

   
    def make_bed_window(self,window, exclusions):
        def form_window(l, window):
            if l[4] == '1':
                tlx_coord = [l[2], str(int(l[5]) - window) , str(int(l[5]) + window), l[0], l[4]]
            else:
                tlx_coord = [l[2], str(int(l[6]) - window) , str(int(l[6]) + window), l[0], l[4]]
            
            return tlx_coord
        
        def form_exclusion_bed(exclusions):
            bed = 'Excluded_Regions.bed'
            with open(bed, 'w') as w:
                for ex in exclusions:
                    c = [ex.split(':')[0], 
                         ex.split(':')[1].split('-')[0],
                         ex.split(':')[1].split('-')[1],]           
                    w.write('\t'.join(c) + '\n')
            return bed

        print(f'Reading TLX file and generating .bed file with a {window} window')


        with open(self.file_name) as r:
            tlx_coords = [x.split('\t') for x in r][1:]
            tlx_windows = ['\t'.join(form_window(x, window)) for x in tlx_coords]
            tlx_bed_uf = '\n'.join(tlx_windows)
        
        ex = form_exclusion_bed(exclusions)
        cmd = f'bedtools subtract -a - -b {ex}'
        bedtools_out = sub.check_output(cmd.split(), input = str.encode(tlx_bed_uf))
        self.bed = bytes.decode(bedtools_out).split('\n')[:-1]

        self.write_bed()  

        return self.bed
    
    def make_bed_full_transloc(self, exclusions):
        def get_bed_coords(l):
            return [l[2], l[5], l[6], l[0]]
        
        def form_exclusion_bed(exclusions):
            bed = 'Excluded_Regions.bed'
            with open(bed, 'w') as w:
                for ex in exclusions:
                    c = [ex.split(':')[0], 
                         ex.split(':')[1].split('-')[0],
                         ex.split(':')[1].split('-')[1],]           
                    w.write('\t'.join(c) + '\n')
            return bed


        print(f'Reading TLX file and generating .bed from the full translocation')


        with open(self.file_name) as r:
            tlx_coords = [x.split('\t') for x in r][1:]
            tlx_windows = ['\t'.join(get_bed_coords(x)) for x in tlx_coords]
            tlx_bed_uf = '\n'.join(tlx_windows)
        
        ex = form_exclusion_bed(exclusions)
        cmd = f'bedtools subtract -a - -b {ex}'
        bedtools_out = sub.check_output(cmd.split(), input = str.encode(tlx_bed_uf))
        self.bed = bytes.decode(bedtools_out).split('\n')[:-1]

        self.write_bed()

        return self.bed  
        
    def write_bed(self):
        if self.bed:
            w = open(self.bed_file, 'w')
            for l in self.bed:
                w.write(l + '\n')
    
    def make_fasta(self, genome):
        if self.bed:
            fasta = self.name + '.fasta'
            cmd = f'bedtools getfasta -fi {genome} -bed {self.bed_file} -fo {self.fasta_file}'
            sub.run(cmd.split())
         
        else:
            print('No bed file specified, make one first')


def get_chr_no(chromosome):

    chr_no = chromosome[3:]

    if chr_no == 'X':
        return 23
    elif chr_no == 'Y':
        return 24
    elif chr_no == 'M':
        return 25
    else:
        try:
            return int(chr_no)
        except:
            return 26

def get_chr_dist(bed):
    
    
    chrom_l = [x.split('\t')[0] for x in bed]
    chrom_l.sort(key = lambda x: get_chr_no(x))
    
    chrom_s = sorted(set(chrom_l), key = lambda x: get_chr_no(x))
    c = [0 for x in chrom_s]
    
    chr_d = {k:v for k,v in zip(chrom_s, c)}

    for chrom in chrom_l:
        chr_d[chrom] += 1
    
    chrom_d = {k:v for k, v in chr_d.items() if not 'gl' in k}
    
    return chrom_d
    
def plot_pie_chart(counts):
    
    for name, exp in counts.items():
        l = [*exp]
        c = exp.values()

        plt.pie(c, labels= l)
        plt.axis('equal')
        plt.tight_layout()
        plt.savefig(name + '.png')
        plt.clf() 

    k = list(counts.keys())
    k.sort(key = lambda x: (x.split('_')[0], x.split('_')[1]))
    data_grouped = [list(g) for k, g in itertools.groupby(k, lambda x: (x.split('_')[0], x.split('_')[1]))]
    
    data_df = pd.DataFrame()
    for group in data_grouped:

        exp = [counts[x] for x in group]
        df = pd.DataFrame(exp).sum(axis=0)
        data_df = pd.concat([data_df, df], axis=1)
    
    get_name = lambda x: x[0].split('_')[0] + '_' + x[0].split('_')[1]
    col = [get_name(x) for x in data_grouped]
    data_df.columns = col
    print(data_df)
    data_df.to_csv('Chromosome_distribution.csv')
    
    


def main():


    chromosome_counts = {}
    for f in glob.glob('*.tlx'):

        tlx = Tlx(f)
        bed = tlx.make_bed_full_transloc(['chr6:74225755-74235755'])
        chromosome_counts[tlx.name] = get_chr_dist(bed)

    
    chromosome_counts_df = pd.DataFrame(chromosome_counts, columns=chromosome_counts.keys())
    chromosome_counts_df.sort_index(axis=1, inplace=True)
    chromosome_counts_df.to_csv('chromosome_distribution.csv')

      
    
    plot_pie_chart(chromosome_counts)

    

if __name__ == '__main__':
    main()


    
    














