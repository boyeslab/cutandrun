import argparse
import glob
import itertools
import os
import socket
import subprocess as sub
import time
import urllib.request
import zipfile

import pandas as pd
import psutil
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.firefox.firefox_profile import FirefoxProfile
from selenium.webdriver.remote.command import Command
from selenium.webdriver.support.ui import Select


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

class RSSsite_interface():
    ''' Interfaces with the RSSsite to identify RSSs in .fasta files '''
    def __init__(self, fasta, output, rss_type):
        self.fasta = fasta
        self.output = output
        self.rss_type = {'12':'12 bases', '23': '23 bases'}.get(rss_type)
    
    def communicate(self):
        ''' Uploads the generated .fasta file to RSSsite and analyses the sequences for potential RSSs '''

        driver = webdriver.Firefox()
        driver.get('http://www.itb.cnr.it/rss/analyze.html')

        select_species_option = Select(driver.find_element_by_name("species"))
        select_species_option.select_by_visible_text('human')

        select_rss_option = Select(driver.find_element_by_name('spacer'))
        select_rss_option.select_by_visible_text(self.rss_type)

        upload_file = driver.find_element_by_name('upfile')
        upload_file.send_keys(self.fasta)

        analyse_sequence_button = driver.find_element_by_xpath('/html/body/table[2]/tbody/tr/td[3]/table/tbody/tr[1]/td/table/tbody/tr[8]/td[2]/input')
        analyse_sequence_button.click()

        time.sleep(5)
        cpu = psutil.cpu_percent()
        while cpu > 20:
            cpu = psutil.cpu_percent()
        
        download_link = driver.find_element_by_link_text('Click here to get the zipped txt tab separated version of this table')
        download_address = download_link.get_attribute('href')

        urllib.request.urlretrieve(download_address, self.output)
        driver.close()

        file_to_extract = zipfile.ZipFile(self.output)
        file_to_extract.extractall()
        files = glob.glob('file*.txt')
        rss_file = files.pop(0)
        os.rename(rss_file, self.output)

        for f in files:
            os.remove(f)
        
        return rss_file

def count_rss(rss_file):
    '''Counts RSSs with a RIC score above the threshold '''
    def remove_fails(rss_file):
        with open(rss_file) as reader:
            rss = [x for x in reader if not 'FAIL' in x.split('\t')[-1]]
            rss = [x for x in rss if x[:3] == 'chr']
        return rss

    rss_passes = remove_fails(rss_file)
    rss_unique = {x.split('\t')[0] for x in rss_passes}
    return len(rss_unique)

def display_data(data, rss_type):
    ''' Groups and desplays data by locus/experiment '''
    def group_data(data):
        ''' Groups the supplied data by locus and experiment''' 

        key_func = lambda x: (x.split('_')[0], x.split('_')[1])
        #key_func =  lambda x: x.split('_')[1]
        data = sorted(iter(data), key = key_func)
        for k, g in itertools.groupby(data, key = key_func):
            yield list(g) 
    
    data_df = pd.DataFrame()
    for group in group_data(data):
        exp = [data[x] for x in group]
        df = pd.DataFrame(exp).sum(axis=0)
        data_df = pd.concat([data_df, df], axis=1, sort = True)
    
    get_name = lambda x: x[0].split('_')[0] + '_' + x[0].split('_')[1]
    col = [get_name(x) for x in group_data(data)]
    data_df.columns = col
    data_df.to_csv(f'Translocations_to_{rss_type}_RSSs')
        
def main():
    def get_args():
        parser = argparse.ArgumentParser()
        parser.add_argument('-w','--window', default = 45, type = int)
        parser.add_argument('-RSS', '--RSS_type', default=12)
        parser.add_argument('-g', '--genome', help = 'Genome in fasta format')
        return parser.parse_args()

    args = get_args()
    experiments = {}
    for tlx_file in sorted(glob.glob("*.tlx")):

        tlx = Tlx(tlx_file)
        tlx.make_bed(exclusions=('chr6:74225755-74235755',), window = args.window)
        tlx.make_fasta(args.genome)

        interface = RSSsite_interface(tlx.fasta_file, tlx.rss_file, args.RSS_type)
        interface.communicate()

        rss_count = count_rss(tlx.rss_file)
        total = len(tlx.bed)

        experiment = {f'{args.RSS_type} RSS': rss_count, 'Not_RSS': (total - rss_count), 'Total': total}
        experiments[tlx.name] = experiment

    display_data(experiments, args.RSS_type)
    
if __name__ == '__main__':
    main()
