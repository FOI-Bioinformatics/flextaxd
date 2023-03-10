#!/usr/bin/env python3 -c
'''
Download source genomes for building databases
'''

import logging
logger = logging.getLogger(__name__)

from modules.functions import download_genomes
from multiprocessing import Process,Manager,Pool
from subprocess import Popen,PIPE,STDOUT
import os,time,zipfile,re,shutil

class DownloadGenomes(object):
    """DownloadGenomes takes a list of genome names (GCF or GCA) and download (if possible) files from NCBI refseq and or genbank"""

    def __init__(self,processes=20,outdir="./",force=False,download_path="./"):
        super(DownloadGenomes, self).__init__()
        self.not_downloaded = []        ## Place holder for files not possible to download
        self.download_map   = []        ## Place holder for files in subprocess
        self.processes = processes      ## Number of simultanous downloads allowed
        self.genome_names = []
        self.genome_path = {}
        self.download_path = download_path
        self.outdir = outdir
        self.force=force
        if self.force:
            logger.info("Download genomes even if given GCF is not available WARNING: this may attempt to download withdrawn genome assemblies!")
        max = 50
        if self.processes > max:
            logger.warning("Maximum {p} parallel downloads allowed!".format(p=max))
            self.processes = 50         ## Limit the number of paralell processes to maximum NCBI connections

    def get_genome_names(self):
        '''
        Returns
            list - genome_names of files annotated in the database matching local file
        '''
        return self.genome_names

    def get_genome_path(self):
        '''
        Returns
            list - genome_paths dictionary of genome_id to location on disk
        '''
        return self.genome_path

    def _split(self,a, n):
        k, m = divmod(len(a), n)
        return list(a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))

    def write_missing(self,missing):
        '''Write missing genomes to file'''
        with open("{outdir}/FlexTaxD.missing".format(outdir=self.outdir), "w") as of:
            for gen in missing:
                print(gen["genome_id"], end="\n", file=of)
        return

    def download_represenatives(self,genome_path,url="https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/gtdb_genomes_reps.tar.gz"):
        '''Specific function for GTDB download using their tar files instead of downloading from NCBI
            Parameters
                str - url
                str - path to genome directory
            Returns
                path
        '''
        logger.debug(genome_path)
        input_file_name = url.split("/")[-1]
        self.location = genome_path.rstrip("/")+"/representatives"
        self.representative_file = url.rsplit("/")[-1]
        '''Check if exists already if so ask to replace or continue without downloading'''
        if os.path.exists(self.location+"/"+input_file_name):
            ans = input("A represenative file already exist, (u)se file, (o)verwrite (c)ancel? (u o,c): ")
            if ans in ["o", "O"]:
                logger.info("Overwrite current progress")
                os.remove(self.location+"/"+"{file}".format(file=input_file_name))
            elif ans.strip() in ["u", "U"]:
                logger.info("Resume database build")
                return input_file_name
            else:
                exit("Cancel execution!")
        logger.info("Download represenative genomes from {url}".format(url=url))
        #as suggested by @SilentGhost the `self.location` and `url` should be separate argument
        args = ['wget', '-nd', '-r', '-l', '1', '-p', '-P', self.location, url]
        logger.debug(" ".join(args))
        logger.info("Waiting for download process to finish (this may take a while)!")
        p = Popen(args, stdout=PIPE)
        (output, err) = p.communicate()
        p_status = p.wait()
        return input_file_name

    def parse_representatives(self,downloaded_file):
        '''Parse representative genomes
            untar the representative genomes directory to access the source genomes
            Parameters
                str - path to tar file
            Returns
                path - path to unzipped folder
        '''
        args = ["tar", "-xf", downloaded_file]
        logger.info("Untar file {f}".format(f=downloaded_file))
        logger.debug(" ".join(args))
        p = Popen(args, stdout=PIPE,cwd=self.location)
        (output, err) = p.communicate()
        p_status = p.wait()
        return self.location

    def download_files(self,files):
        '''Download list of GCF and or GCA files from NCBI

        Parameters
            list - list of GCF/GCA ids

        Returns
            list - list of files not downloaded
        '''
        try:
            logger.info("Downloading {files} files".format(files=len(files)))
            if len(files) < self.processes:
                self.processes = len(files)
            self.download_map = self._split(files,self.processes)
            logger.info("Using {np} parallel processes to download files".format(np=self.processes))

            '''function to run download of genomes in paralell'''
            jobs = []
            manager = Manager()
            added = manager.Queue()
            fpath = manager.Queue()
            missing = manager.Queue()
            for i in range(self.processes):
                p = Process(target=download_genomes, args=(self.download_map[i],added,fpath,missing,self.force))
                p.daemon=True
                p.start()
                jobs.append(p)
            for job in jobs:
                job.join()
            self.added = added.qsize()
            if not any(proc.is_alive() for proc in jobs):
                logger.info('All download processes completed')
            elif len(added.qsize()) % 100 == 0:
                print("{n} downloaded".format(added.qsize()),end="\r")

            count = 0
            while True:
                if self.added == 0: ## Check if no genome was succesfully downlaoded if no break
                    logger.info("None of listed genomes could be downloaded! Files not downloaded will be printed to {outdir}/FlexTaxD.missing".format(outdir=self.outdir.rstrip("/")))
                    self.write_missing(files)
                    break
                l = len(self.genome_names)
                gen = added.get()
                if gen:
                    self.genome_names += gen
                    self.genome_path[gen] = fpath.get()
                    if l == len(self.genome_names):
                        break
                else:
                    self.not_downloaded += missing.get()
                    count+=1
                if added.qsize():
                    break
            if len(self.not_downloaded) > 0:
                self.write_missing(self.not_downloaded)
        except KeyboardInterrupt:
            logger.info("Program was interrupted by user: cleaning up subprocesses!")
        finally: ## Make sure all sub-processes are ended even if program is forced to quit
            if any(proc.is_alive() for proc in jobs):
                for p in jobs:
                    print(p)
                    p.kill()
            time.sleep(1)
            if any(proc.is_alive() for proc in jobs):
                logger.error("Could not stop all subprocesses check your process manager and end them manually!")
        return self.not_downloaded

    def download_from_file(self,inputfile,exclude="",all=False):
        '''Download files with ncbis download script using an accession input file'''
        if not all:
            exclude = " ".join([exclude,"--exclude-gff3 --exclude-protein --exclude-rna --exclude-genomic-cds"])
        args = ["genome", "accession" ,"--inputfile", inputfile, exclude]
        p = Popen(args, stdout=PIPE,cwd=self.location)
        (output, err) = p.communicate()
        p_status = p.wait()
        return self.location

    def download_genomes(self,genomes_list):
        ### TERMINOLOGY:
            # self.outdir here is the temporary directory
            # self.download_path here is the --genomes_path directory
            
        ## Write missing genomes to file (to be imported by NCBI datasets software)
        logger.info('Writing missing genomes, N='+str(len(genomes_list)))
        self.write_missing(genomes_list)
        ##/
        ## Compile CMD
        ncbi_datasets_cmd = ['datasets','download','genome','accession',
                            '--inputfile',self.outdir+'/'+'FlexTaxD.missing',
                            '--filename',self.outdir+'/'+'ncbi_download.zip']
        logger.info('Attempting download using NCBI datasets software')
        try:
            ncbi_datasets_process = Popen(ncbi_datasets_cmd,stdout=PIPE,stderr=STDOUT)
            while ncbi_datasets_process.poll() is None:
                line_raw = ncbi_datasets_process.stdout.readline() # bufsize?
                line = line_raw.decode()
                line = line.strip('\n')
                if not line: continue
            
                # write status update lines with carriage return
                if line.startswith('Downloading') or (line.find('Collecting') != -1 and line.find('records') != -1):
                    print('[NCBI-DATASETS] ' + line)
                else:
                    print('[NCBI-DATASETS] ' + line)
        except:
            print('Failed to download from NCBI using "datasets"-software. You can try to run the command manually:')
            print(' '.join(ncbi_datasets_cmd))
            logger.info('NCBI datasets download failed!')
            logger.info('CMD: '+' '.join(ncbi_datasets_cmd))
            exit('Terminated.')
        ##/
        
        ## Unpack NCBI ZIP
        logger.info('Attempting to unzip NCBI zip-file')
        with zipfile.ZipFile(self.outdir+'/'+'ncbi_download.zip','r') as zip_fo:
            zip_fo.extractall(self.outdir+'/'+'ncbi_download')
            
        logger.info('Unpacking genomes')
        for path,dirnames,filenames in os.walk(self.outdir+'/'+'ncbi_download'):
            for dirname in dirnames:
                # Get file from directory
                genome_fa_path = None
                for file_ in os.listdir(path+'/'+dirname):
                    if file_.endswith('.fna'):
                        genome_fa_path = path+'/'+dirname+'/'+file_
                        break
                #/
                # Move file to "genomes"-directory (rename to GCx_123456789.v, where X is A/F and v is a digit 1-9)
                if genome_fa_path != None:
                    genome_fa_new_name = None
                    regex_pattern = r".*GC[A|F]_\d{9}\.\d*" # match file on format GCx_123456789.1
                    match = re.search(regex_pattern,os.path.basename(genome_fa_path))
                    if match:
                        matched_string = match.group()
                        stripped_string = re.sub(f".*({regex_pattern}).*", r"\1", matched_string)
                        genome_fa_new_name = stripped_string+'.fna' # should be formatted as GCX_123456789.1
                    
                    if genome_fa_new_name:
                        os.rename(genome_fa_path,self.download_path+'/'+genome_fa_new_name) # move the file
                        os.system('gzip '+self.download_path+'/'+genome_fa_new_name) # gzip the file
                    else:
                        logger.debug('WARNING: was unable to move genome for '+str(dirname)+', attempted new path:'+str(genome_fa_new_name))
                else:
                    if not dirname in ('data','ncbi_dataset',): # skip known "trigger"-folders in the ZIP-file
                        logger.debug('WARNING: was unable to move genome for '+str(dirname)+', attempted genome path:'+str(genome_fa_path))
                #/
        ##/
        ## Clean up download-files
        logger.info('Clearing download files')
        os.remove(self.outdir+'/'+'ncbi_download.zip')
        shutil.rmtree(self.outdir+'/'+'ncbi_download')
        ##/

    def run(self, files, representative=False,url="https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/gtdb_genomes_reps.tar.gz"):
        '''Download list of GCF and or GCA files from NCBI or download represenative genomes

        Parameters
            list - list of GCF/GCA ids
            bool - download representative genomes instead of list
            url - e.g., https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/gtdb_genomes_reps.tar.gz

        Returns
            list - list of files not downloaded
        '''
        if representative:
            logger.info("Download GTDB Representative genomes")
            tarfile = self.download_represenatives(genome_path=self.download_path, url=url)
            folder_path = self.parse_representatives(tarfile)
            '''Process the downloaded folder'''
            return folder_path,[]
        else:
            if len(files) > 0:
                logger.info("Download missing genomes")
                not_downloaded = self.download_files(files)
                return False,not_downloaded
