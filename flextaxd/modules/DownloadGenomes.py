#!/usr/bin/env python3 -c
'''
Download source genomes for building databases
'''

import logging
logger = logging.getLogger(__name__)

from modules.functions import download_genomes
from multiprocessing import Process,Manager,Pool
from subprocess import Popen,PIPE
import os,time

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

	def download_represenatives(self,genome_path,url="https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/bac120_marker_genes_reps.tar.gz"):
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
				os.remove("{file}".format(file=input_file_name))
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
