#!/usr/bin/env python3 -c
'''
Download source genomes for building databases
'''

import logging
logger = logging.getLogger(__name__)

from modules.functions import download_genomes
from multiprocessing import Process,Manager,Pool

class DownloadGenomes(object):
	"""DownloadGenomes takes a list of genome names (GCF or GCA) and download (if possible) files from NCBI refseq and or genbank"""

	def __init__(self,processes=20,outdir="./",force=False):
		super(DownloadGenomes, self).__init__()
		self.not_downloaded = []        ## Place holder for files not possible to download
		self.download_map   = []        ## Place holder for files in subprocess
		self.processes = processes      ## Number of simultanous downloads allowed
		self.genome_names = []
		self.genome_path = {}
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

	def run(self, files):
		'''Download list of GCF and or GCA files from NCBI

		Parameters
			list - list of GCF/GCA ids

		Returns
			list - list of files not downloaded
		'''
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
		count = 0
		while True:
			if added.qsize() == 0: ## Check if no genome was succesfully downlaoded if no break
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
		if len(self.not_downloaded) > 0:
			self.write_missing(self.not_downloaded)
		return self.not_downloaded
