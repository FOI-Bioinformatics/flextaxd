#!/usr/bin/env python3 -c
'''
Process directory
'''

import logging,os
from .database.DatabaseConnection import DatabaseFunctions
logger = logging.getLogger(__name__)

class ProcessDirectory(object):
	"""ProcessDirectory matches database entries to files on disk
		Potential development is to store annotation to disk location for faster
		processing an option can allow full reprocess of a input directory
	"""

	def __init__(self, database,limit=False):
		super(ProcessDirectory, self).__init__()
		self.database = DatabaseFunctions(database)
		self.genome_id_dict = self.database.get_genomes(self.database , limit=limit)
		self.ref_ext = [".fna"]
		self.oth_ext = [".fasta",".fa"]
		self.ext = self.ref_ext+self.oth_ext
		self.genome_names = []
		self.genome_path_dict = {}
		self.files = []

	def get_genome_names(self):
		'''
		Returns
			list - genome_names of files annotated in the database matching local file
		'''
		return self.genome_names

	def get_genome_path_dict(self):
		'''
		Returns
			list - genome_path_dicts dictionary of genome_id to location on disk
		'''
		return self.genome_path_dict

	def get_files(self):
		'''
		Returns
			list - list of dictionaries with accession and path
		'''
		return self.files

	def get_GCF_taxid(self,genome_name):
		'''If file contains GCF in the start and fna in the end it is most likely
			from refseq try match do genomeid2taxid dictionary
			otherwise the name should match the full name of the sequence file minus fa,fasta,fna

		Parameters
			str    - GCF/GCA
		------
		Returns
			int 	- taxid
		'''
		try:
			taxid = self.genome_id_dict[genome_name.strip()]
			return taxid
		except KeyError:
			return False

	def check_GCF(self,genome_name):
		'''The file was neither a refseq file nor a custom fasta or fa file,
			try if the file is a genbank file (GCA instead of GCF although not 100%)
		Parameters
			str    - GCF/GCA
		------
		Returns
			int 	- taxid
		'''
		try:
			taxid = self.genome_id_dict[genome_name]
			return taxid
		except KeyError:
			return False

	def find_local(self,genome_name,fname):
		'''Does the file have a GCF/GCA start, ends with fna but is still a custom named genome
		Parameters
			str    - GCF/GCA
		------
		Returns
			int 	- taxid
		'''
		#### Swap back namechange if no success
		try:
			taxid = self.genome_id_dict[fname.rsplit(".",1)[0]] ## strip .fa .fasta or .fna from base filename
			return taxid
		except KeyError:
			return False

	def find_local_fasta(self,genome_name,fname):
		'''Does the file have a GCF/GCA start ends with fna but does not have fna in database name and is still a custom named genome
		Parameters
			str    - GCF/GCA
		------
		Returns
			int 	- taxid
		'''
		try:
			taxid = self.genome_id_dict[fname] ## strip .fa .fasta or .fna from base filename
			return taxid
		except KeyError:
			'''This file had no match in the reference folder, perhaps it is not annotated in the database'''
			self.notused.add(genome_name)
			logger.debug("#Warning {gcf} could not be matched to a database entry!".format(gcf=genome_name.strip()))
			return False

	def process_file(self,file,fname,root):
		'''Parameters
			str    - name of file
			str    - path to file location
		------
		Returns
			boolean - true if file existed
			'''

		'''The bulk of genomes is expected to come from official sources therefore search '''
		if (fname.startswith("GCF_") or fname.startswith("GCA_")) and fname.endswith(tuple(self.ref_ext)):
			## Rename file to only GCF_name as this is what is expected to be present in the database for bulk genomes (NCBI/GTDB)
			genome_name = fname.split("_",2)
			genome_name = genome_name[0]+"_"+genome_name[1]
		else:
			genome_name = fname.split(".fna")[0] ## strip .fna from name
		taxid = self.get_GCF_taxid(genome_name)

		'''If the GCF genome is not present, check if the file starts with GCF/GCA but is a custom filename'''
		if not taxid:
			taxid = self.find_local(genome_name,fname)
		'''If the file is still not matching a database entry use the complete name (including .fasta/.fna/.fa)'''
		if not taxid:
			taxid = self.find_local_fasta(genome_name,fname)

		'''In official sources there is sometimes a file called from_genomic.fna make sure this file does not get included in the file list'''
		if not file.strip(".gz").endswith("from_genomic.fna"):
			filepath = os.path.join(root, file)  ## Save the path to the file
			#logger.debug("File added")
			self.files.append(filepath)
			self.genome_names.append(genome_name.strip())
			self.genome_path_dict[genome_name.strip()] = filepath
		return True

	def walk_directory(self,folder_path):
		'''Use os walk to find all existing fasta genomes in path
		Parameters
			str     - path to genomes directory
			dict    - dictionary with database genome annotations
		------
		Returns
			list - list of genomes found in directory structure
			list - list of genome id´s not in the directory structure that can optionally be downloaded
		'''
		count = 0
		self.notused = set()
		download_files = []
		logger.info("Process genome path ({path})".format(path=folder_path))
		for root, dirs, files in os.walk(folder_path,followlinks=True):
			for file in files:
				fname = file.strip(".gz") ## remove gz if present
				if fname.endswith(tuple(self.ext)):
					if count % 1000 == 0:
						print("Processed {count} genomes".format(count=count), end="\r")
					if self.process_file(file,fname,root):
						count +=1
				elif file == "MD5SUMS" or file.endswith(".txt"):
					pass
				else:
					logger.debug("#Warning {gcf} does not have a valid file ending".format(gcf=file))
		logger.info("Processed {count} genomes".format(count=count))
		self.files = list(set(self.files))
		self.genome_names = list(set(self.genome_names))
		return self.files, self.genome_names

	def process_folder(self,folder_path):
		'''Walk through folder and match genomes to database entries, database entries with no matching file be downloaded'''
		logger.info("Number of genomes annotated in database {n}".format(n=len(self.genome_id_dict)))
		self.files, self.genome_names = self.walk_directory(folder_path)
		existing_genomes = []
		download_files = []
		for genome in self.genome_names:
			existing_genomes.append(genome)
		for file_not_present in set(self.genome_id_dict.keys()) - set(existing_genomes):
			download_files.append({"genome_id":file_not_present,"outdir":folder_path.rstrip("/")+"/downloads"})
		return self.files,download_files
