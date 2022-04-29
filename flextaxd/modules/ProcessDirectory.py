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
		#logger.debug(self.genome_id_dict)

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

	def get_taxid(self,genome_name):
		'''If file contains GCF in the start and fna in the end it is most likely
			from refseq try match do genomeid2taxid dictionary
			otherwise the name should match the full name of the sequence file minus fa,fasta,fna

		Parameters
			str    - genome_name
		------
		Returns
			int     - taxid
		'''
		try:
			taxid = self.genome_id_dict[genome_name]
			return taxid
		except KeyError:
			return False

	def is_gcf_gca(self,fname,debug=False):
		'''Paramterers
			str     - File name

		------
		Returns
			str     - GCF name
			boolean - false if not GCF/GCA
		'''
		try:
			GCX,END,REST = fname.split("_",2)  ## If a name contains anything after the GCF number remove this if split by _
			if debug:
				logger.debug("[{} {} {}]".format(GCX,END,REST))
			NUM,version = END.split(".",1)
			if debug:
				logger.debug("[{} {}]".format(NUM,version))
			if GCX.startswith(("GCF","GCA")):                        ## Must start with GCF/GCA
				if len(NUM) == 9 and NUM.isdigit():                    ## All true GCF/GCA names have 9 digits
					if len(version) <= 2 and version.isdigit():      ## version number after . is 1-99 OBS -> Will have to be updated if a genome reach version number higher than 99
						if fname.endswith(tuple(self.ref_ext)):         ## A genome downloaded from refseq or genbank will end with .fna
							genome_name = "{GCX}_{NUM}.{version}".format(GCX=GCX,NUM=NUM,version=version)
							return genome_name
		except:        ## If the above is not true it is not a GCF file name return False
			pass
		return False

	def find_local(self,fname):
		'''Check if file is a custom genome defined without extension
		Parameters
			str    - filename
		------
		Returns
			int     - taxid
		'''
		fname = fname.rsplit(".",1)[0]  ## If a name contains anything after the GCF number remove this if split by _
		taxid = self.get_taxid(fname)
		return taxid,fname

	def find_local_fasta(self,fname):
		'''Check if file is a custom genome defined with extension
		Parameters
			str    - filename
		------
		Returns
			int     - taxid
		'''
		taxid = self.get_taxid(fname)
		if not taxid:
			'''This file had no match in the reference folder, perhaps it is not annotated in the database'''
			self.is_gcf_gca(fname,True)
			self.notused.add(fname)
			logger.debug("#Warning {gcf} could not be matched to a database entry!".format(gcf=fname.strip()))
		return taxid,fname

	def process_file(self,file,fname,root,taxid=False):
		'''Parameters
			str    - name of file
			str    - path to file location
		------
		Returns
			boolean - true if file was processed
			'''
		'''The bulk of genomes is expected to come from official sources'''
		genome_name = self.is_gcf_gca(fname)
		if genome_name:
			taxid = self.get_taxid(genome_name)
		'''If the file is not a GCF or GCA file check if the file starts with GCF/GCA but is a still a custom filename'''
		if not taxid:
			taxid,genome_name = self.find_local(fname)
		'''If the file is still not matching a database entry use the complete name (including .fasta/.fna/.fa)'''
		if not taxid:
			taxid,genome_name = self.find_local_fasta(fname)
		'''In official sources there is sometimes a file called from_genomic.fna; make sure this file does not get included in the file list'''
		if not file.strip(".gz").endswith("from_genomic.fna") and taxid:
			filepath = os.path.join(root, file)  ## Save the path to the file
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
		if not folder_path:
			raise IOError("Parameter --genomes_path was not set".format(folder_path))
		logger.info("Process genome path ({path})".format(path=folder_path))
		logger.debug("Extensions valid: ({f})".format(f=self.ext))
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
		if len(self.genome_names) == 0:
			logger.info("# WARNING: No genomes from the input folder was added! Are your genomes annotated properly in the database?")
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
