#!/usr/bin/env python3 -c

'''
Read NCBI taxonomy dmp files (nodes or names) and holds a dictionary
'''

from .ReadTaxonomy import ReadTaxonomy
from gzip import open as zopen
import zlib
import os
import logging
logger = logging.getLogger(__name__)

class ReadTaxonomyNCBI(ReadTaxonomy):
	"""docstring for ReadTaxonomyNCBI."""
	def __init__(self, taxonomy_file=False, database=False,**kwargs):
		super(ReadTaxonomyNCBI, self).__init__(database=database)
		self.taxonomy_file = taxonomy_file
		self.names_dmp = taxonomy_file.replace("nodes","names")
		self.names = {}
		self.length = 0
		self.ids = 0
		self.accessionfile = False

	def write_missing(self,missing):
		'''Write missing genomes to file'''
		with open("FlexTaxD.not_added", "w") as of:
			for gen in missing:
				print(gen, end="\n", file=of)
		return

	def set_accession_file(self,file):
		self.accessionfile = file

	def parse_taxonomy(self,treefile=False):
		'''Parse taxonomy information'''
		if hasattr(self, "names_dmp"):
			logger.info("Parse names file {}".format(self.names_dmp))
			self.read_names(self.names_dmp)
			self.length = self.database.num_rows("nodes")
		if self.taxonomy_file:
			logger.info("Parse nodes file {}".format(self.taxonomy_file))
			self.read_nodes(self.taxonomy_file)
			self.ids = self.database.num_rows("tree")

	def read_nodes(self, taxfile):
		'''Read NCBI node file and store a dictionary with names in object'''
		with open(taxfile, "r") as _taxfile:
			for taxonomy_row in _taxfile:
				data = taxonomy_row.strip().split("\t|\t")
				child,parent,rank = data[0],data[1],data[2]
				if rank == "None" or rank == None:
					rank = "no rank"
				if rank not in self.rank:
					self.add_rank(rank)
				lev = self.rank[rank]
				self.database.add_link(child=data[0],parent=data[1],rank=lev)
			self.database.commit()
		return

	def read_names(self, taxfile):
		'''Read NCBI node file and store a dictionary with nodes in object'''
		with open(taxfile, "r") as _taxfile:
			for taxonomy_row in _taxfile:
				data = taxonomy_row.strip().split("\t|\t")
				taxid = data[0]
				name = data[1]
				_type = False
				if len(data) > 3:
					_type = data[3].rstrip("|\t")
				if _type == "scientific name" or not _type:
					self.add_node(name, id=taxid)
			self.database.commit()
		return

	def parse_genebank_file(self,filepath,filename):
		logger.debug("Parse file {filename}".format(filename=filename))
		genebankid = filename.split("_",2)
		genebankid = genebankid[0]+"_"+genebankid[1]
		f = zopen(filepath,"r")
		refseqid = f.readline().split(b" ")[0].lstrip(b">")
		f.close()
		self.refseqid_to_GCF[refseqid] = genebankid
		return

	def parse_genomeid2taxid(self, genomes_path,annotation_file):
		'''To allow NCBI databases to be build from scratch the sequences names needs to be stored in the database,
			this function parses the accession2taxid file from NCBI to speed up the function and reduce the amount
			of stored datata only sequences in input genomes_path will be fetched
		'''
		logger.info("Parsing ncbi accession2taxid, genome_path: {dir}".format(dir = genomes_path))
		self.refseqid_to_GCF = {}
		for root, dirs, files in os.walk(genomes_path,followlinks=True):
			for filename in files:
				if filename.strip(".gz").endswith(".fna"):
					filepath = os.path.join(root, filename)
					self.parse_genebank_file(filepath,filename)
		logger.info("genomes folder read, {n} sequence files found".format(n=len(self.refseqid_to_GCF)))
		if not annotation_file.endswith("accession2taxid.gz"):
			raise TypeError("The supplied annotation file does not seem to be the ncbi nucl_gb.accession2taxid.gz")
		annotated_genome = set()
		try:
			with zopen(annotation_file,"r") as f:
				headers = f.readline().split(b"\t")
				for row in f:
					if row.strip() != "": ## If there are trailing empty lines in the file
						if len(row.split(b"\t")) > 2:
							try:
								refseqid,taxid = row.split(b"\t")[1:3]
							except:
								logger.info(row)
								logger.info(row.split(b"\t"))
								if len(annotated_genome) > 0:
									logger.info("Potential error in last row?")
								else:
									logger.info("Error on first line in annotation file, check format!")
							try:
								genebankid = self.refseqid_to_GCF[refseqid]
								self.database.add_genome(genome=genebankid,_id=taxid.decode("utf-8"))
								annotated_genome.add(refseqid)
							except KeyError:
								pass
				self.database.commit()
		except zlib.error as e:
			logger.info("Error in annotation file {e}".format(e=e))
		missing = set(self.refseqid_to_GCF.keys()) - annotated_genome
		missing = [self.refseqid_to_GCF[m] for m in missing] ## Translate to GCF ids
		if logging.root.level <=20: ## Equal to --verbose
			logger.info("Printing non added genome id´s (GCF) to ./FlexTaxD.not_added")
			self.write_missing(missing)
		logger.debug(missing)  ## If debug also print genomes to terminal
		logger.info("Genomes not matching any annotation {len}".format(len=len(missing)))
		return missing
