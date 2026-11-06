#!/usr/bin/env python3 -c

'''
Read NCBI taxonomy dmp files (nodes or names) and holds a dictionary
'''

from .ReadTaxonomy import ReadTaxonomy
from . import functions
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
		self.force_multisource = kwargs["force_multisource"]

	def write_missing(self,missing):
		'''Write missing genomes to file, nothing is written when nothing is missing'''
		if not missing:
			return False
		with open("FlexTaxD.not_added", "w") as of:
			for gen in missing:
				print(gen, end="\n", file=of)
		return "FlexTaxD.not_added"

	def sequence_id(self,header):
		'''The id of a sequence is the first white space separated word of its FASTA header, the
			header may hold no description and end with the line ending'''
		words = header.split()
		if not words:
			return False
		return words[0].lstrip(b">")

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

	def parse_genebank_file(self,filepath,filename,reference="refseq"):
		logger.debug("Parse file {filename}".format(filename=filename))
		genebankid = filename.split("_",2)
		genebankid = genebankid[0]+"_"+genebankid[1]
		'''The sequence names are compared as bytes, so the file is read binary, compressed or not'''
		with self.zopen(filepath,"rb") as f:
			refseqid = self.sequence_id(f.readline())
		if not refseqid:
			logger.warning("No sequence id in the first line of {file}, the genome stays unannotated".format(file=filename))
			return False
		self.refseqid_to_GCF[refseqid] = genebankid
		self.refseqid_reference[refseqid] = reference
		return

	def parse_nt_file(self,filepath,filename,reference="nt"):
		'''Parse large file with all nt entries'''
		logger.info("Processing {file}".format(file=filename))
		with self.zopen(filepath,"rb") as f:
			for row in f:
				if row.startswith(b">"):
					seqid = self.sequence_id(row)
					if not seqid:
						continue
					self.refseqid_to_GCF[seqid] = seqid.decode("utf-8")
					self.refseqid_reference[seqid] = reference
					#self.database.add_genome(genome=seqid,_id=taxid.decode("utf-8"),reference="nt")
		return

	def get_annotation_columns(self,headers):
		'''Find the columns holding the sequence id and the taxid in an NCBI accession2taxid file,
			the header of the file decides the layout:
				accession | accession.version | taxid | gi
				accessions | accessions.versioned | accessions.nonredundant_GI_seqVersion | taxid
				gi | accession.version | taxid
			The versioned accession column is the one that matches the sequence names in the genome
			files, the gi column must never be read as a taxid.
		------
		Returns
			tuple - (accession column, taxid column)
		'''
		accession_col,taxid_col = False,False
		for col,name in enumerate([header.strip().lower() for header in headers]):
			if name.endswith(b"taxid"):
				taxid_col = col
			elif (name.startswith(b"accession") and b"version" in name and b"nonredundant" not in name
					and accession_col is False):
				accession_col = col
		if accession_col is False:
			for col,name in enumerate([header.strip().lower() for header in headers]):
				if name.startswith(b"accession"):
					accession_col = col
					break
		if taxid_col is False or accession_col is False:
			raise TypeError("Could not find the accession and taxid columns in the header of the annotation file: {headers}".format(headers=b"\t".join(headers)))
		return (accession_col,taxid_col)

	def parse_genomeid2taxid(self, genomes_path,annotation_file,reference="refseq"):
		'''To allow NCBI databases to be build from scratch the sequences names needs to be stored in the database,
			this function parses the accession2taxid file from NCBI to speed up the function and reduce the amount
			of stored datata only sequences in input genomes_path will be fetched
		'''
		file_endings = (".fna",".fa",".fasta")
		logger.info("Parsing ncbi accession2taxid, genome_path: {dir}".format(dir = genomes_path))
		self.refseqid_to_GCF = {}
		self.refseqid_reference = {}
		for root, dirs, files in os.walk(genomes_path,followlinks=True):
			for filename in files:
				base_name = functions.strip_gz(filename)
				if base_name.endswith(file_endings):
					filepath = os.path.join(root, filename)
					if filename.startswith("GC") and base_name.endswith(".fna"):
						self.parse_genebank_file(filepath,filename,reference=reference)
					elif filename.startswith("nt") or self.force_multisource:
						if self.force_multisource:
							self.parse_nt_file(filepath,filename,reference="nt")
						else:
							if filename.endswith(".gz"): 
								compressed=True 
							else: 
								compressed=False
							stats = os.stat(filepath)
							logger.info("Parsing nt archive, filesize: {filesize}Mb, compressed: {compressed}".format(compressed=compressed,filesize=(stats.st_size / (1024 * 1024))))
							self.parse_nt_file(filepath,filename,reference="nt")
					else:
						self.parse_genebank_file(filepath,filename,reference=reference)
		logger.info("genomes folder read, {n} sequence files found".format(n=len(self.refseqid_to_GCF)))
		if not annotation_file.endswith("accession2taxid.gz"):
			raise TypeError("The supplied annotation file does not seem to be the ncbi nucl_gb.accession2taxid.gz")
		annotated_genome = set()
		try:
			with self.zopen(annotation_file,"rb") as f:
				headers = f.readline().split(b"\t")
				accession_col,taxid_col = self.get_annotation_columns(headers)
				logger.info("Annotation file columns: accession {acc}, taxid {tax}".format(acc=headers[accession_col].strip().decode("utf-8"),tax=headers[taxid_col].strip().decode("utf-8")))
				for row in f:
					if row.strip() == b"": ## If there are trailing empty lines in the file
						continue
					columns = row.rstrip(b"\n").split(b"\t")
					if len(columns) <= max(accession_col,taxid_col):
						logger.info(row)
						logger.info(columns)
						if len(annotated_genome) > 0:
							logger.info("Potential error in last row?")
						else:
							logger.info("Error on first line in annotation file, check format!")
						continue
					refseqid = columns[accession_col].strip()
					taxid = columns[taxid_col].strip()
					try:
						genebankid = self.refseqid_to_GCF[refseqid]
						self.database.add_genome(genome=genebankid,_id=taxid.decode("utf-8"),reference=self.refseqid_reference.get(refseqid,reference))
						annotated_genome.add(refseqid)
					except KeyError:
						pass
				self.database.commit()
		except zlib.error as e:
			logger.info("Error in annotation file {e}".format(e=e))
		missing = set(self.refseqid_to_GCF.keys()) - annotated_genome
		missing = [self.refseqid_to_GCF[m] for m in missing] ## Translate to GCF ids
		if logging.root.level <=20 and missing: ## Equal to --verbose
			logger.info("Printing non added genome id´s (GCF) to ./FlexTaxD.not_added")
			self.write_missing(missing)
		logger.debug(missing)  ## If debug also print genomes to terminal
		logger.info("Genomes not matching any annotation {len}".format(len=len(missing)))
		return missing
