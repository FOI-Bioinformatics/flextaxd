#!/usr/bin/env python3 -c

'''
Read NCBI taxonomy dmp files (nodes or names) and holds a dictionary
'''

from .ReadTaxonomy import ReadTaxonomy
from gzip import open as zopen
import os

class ReadTaxonomyNCBI(ReadTaxonomy):
	"""docstring for ReadTaxonomyNCBI."""
	def __init__(self, taxonomy_file=False, database=False):
		super(ReadTaxonomyNCBI, self).__init__(database=database)
		self.taxonomy_file = taxonomy_file
		self.names_dmp = taxonomy_file.replace("nodes","names")
		self.names = {}
		self.length = 0
		self.ids = 0
		self.accessionfile = False
		self.add_rank("no rank",ncbi=True)

	def set_accession_file(self,file):
		self.accessionfile = file

	def parse_taxonomy(self,treefile=False):
		'''Parse taxonomy information'''
		if hasattr(self, "names_dmp"):
			print("Parse names file {}".format(self.names_dmp))
			self.read_names(self.names_dmp)
			self.length = self.database.num_rows("nodes")
		if self.taxonomy_file:
			print("Parse nodes file {}".format(self.taxonomy_file))
			self.read_nodes(self.taxonomy_file)
			self.ids = self.database.num_rows("tree")

	def read_nodes(self, taxfile):
		'''Read NCBI node file and store a dictionary with names in object'''
		with open(taxfile, "r") as _taxfile:
			for taxonomy_row in _taxfile:
				data = taxonomy_row.strip().split("\t|\t")
				child,parent,rank = data[0],data[1],data[2]
				if rank not in self.rank:
					self.add_rank(rank,ncbi=True)
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
					self.database.add_node(name, id=taxid)
			self.database.commit()
		return

	def parse_genebank_file(self,filepath,filename):
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
		if self.verbose: print("Parsing ncbi accession2taxid")
		self.refseqid_to_GCF = {}
		for root, dirs, files in os.walk(genomes_path):
			for filename in files:
				if filename.endswith(".fna.gz"):
					filepath = os.path.join(root, filename)
					self.parse_genebank_file(filepath,filename)
		if self.verbose: print("genomes folder read {n} sequence files found".format(n=len(self.refseqid_to_GCF)))
		with zopen(annotation_file,"r") as f:
			headers = f.readline().split(b"\t")
			for row in f:
				if row.strip() != "": ## If there are trailing empty lines in the file
					refseqid,taxid = row.split(b"\t")[1:3]
					try:
						genebankid = self.refseqid_to_GCF[refseqid]
						self.database.add_genome(id=taxid.decode("utf-8"),genome=genebankid)
					except:
						pass
			self.database.commit()
		return
