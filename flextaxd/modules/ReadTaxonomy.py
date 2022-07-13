#!/usr/bin/env python3 -c

'''
Read SILVA taxonomy dmp files (nodes or names) and holds a dictionary
'''

from .database.DatabaseConnection import DatabaseFunctions
import logging
import gzip
import sys
logger = logging.getLogger(__name__)

class InputError(Exception):
	"""Exception raised for errors in the input."""
	def __init__(self, message):
		#self.expression = expression
		self.message = message

class ReadTaxonomy(object):
	"""docstring for ReadTaxonomy."""
	def __init__(self, taxonomy_file=False, taxonomy_name=False, database=False,verbose=False,**kwargs):
		super(ReadTaxonomy, self).__init__()
		### Connect to or create database
		if database:
			self.database = DatabaseFunctions(database,verbose=verbose)
		else:
			raise InputError("No database was provided to ReadTaxonomy, abort!")
		self.taxonomy_file = taxonomy_file
		self.verbose = verbose
		self.taxonomy = {}
		self.rank = {}
		self.levelDict = {}
		self.sep = "\t"
		self.length = 0
		self.ids = 0
		self.qiime = 0

		## Add base node
		self.root = self.add_node("root")
		self.add_rank("no rank")
		## Add basic link
		self.add_link(child=1, parent=1,rank="no rank")

	def set_qiime(self,opt):
		self.qiime = opt

	def zopen(self, path,*args, **kwargs):
		'''Redefine open to handle zipped files automatically'''
		if path.endswith(".gz"):
			if str(*args) in ["r","w"] and sys.version_info.major >= 3:
				## Python version three and above interprets binary formats as binary, add t (rt) to get text returned
				args = (str(*args)+"t",)
			else:
				args = ("rt")
			return gzip.open(path,*args,**kwargs)
		return open(path,*args,**kwargs)

	#@staticmethod
	def parse_taxonomy(self,treefile=False):
		'''Parse taxonomy information'''
		self.read_nodes(treefile=treefile)

	def set_separator(self,sep):
		'''Change separator from default (\t)'''
		self.sep = sep

	def add_rank(self, rank,qiime=False):
		'''Insert rank into database'''
		try:
			return self.rank[rank]
		except KeyError:
			pass
		if rank:
			if qiime:
				rank_i = self.database.add_rank(self.levelDict.get(rank))
			else:
				rank_i = self.database.add_rank(rank)
			self.rank[rank] = rank_i
		else:
			return False
		return rank_i

	def add_link(self, child=None, parent=None,rank="no rank"):
		'''Add relationship in tree'''
		#print(child,parent,rank)
		if not child and parent:
			raise InputError("A link requires both a child and a parent!")
		self.database.add_link(child=child,parent=parent,rank=self.rank[rank])
		self.ids+=1

	def add_node(self, description,id=False):
		'''Add node to tree
			extend databaseFunction add_node function using self.taxonomy
			dict to keep track of which nodes were already added
		'''
		logger.debug("Node desc: {desc} id: {id}".format(desc=description,id=id))
		self.taxid_base = self.database.add_node(description,id)
		self.taxonomy[description] = self.taxid_base
		return self.taxid_base


	def read_nodes(self, treefile=False):
		'''Read a tab separated node file and store a dictionary with names in object'''
		logger.info("Read nodes in taxonomy file {}".format(treefile))
		swap = False
		rank = "no rank"  #Base rank if rank is not used
		if not treefile:
			treefile = self.taxonomy_file
		with self.zopen(treefile, "r") as _treefile:
			headers = _treefile.readline().strip().split(self.sep)
			if "parent" not in headers or "child" not in headers:
				logger.debug("Headers:  {h} separator [{sep}]".format(h=headers,sep=self.sep))
				raise InputError("Your input tree file does not contain the headers to specify child and parent!")
			if headers[0] == "parent":
				swap = True
			logger.debug("Swap: {swap}".format(swap=swap))
			for tree_row in _treefile:
				data = tree_row.strip().split(self.sep)
				logger.debug(data)
				if len(data) > 2:
					rank = data.pop().strip()
					if rank.strip() != "":
						self.add_rank(rank)
				if swap:
					data[0],data[1] = data[1],data[0]
				'''Add nodes into database'''
				if data[0] == "":
					'''Check for empty rows'''
					continue
				for node in data:
					logger.debug(node)
					node = node.strip()
					try:
						self.taxonomy[node]
					except KeyError:
						logger.debug("Add node: {node}".format(node=node))
						self.taxonomy[node] = self.add_node(node)
						self.ids +=1
				#logger.debug("Add link: {parent}-{child}".format(parent=data[1].strip(),child=data[0].strip()))
				self.add_link(child=self.taxonomy[data[0].strip()],parent=self.taxonomy[data[1].strip()],rank=rank)
				self.length +=1
			self.database.commit()

	def parse_genomeid2taxid(self,genomeid2taxid,reference=False):
		'''Parse file that annotates genome_id´s to nodes in the tree'''
		nodeDict = self.database.get_nodes()
		_ref = reference
		with self.zopen(genomeid2taxid,"rt") as f:
			headers = f.readline().strip().split("\t")
			for row in f:
				if row.strip() != "": ## If there are trailing empty lines in the file
					try:
						genomeid,taxid = row.strip().split("\t")
					except:
						if not _ref: ## override if there is a reference in file, and use given ref
							genomeid,taxid,reference = row.strip().split("\t")
						else:
							genomeid,taxid,override = row.strip().split("\t")
					try:
						self.database.add_genome(genome=genomeid.strip(),_id=nodeDict[taxid.strip()],reference=reference)
					except KeyError:
						logger.warning("# WARNING: {taxid} not found in the database".format(taxid=taxid))
			self.database.commit()
		return
