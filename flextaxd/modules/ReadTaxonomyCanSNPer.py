#!/usr/bin/env python3 -c

'''
Read CanSNPer taxonomy holds a dictionary with taxonomy tree and name translation
'''

__version__ = "1.0"
__author__ = "David Sundell"
__credits__ = ["David Sundell"]
__license__ = "GPLv3"
__maintainer__ = "FOI bioinformatics group"
__email__ = "bioinformatics@foi.se"
__date__ = "2020-01-17"
__status__ = "Production"

from .ReadTaxonomy import ReadTaxonomy
import logging
logger = logging.getLogger(__name__)

class ImportFormatError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)


class ReadTaxonomyCanSNPer(ReadTaxonomy):
	"""docstring for ReadTaxonomyCanSNPer."""
	def __init__(self, taxonomy_file=False, database=".canSNPdb",  taxid_base=1,root_name=False,rank="family", verbose=False,**kwargs):
		super(ReadTaxonomyCanSNPer, self).__init__(taxonomy_file=taxonomy_file, database=database,verbose=verbose,**kwargs)
		self.input = taxonomy_file
		self.taxonomy = {}
		self.taxid_num = taxid_base
		## Initiate database
		logger.info(taxonomy_file)
		logger.debug(root_name)
		if not root_name:
			logger.info("Fetching root name from file")
			root_name = self.get_root_name(taxonomy_file)
		logger.info("Adding, cellular organism node")
		c_i = self.add_node("cellular organisms")
		logger.info("Adding root node {node}!".format(node=root_name))
		root_i = self.add_node(root_name)
		self.taxid_num =taxid_base ## reset
		logger.debug("Adding ranks!")
		if rank != "no rank":
			self.add_rank("no rank")
		self.add_rank(rank)
		if root_name != "cellular organisms":
			self.add_link(c_i,1)
		if root_name != "root":
			self.add_link(root_i,c_i,rank=rank) ## root is always added as top parent if not one force 1
			self.taxonomy[root_i] = 1  ## Add local link to root
		self.names = {}
		self.root = root_i
		self.length = 0
		self.ids = 0

		self.database.commit()
		logger.debug("Root: {root} link: [{base}]".format(root=self.root, base=self.taxid_num))

	def get_root_name(self,fin):
		'''Get the root name of the CanSNP file unless root is specified'''
		with self.zopen(fin) as f:
			firstline = f.readline()
			firstline = firstline.strip().replace("\t",";")
			root = firstline.split(";") ## Root should be either the leftmost annotation in the tree or the only annotation on that row
			if isinstance(root,list):
				root = root[0]
		return root

	def add_SNP(self,nodes,i):
		'''Add name of node to database'''
		name = nodes[i].strip()
		logger.debug("Parent did not exit Add parent: {name}".format(name=name))

		#i-=1 ## Check next parent
		try:
			parent_i = self.taxonomy[nodes[i-1]]  ## check if current nodes parent exists
		except KeyError:                          ## Parent did not exist, keep walking up the tree and add all parents until root!
			if i < -len(nodes):
				return self.root                ## The root has been reached return index of root
			else:                                ## Parent didn't exist again, add parent to this node and then add the link to that parent
				parent_i = self.add_SNP(nodes,i-1)## Add node of parent
		node_i = self.add_node(name)  ## Add node
		if node_i:
			self.taxid_num += 1
		self.add_link(node_i, parent_i)            ## Add link to next parent
		return node_i     ##  index of child

	def parse_taxonomy(self):
		'''Retrieve node description from CanSNPer formatted tree'''
		logger.info("Parse CanSNP tree file")
		with self.zopen(self.taxonomy_file,"r") as f:
			for row in f:
				row = row.strip().replace("\t",";")  ## Also accept tab separated tree files
				nodes = row.strip().split(";")  ## get node and all its parents in a list
				logger.debug(nodes)
				child = nodes[-1]  ## get name of child node
				'''If the tree was not properly formatted an a parent is missing make sure that function works anyway by adding any parent node above child'''
				try:
					pname = nodes[-2].strip()
					parent_i = self.taxonomy[pname]  ## Check if parent of child node exists
				except KeyError: ## parent node does not exist, add parent
					parent_i = self.add_SNP(nodes,-2)
				except IndexError: ## Should be first row with only one node (parent)
					self.taxonomy[child] = self.root
					parent_i = self.root
					continue
				'''Now all parents exists, add new child node and add the link'''
				child_i = self.add_node(child)
				self.database.add_link(child_i,parent_i)
		self.database.commit()
		self.length = self.taxid_num - self.root                ## Check number of new nodes added
		logger.info("New taxonomy ids assigned {taxidnr}".format(taxidnr=self.length))
