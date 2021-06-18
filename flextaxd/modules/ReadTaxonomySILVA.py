#!/usr/bin/env python3 -c

'''
Read SILVA taxonomy holds a dictionary with taxonomy tree and name translation
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


class ReadTaxonomySILVA(ReadTaxonomy):
	"""docstring for ReadTaxonomySILVA."""
	def __init__(self, taxonomy_file=False, database=".silva",  taxid_base=1,root_name="root",rank="family", verbose=False,**kwargs):
		super(ReadTaxonomySILVA, self).__init__(taxonomy_file=taxonomy_file, database=database,verbose=verbose,**kwargs)
		self.input = taxonomy_file
		self.taxid_num = taxid_base
		## Initiate database
		logger.info(taxonomy_file)
		self.length = 0
		self.database.commit()
		logger.debug("Root: {root} link: [{base}]".format(root=self.root, base=self.taxid_num))

	def parse_taxonomy(self):
		'''Retrieve node description from SILVA formatted tree'''
		logger.info("Parse SILVA tree file")
		with self.zopen(self.taxonomy_file,"r") as f:
			for row in f:
				tree,info = row.strip().rsplit(";",1) ## separate tree from info columns
				nodes = tree.strip().split(";")  ## get node and all its parents in a list.
				info = info.lstrip("\t")
				try:
					taxid,rank,_ = info.split("\t",2)
				except ValueError:
					taxid,rank = info.split("\t")
				child = nodes[-1]  ## get name of child node
				self.add_rank(rank)
				'''If the tree was not properly formatted an a parent is missing make sure that function works anyway by adding any parent node above child'''
				try:
					pname = nodes[-2].strip()
					parent_i = self.taxonomy[pname]  ## Check if parent of child node exists
				except IndexError: ## Should be first row with only one node (parent)
					child_i = self.add_node(child,id=taxid)
					self.add_link(child_i,self.root,rank=rank)
					continue
				'''Now all parents exists, add new child node and add the link'''
				child_i = self.add_node(child,id=taxid)
				self.add_link(child_i,parent_i,rank=rank)
		self.database.commit()
		self.length = len(self.taxonomy)                ## Check number of new nodes added
		self.taxid_num = len(self.taxonomy)
		logger.info("New taxonomy ids assigned {taxidnr}".format(taxidnr=self.length))
