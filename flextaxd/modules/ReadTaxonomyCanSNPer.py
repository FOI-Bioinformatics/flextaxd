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
__date__ = "2019-03-11"
__status__ = "Production"

from .ReadTaxonomy import ReadTaxonomy

class ReadTaxonomyCanSNPer(ReadTaxonomy):
	"""docstring for ReadTaxonomyCanSNPer."""
	def __init__(self, taxonomy_file=False, database=".canSNPdb",  taxid_base=1,root_name="Francisella",verbose=False):
		super(ReadTaxonomyCanSNPer, self).__init__(taxonomy_file=taxonomy_file, database=database,verbose=False)
		self.input = taxonomy_file
		self.taxonomy = {}
		self.taxid_base = taxid_base
		## Initiate database
		root_i = self.add_node(root_name)
		self.taxid_base -=1 ## reset
		self.add_rank("no rank",ncbi=True)
		self.add_link(root_i, root_i,rank="no rank")
		#self.database.commit()
		self.names = {}
		self.root = root_i
		self.length = 0
		self.ids = 0

	def add_SNP(self,nodes,i):
		'''Add name of node to database'''
		name = nodes[i].strip()
		node_i = self.add_node(name)  ## Add node
		#i-=1 ## Check next parent
		try:
			parent_i = self.taxonomy[nodes[i-1]]  ## check if current nodes parent exists
		except KeyError:  						## Parent did not exist, keep walking up the tree and add all parents until root!
			if i < -len(nodes):
				return self.root				## The root has been reached return index of root
			else:								## Parent didn't exist again, add parent to this node and then add the link to that parent
				parent_i = self.add_SNP(nodes,i-1)## Add node of parent
		self.add_link(node_i, parent_i)			## Add link to next parent
		return node_i 	##  index of child

	def parse_taxonomy(self):
		'''Retrieve node description from CanSNPer formatted tree'''
		if self.verbose: print("Parse CanSNP tree file")
		with open(self.taxonomy_file,"r") as f:
			for row in f:
				nodes = row.strip().split(";")  ## get node and all its parents in a list
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
		self.database.commit()									## Commit changes to database
		self.length = self.taxid_base - self.root				## Check number of new nodes added
		print("New taxonomy ids assigned {taxidnr}".format(taxidnr=self.length))
