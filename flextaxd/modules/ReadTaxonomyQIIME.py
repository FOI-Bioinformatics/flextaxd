#!/usr/bin/env python3 -c

'''
Read QIIME formatted taxonomy files and holds a dictionary with taxonomy tree and name translation
'''

from .ReadTaxonomy import ReadTaxonomy

class ReadTaxonomyQIIME(ReadTaxonomy):
	"""docstring for ReadTaxonomyQIIME."""
	def __init__(self, taxonomy_file=False, names_dmp=False, database=False,  taxid_base=1):
		super(ReadTaxonomyQIIME, self).__init__(database=database)
		self.input = taxonomy_file
		self.taxonomy = {}
		self.names = {}
		self.taxid_base = taxid_base
		self.taxonomy = {}
		self.length = 0
		self.ids = 0
		self.levelDict = {
				"n": "no rank",
				"sk": "superkingdom",
				"k": "kingdom",
				"d": "domain",
				"p": "phylum",
				"c": "class",
				"o": "order",
				"f": "family",
				"g": "genus",
				"s": "species"
		}

		### Add root name these manual nodes are required when parsing the GTDB database!
		self.add_node("root")
		self.add_node("cellular organism")
		self.add_node("Bacteria")

		self.add_rank("n")
		self.add_rank("sk")
		## Add basic links
		self.add_link(1, 1,rank="n")
		self.add_link(2, 1,rank="n")
		self.add_link(3, 2,rank="sk")

	def parse_taxonomy(self):
		'''Parse taxonomy information'''
		self.qiime_to_tree()

	def parse_tree(self,tree,current_i=0):
		'''The taxonomy tree does not exist in the standard nomenclature, add a new tree'''
		### Add parent description if not exist
		level,description = self.parse_description(tree,current_i)
		if description.strip() == "":
			return False
		if level not in self.taxonomy:
			self.add_rank(level)
		if current_i == len(tree)-1:
			'''Top parent reached return top parent id'''
			try:
				return self.taxonomy[description]
			except KeyError:
				node_i = self.add_node(description)
			return node_i
		parent_i = self.parse_tree(tree,current_i+1)
		try:
			'''If parent is already in the database return parent ID'''
			return self.taxonomy[description]
		except KeyError:
			''' Add current node to names file '''
			node_i = self.add_node(description)

		'''When all parents exist add current relation to tree file'''
		self.add_link(node_i,parent_i,rank=level)
		'''return new tax_i'''
		return node_i

	def parse_description(self,tree,current_i):
		'''Retrieve node description from QIIME formatted tree'''
		current_level = tree[current_i]
		level, description = current_level.split("__")
		return level,description

	def qiime_to_tree(self, sep="\t"):
		'''Read the qiime format file and parse out the relation tree (nodes.dmp)'''
		self.sep = sep
		self.tree = set()
		taxid_start = self.taxid_base
		with open(self.input) as f:
			'''Each row defines a GCF genome file connected to a tree level'''
			for row in f:
				if row.strip() != "":  ## If there are trailing empty lines in the file
					data = row.strip().split("\t")
					try:
						genome_id = data[0].split("_",1)[1]   ## Genome ID
					except IndexError:
						pass
					### Walk through tree and make sure all nodes back to root are annotated!
					taxonomy = list(reversed(data[-1].split(";")))
					taxonomy_i = self.parse_tree(taxonomy)
					if taxonomy_i:
						self.database.add_genome(id=taxonomy_i,genome=genome_id)
		self.database.commit()
		self.length = self.taxid_base - taxid_start
		print("New taxonomy ids assigned {taxidnr}".format(taxidnr=self.length))
