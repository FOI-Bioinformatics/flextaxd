'''
Module to read and write newick trees

'''
from flextaxd.modules.database.DatabaseConnection import ModifyFunctions
from io import StringIO
import importlib

'''Temporary fix for conda that refuses to select the correct version of ete3 during test installation.
	It fails due to faces not being available in that ete3 version on import, but it works when ete3 is
	 being installed manually using conda.
'''

import sys
import logging
logger = logging.getLogger(__name__)


__version__ = "0.0.1"
__author__ = "David Sundell"
__credits__ = ["David Sundell"]
__license__ = "GPLv3"
__maintainer__ = "FOI bioinformatics group"
__email__ = ["bioinformatics@foi.se", "david.sundell@foi.se"]
__date__ = "2020-10-15"
__status__ = "Production"
__partof__ = "FlexTaxD"

class VisualisationError(Exception):
	"""Exception raised for errors in the input."""
	def __init__(self, message):
		#self.expression = expression
		self.message = message

class NewickNode(object):
	"""The NewickNode class stores the information of a taxonomy node
			ID
			name
			children
			parent
		The purpose of this class is to allow a fast printout of a newick tree.
		All nodes in a newick tree knows it's children and parent, so by selecting a
		print option for the class objects can inheritly print all its children or the lineage (parents).

		main functions
			add_child  ## Add a newick node object reference as a child of the node
			set_print  ## Set the print type  (name, id, lineage or newick <- default)
	"""

	def __init__(self, id, name, parent=False):
		super(NewickNode, self).__init__()
		self.id         = id            	## Node id
		self.name       = name          	## Node name
		self.parent     = parent        	## Parent id False for root
		self.children   = set()         	## Set of newick children
		self.__class__.print_opt = "newick" ## The default behaviour of this class is to print out a
											## 		newick tree from the given node (as root)

	def __str__(self):
		'''The print function of NewickNode allows any node to print all its children in newick format or the lineage to root
			There is also an option to print just the name of a node
		'''
		if self.__class__.print_opt == "name":
			return "{name}".format(name=self.name)
		elif self.__class__.print_opt == "lineage":
			if self.parent:
				return "{parent};{name}".format(name=self.name, parent=self.parent)
			else:
				return "root"
		elif self.__class__.print_opt == "newick":
			if len(self.children) > 0 and self.parent:
				return "({children}){name}".format(name=self.name,children=",".join(str(child) for child in self.children))
			elif not self.parent:  #Only the root node has this property
				return "(({children}){parent})ROOT;".format(children=",".join(str(child) for child in self.children),parent=self.name)
			else:
				return "{name}".format(name=self.name)
		else:
			return "{id}: {name}; children: {nchildren} ".format(id=self.id, name=self.name, nchildren = len(self.children))

	def __repr__(self):
		return "NewickNode()"

	def add_child(self,child):
		'''Add a NewickNode object as child'''
		self.children.add(child)
		return

	def set_print(self,_type):
		'''Set the variable that controls the print style
			name
			newick - the complete subtree in newick format
			lineage
		'''
		self.__class__.print_opt = _type


class NewickTree(object):
	"""The NewickTree class parses a CanSNP or a FlexTaxD database and prints the database as a newick tree

		print styles
			name
			newick - the complete subtree in newick format
			lineage - prints all parents of a node up to root

	"""

	def __init__(self, database,name="newick",outdir="./",taxid=False):
		super(NewickTree, self).__init__()
		self.database = ModifyFunctions(database) ## Initiate database connection with CanSNPdbFunctions
		if taxid:
			self.taxid = self.database.get_id(taxid)
		self.nodeDict = {}							## Dictionary to store references to all newick nodes
		self.c_p_set = set()
		self.tree_file = "{outdir}/{name}_tree.pdf".format(outdir=outdir.rstrip("/"),name=name) ## output file
		self.tmp_tree = "{outdir}/.newick"
		## Build the newick tree
		self.newickTree = str(self.build_tree(taxid=self.taxid))

	def __repr__(self):
		return "NewickTree()"

	def print(self,type="newick"):
		'''Description function to print tree output
		Parameters
			str 		- type
			Formats
				newick 		- newick format B,(A,C,E),D);
				newick_vis	- newick format as ascii tree
				trees		- newick format as phylogenetic tree (using matplotlib)
		------
		Returns
			boolean 	- True
		'''
		if type == "newick":
			print(self.newickTree)
			return

		'''Local import allows default newickTree output to be independent of non standard python libraries'''
		exists = importlib.find_loader('Bio')
		if not exists:
			raise VisualisationError("Visualisations other than newick requires biopython package (conda install biopython)!")

		from Bio import Phylo
		self.phylo = Phylo.read(StringIO(self.newickTree), "newick")
		if type == "newick_vis":
			Phylo.draw_ascii(self.phylo)
		import matplotlib.pylab as pylab
		if type == "tree":
			Phylo.draw(self.phylo)
		#fig.set_size_inches(5,5)
		pylab.savefig("flextaxd.vis.png",)
		return True

	def get_tree(self,table="tree",taxid=False):
		'''Parameters
			table 	- table in database (default tree)
			taxid	- parent taxid
		Function that returns the whole tree in the database the script expects
			the tree to be rooted at the lowest value and that the root has itself as parent

		------
		Returns
			lists	- List of links in tree, list of nodes in tree
			or
			list	- List of links in tree, False
		'''
		if taxid:
			nodes = self.database.get_children([taxid])
			return self.database.get_links(nodes,order=True),nodes
		else:
			SELECT = "SELECT parent,child,rank_i FROM {table} ORDER BY child ASC".format(table=table)
			return self.database.query(SELECT).fetchall(),False

	def get_nodes(self, names=False,col=False):
		'''Retrieve the whole node info table of the database to decrease the number of database calls!

		------
		Returns
			dict - dictionary with node id to name and name to node id
		'''
		nodeDict = {}
		QUERY = '''SELECT id,name FROM nodes'''
		if names:
			QUERY += " WHERE id in({nodes})".format(nodes=",".join(map(str,names)))
		logger.debug(QUERY)
		for node in self.database.query(QUERY).fetchall():
			nodeDict[node[0]] = node[1]
			if col == 1:
				continue
			nodeDict[node[1]] = node[0]
		return nodeDict

	def get_parent(self,name):
		'''return parent'''
		#QUERY = '''SELECT parent,child,rank FROM tree LEFT JOIN rank on (tree.rank_i = rank.rank_i) WHERE child = "{node}"'''.format(node=name)
		QUERY = '''SELECT parent,child,rank_i FROM tree WHERE child = "{node}"'''.format(node=name)
		res = self.database.query(QUERY).fetchone()
		return res

	def get_child(self,name):
		'''return child'''
		#QUERY = '''SELECT child,child,rank FROM tree LEFT JOIN rank on (tree.rank_i = rank.rank_i) WHERE child = "{node}"'''.format(node=name)
		QUERY = '''SELECT parent,child,rank_i FROM tree WHERE parent = "{node}"'''.format(node=name)
		res = self.database.query(QUERY).fetchone()
		return res

	def new_node(self,child,nodes,parent):
		'''Function that adds a new node to the newick tree'''
		try:
			if len(self.c_p_set & set([(child,parent)])) == 0:
				node = NewickNode(child, nodes[child], self.nodeDict[parent])  		## this works also for root as root has itself as child
				'''Make sure link to parent was not made before'''
				self.c_p_set |= set([(child,parent)])
				self.c_p_set |= set([(parent,child)])
				self.nodeDict[child] = node
			else:
				logger.debug("Link {child}-{parent} already exists, retrieve node!".format(child=child,parent=parent))
				node = self.nodeDict[child]										## add a reference to the node so that children can be added
			pnode = self.nodeDict[parent]										## add child to the parent node
			pnode.add_child(node)
		except KeyError:
			logger.debug("Error in adding NewickNode parent: {parent} does not exist, trying to add parent".format(parent=parent))
			t_parent,t_child,rank = self.get_parent(parent)
			logger.debug("Adding parent: {parent} of parent: {child}".format(parent=t_parent,child=t_child))
			self.new_node(t_child,nodes,t_parent)
			self.new_node(child,nodes,parent)
			return
		logger.debug("NewickNode p:{parent} c: {child} added".format(parent=parent,child=child))
		return

	def build_tree(self,taxid=False):
		'''Build newick tree from database
			This function walks through a database of nodes and creates NewickNode objects
			self-aware of their decending newick tree or their parent lineage,

			Parameters: Select a taxid on which to start from instead of root

			Returns: The root of the tree, however all nodes are accesible from the
						NewickTree nodeDict by their node name
		'''
		tree,nodes = self.get_tree(taxid=taxid)
		nodes = self.get_nodes(nodes & set([self.taxid]))
		logger.debug("Nodes: {n} Links: {l}".format(n=len(nodes),l=len(tree)))
		logger.debug([nodes,tree])

		'''Add root node as incoming taxid'''
		if taxid:
			parent,child,rank = self.get_child(taxid)
			root = NewickNode(self.taxid, nodes[self.taxid], False)					## Create the root node
			self.nodeDict[self.taxid] = root										## Add the root node to the node dictionary
			self.nodeDict["root"] = root									## Also add this reference as "root"
			newickTree = root  												## The master parent will contain the full tree
		for parent,child,rank in tree:
			if parent == child:  ## root
				root = NewickNode(child, nodes[child], False)					## Create the root node
				self.nodeDict[child] = root										## Add the root node to the node dictionary
				self.nodeDict["root"] = root									## Also add this reference as "root"
				newickTree = root  												## The master parent will contain the full tree
				continue
			self.new_node(child,nodes,parent)

		## The newickTree is the same as the master parent (containing the full tree, the root node is defined on row 135)
		logger.debug("Tree complete, return newickTree")
		#logger.debug(newickTree)
		return newickTree
