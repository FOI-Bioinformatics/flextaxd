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

try:
	import inquirer
except:
	Raise(ImportError("The package inquirer is required for flextaxd visualisation"))


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

	def __init__(self, database,name="newick",outdir="./",taxid=False,maxdepth=3,label_size=None,vis_clip_labels=False):
		super(NewickTree, self).__init__()
		self.database = ModifyFunctions(database) ## Initiate database connection with CanSNPdbFunctions
		if taxid:
			logger.info("Visualise {taxid}".format(taxid=taxid))
			self.taxid = self.database.get_id(taxid)
		self.nodeDict = {}							## Dictionary to store references to all newick nodes
		self.c_p_set = set()
		self.tree_file = "{outdir}/{name}_tree.pdf".format(outdir=outdir.rstrip("/"),name=name) ## output file
		self.tmp_tree = "{outdir}/.newick"
		self.taxonomy = self.database.get_nodes()
		## Build the newick tree
		self.maxdepth = maxdepth
		self.label_size = label_size
		self.vis_clip_labels = vis_clip_labels
		self.newickTree = str(self.build_tree(taxid=self.taxid,maxdepth=self.maxdepth))

		## Check difference between database and constructed tree
		# Get nodes from database
		db_index_nodes = self.get_nodes(col=1)
		db_nodes = set()
		for index,node in db_index_nodes.items():
			if node == 'cellular organisms': continue # do not add this arbitrary taxonomy
			db_nodes.add(node.lower())
		#/
		# Get nodes from tree
		"""
		# Ete3 returns full node names. Bio-phylo returns only the last work after splitting for spaces. But we dont use Ete3 anywhere else, so I stick with Bio-phylo for now.
		#                               This means exact-matches between nodes in database and tree cannot be made. We have to resort to comparing number of nodes.
		from ete3 import Tree
		tree = Tree(self.newickTree,format=1)
		num_nodes = 0
		tree_nodes = set()
		for node in tree.traverse():
			tree_nodes.add(node.name.lower())
			num_nodes += 1
		"""
		from Bio import Phylo
		tree = Phylo.read(StringIO(self.newickTree), "newick")
		tree_nodes = set()
		for node in tree.find_clades():
			try:
				tree_nodes.add(node.name.lower())
			except:
				logger.info('Warning: node did not have a name: '+str(node))
		#/
		#print(tree_nodes.difference(db_nodes))
		#print(db_nodes.difference(tree_nodes))
		# Check if identical
		if not len(db_nodes) == len(tree_nodes):
			print('Warning: there is a difference between the database and newick tree. The displayed tree may lack nodes. This error-message is normal if you are not visualising from root.')
			logger.debug('Inconsistent datasets returned from database and the contructed Newick tree.\n'+','.join(map(str,db_nodes))+'\n'+','.join(map(str,tree_nodes)))
		#/
		##/

	def __repr__(self):
		return "NewickTree()"

	def set_max_depth(self,depth):
		'''Change the object maxdepth'''
		self.maxdepth = depth
		return self.maxdepth

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
			raise VisualisationError("Visualisations other than newick requires biopython package (mamba install biopython)!")
		from Bio import Phylo
		self.phylo = Phylo.read(StringIO(self.newickTree), "newick")
		if type == "newick_vis":
			Phylo.draw_ascii(self.phylo)
		exists = importlib.find_loader('matplotlib')
		if not exists:
			raise VisualisationError("Visualisations using newick tree requires the matplotlib package (mamba install matplotlib-base)!")
		import matplotlib.pylab as pylab
		import matplotlib
		matplotlib.rcParams['axes.spines.top'] = False
		matplotlib.rcParams['axes.spines.right'] = False
		if self.label_size:
			matplotlib.rc('font', size=self.label_size)
		if type == "tree":
			if self.vis_clip_labels:
				Phylo.draw(self.phylo)
			else:
				Phylo.draw(self.phylo,label_func=lambda leaf: leaf.name) # Replace Biopython label function that returns the full name (default, clip >40)
		#fig.set_size_inches(5,5)
		pylab.savefig("flextaxd.vis.png",)
		return True

	def double_opts_vis(self,links,taxid="name",full=False):
		'''vis double opt'''
		import inquirer
		tn = self.database.get_nodes()
		if taxid != "name":
			taxid = tn[taxid]
		parents = [tn[x[0]] for x in links]
		questions = [
		  inquirer.List('parent',
		                message="The node {name} has multiple optional parents, select which line to visualise: ".format(name=taxid),
		                choices=parents,
		            ),
		]
		answers = inquirer.prompt(questions)
		selected = tn[answers["parent"]] ## Return back to id
		taxids = []
		unique = set()
		for link in links:
			if link[0] == selected:
				rank = link[2]
				taxids.append(link[1])
				'''Update incoming taxid which may not be correct'''
				self.taxid = link[1]
			unique.add(link[1])
		rank = rank +1
		if len(unique) > len(taxids):
			rank = False
		return rank,taxids ## Parent rank is selected, child rank is plus one

	def get_tree(self,table="tree",taxid=False,maxdepth=3):
		'''Parameters
			table 	- table in database (default tree)
			taxid	- parent taxid
			depth	- depth
		Function that returns the whole tree in the database the script expects
			the tree to be rooted at the lowest value and that the root has itself as parent

		------
		Returns
			lists	- List of links in tree, list of nodes in tree
			or
			list	- List of links in tree, False
		'''
		selected = False
		if taxid:
			if maxdepth == 0:
				maxdepth = 1000  ## It is not reasonable to expect trees with more than 1000 levels, if so bug has to be raised
			'''Check if double parent'''
			links = self.database.get_links([taxid],order=True)
			nodes = self.database.get_node(self.taxonomy[taxid])
			find_tax = [taxid]
			if len(nodes) > len(links):
				links = self.database.get_links(nodes,order=True)
			if len(links) > 1:
				selected,find_tax = self.double_opts_vis(links,taxid)
			nodes = self.database.get_children(find_tax,maxdepth=maxdepth,selected=selected)
			if len(nodes) == 0:
				raise VisualisationError("Given node has no children")
			links = self.database.get_links(nodes,order=True)
			return links,nodes
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
		logger.debug(QUERY)
		res = self.database.query(QUERY).fetchone()
		return res

	def get_child(self,name,rank_i=False):
		'''return child'''
		#QUERY = '''SELECT child,child,rank FROM tree LEFT JOIN rank on (tree.rank_i = rank.rank_i) WHERE child = "{node}"'''.format(node=name)
		QUERY = '''SELECT parent,child,rank_i FROM tree WHERE parent = "{node}"'''.format(node=name)
		if rank_i:
			QUERY = '''SELECT parent,child,rank_i FROM tree WHERE parent = "{node}" and rank_i ={rank}'''.format(node=name,rank=rank_i)
		logger.debug(QUERY)
		res = self.database.query(QUERY).fetchone()
		return res

	def new_node(self,child,nodes,parent,link=1):
		'''Function that adds a new node to the newick tree'''
		if self.link_exists == (child,parent,link): ## This has already been tried return directly
			return False
		try:
			if len(self.c_p_set & set([(child,parent,link),(parent,child,link)])) == 0: ## If link does not exist
				try:
					node = NewickNode(child, nodes[child], self.nodeDict[parent])  		## this works also for root as root has itself as child
				except KeyError:
					logger.debug("Unable to add NewickNode for child={child} and parent={parent}. Returning -1 as flag to add later".format(child=child,parent=parent))
					return -1 # all returns with -1 will be re-added in a secondary loop
				'''Make sure link to parent was not made before'''
				self.c_p_set |= set([(child,parent,link)])
				self.c_p_set |= set([(parent,child,link)])
				self.nodeDict[child] = node
				self.link_exists=False
			else:
				logger.debug("Link {parent}-{child} already exists, retrieve node!".format(child=child,parent=parent))
				self.link_exists = (child,parent,link)
				node = self.nodeDict[child]										## add a reference to the node so that children can be added
			pnode = self.nodeDict[parent]										## add child to the parent node
			pnode.add_child(node)
		except KeyError:
			logger.debug("Error in adding NewickNode of {child} parent: {parent} does not exist, trying to add parent".format(parent=parent,child=child))
			try:
				t_parent,t_child,rank = self.get_parent(parent)
			except KeyError:
				raise Error("Something is wrong")
			logger.debug("Adding parent: [{parent},{pname}] of child: [{child},{cname}]".format(parent=t_parent,pname=self.taxonomy[t_parent],cname=self.taxonomy[t_child],child=t_child))
			if self.new_node(t_child,nodes,t_parent,link=link-1): ## Add the missing parent
				self.new_node(child,nodes,parent,link)	  ## Try again to add the node
			else:
				self.added_parent = child
				return False
			return True
		logger.debug("NewickNode p:{parent} c: {child} added".format(parent=parent,child=child))
		return True

	def unique_indexes(self,nodes):
		'''Check duplicated indexes and give them unique IDs before print'''
		QUERY = "SELECT child FROM tree WHERE child in ({nodes}) GROUP BY child HAVING count(parent) > 1 ".format(nodes=",".join(map(str,nodes)))  ## Thanks to andrewjmc@github for this suggestion
		child_w_dpi = self.database.query(QUERY).fetchall()  ## Fetch all conflicting links and give them unique index before printing
		child_w_dpi = list(*child_w_dpi)
		return child_w_dpi

	def fix_names(self,nodes,tn,tr):
		'''Return names instead of taxid'''
		nnodes = []
		for x in range(len(nodes)):
			fixnames = list(nodes[x])
			fixnames[0] = tn[fixnames[0]]
			fixnames[1] = tn[fixnames[1]]
			fixnames[2] = tr[fixnames[2]]
			nnodes.append(fixnames)
		return nnodes


	def double_vis_path(self,duplicates,taxid,nodes):
		'''The range in the tree has two identical nodes, ask user to resolve which path goes where'''
		import inquirer,random
		tn = self.database.get_nodes()
		tr = self.database.get_rank()
		children = self.database.get_links(self.database.get_children([taxid],maxdepth=1))
		parents = duplicates[taxid]
		if len(children) > 2:
			raise VisualisationError("Names occuring three times in the same tree are not taken care of at this time, export function does however!")
		children_N = self.fix_names(children,tn,tr)
		parents_N = self.fix_names(parents,tn,tr)
		selto = parents_N[0]
		default = children_N[0]
		if taxid != "name":
			taxid = tn[taxid]
		questions = [
		  inquirer.List('parent',
		                message="The node {name} has two paths, select the correct path to ({p1})".format(name=taxid,p1=selto,children=children_N),
		                choices=children_N,
		            ),
		]
		selected = inquirer.prompt(questions)["parent"]
		n = random.randint(10000000,10100000)
		if selected != default:
			children[0],children[1] = children[1],children[0]
		## Change id for first pair
		pc = list(parents[0])
		pc[1] = n
		parents[0] = tuple(pc)
		pc = list(children[0])
		pc[0] = n
		children[0] = tuple(pc)
		nodes[n] = taxid
		return [*children,*parents],nodes

	def duplicate_parents_check(self,nodes,tree,taxid):
		'''Check if any node have optional parents
			Parameters:
				nodes - the list of nodes of interest
			Returns:
				list - list of nodes with optional parents
		'''
		duplicates = self.unique_indexes(nodes)
		optional_parents = {}
		for child in duplicates:
			optional_parents[child] = self.database.get_parent(child,all=True)
		return optional_parents

	def fix_tree(self,tree,duplicates,nodes):
		'''Change index for one all but one duplicate to make sure tree is printed correctly'''
		update_tree = []
		dup = [duplicates[k] for k in duplicates.keys()]
		for link in tree:
			if len(set([link]) & set(*dup)) > 0:
				pass
			else:
				update_tree.append(link)
		changes = []
		for child in duplicates.keys():
			changes,nodes = self.double_vis_path(duplicates,child,nodes)
			update_tree += changes
		tree = sorted(update_tree,key=lambda x:x[2])
		return tree,nodes

	def build_tree(self,taxid=False,maxdepth=3,check_parent=False):
		'''Build newick tree from database
			This function walks through a database of nodes and creates NewickNode objects
			self-aware of their decending newick tree or their parent lineage,

			Parameters:
				taxid - Select a taxid on which to start from instead of root
				depth - the number of levels downstream to visualise, default(5) to avoid too large trees for visualisation
			Returns: The root of the tree, however all nodes are accesible from the
						NewickTree nodeDict by their node name
		'''
		tree,nodes = self.get_tree(taxid=taxid,maxdepth=maxdepth)
		duplicates = self.duplicate_parents_check(nodes - set([self.taxid]),tree,taxid)
		nodes = self.get_nodes(nodes | set([self.taxid]),col=1)
		if len(duplicates) > 0:
			tree,nodes = self.fix_tree(tree,duplicates,nodes)
		logger.debug("Nodes: {n} Links: {l}".format(n=len(nodes),l=len(tree)))
		logger.debug([nodes,tree])
		self.added_parent = False
		self.link_exists=False
		'''Add root node as incoming taxid'''
		if taxid:
			parent,child,rank = self.get_child(taxid)
			root = NewickNode(self.taxid, nodes[self.taxid], False)					## Create the root node
			self.nodeDict[self.taxid] = root										## Add the root node to the node dictionary
			self.nodeDict["root"] = root									## Also add this reference as "root"
			newickTree = root  												## The master parent will contain the full tree

		post_add = [] # keep track of which entries in tree failed due to error when adding
		for parent,child,rank in tree:
			if parent == child:  ## root
				root = NewickNode(child, nodes[child], False)					## Create the root node
				self.nodeDict[child] = root										## Add the root node to the node dictionary
				self.nodeDict["root"] = root									## Also add this reference as "root"
				newickTree = root  												## The master parent will contain the full tree
				continue
			new_node_status = self.new_node(child,nodes,parent,link=rank)
			if new_node_status == -1:
				post_add.append([parent,child,rank])

		# Try to add all failed nodes
		if post_add:
			logger.debug('Had '+str(len(post_add))+' nodes that failed to add. Will try to rescue them now...')
		post_adds_done = set()
		maxiterations = (len(post_add)*len(post_add))+len(post_add) # set a maximum number of attempts to add (number of post-adds squared plus a margin of the number of post-adds)
		iterations_done = 0
		while not len(post_adds_done) == len(post_add):
			for enum,(parent,child,rank) in enumerate(post_add):
				new_node_status = self.new_node(child,nodes,parent,link=rank) # try to add node
				post_adds_done.add(enum) # mark done as added
				if new_node_status == -1: # if it fails...
					post_adds_done.remove(enum) #           ... then remove mark as done
			
			iterations_done += 1
			if iterations_done > maxiterations: # check if we have performed the maximum number of add attempts
				print('Warning: Many nodes returned an error when building Newick-tree. Reached maximum number of recues, aborting...')
				logger.debug('Reached limit when adding nodes, aborted. Limit='+str(maxiterations))
				break
		#/

		## The newickTree is the same as the master parent (containing the full tree, the root node is defined on row 135)
		logger.debug("Tree complete, return newickTree")
		#logger.debug(newickTree)
		return newickTree
