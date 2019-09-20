#!/usr/bin/env python3 -c

'''
Modify tree from database or source
'''

from .database.DatabaseConnection import ModifyFunctions
#from .database.database import database

class InputError(Exception):
	"""Exception raised for errors in the input."""

	def __init__(self, message):
		#self.expression = expression
		self.message = message

class ModifyTree(object):
	"""docstring for ModifyTree."""
	def __init__(self, database=".taxonomydb", mod_database=False, mod_file=False, separator="\t",verbose=False,parent=False,replace=False):
		super(ModifyTree, self).__init__()
		self.verbose = verbose
		if self.verbose: print("Modify Tree")
		self.sep = separator
		self.taxonomy = {}
		self.names = {}

		self.parent=parent
		self.replace = replace
		self.ncbi_order = True
		self.mod_genomes =False

		### Connect to or create database
		self.taxonomydb = ModifyFunctions(database,verbose=verbose)
		self.rank= self.taxonomydb.get_rank(col=2)
		## Save all nodes in the current database
		self.nodeDict = self.taxonomydb.get_nodes()

		self.taxid_base = self.taxonomydb.get_taxid_base()
		if mod_database:
			self.modsource = self.parse_mod_database(ModifyFunctions(mod_database,verbose=verbose))
		elif mod_file:
			self.modsource = self.parse_mod_file(mod_file)
		else:
			raise InputError("No modification source could be found both mod_database and mod_file file are empty!")

	def add_link(self, child, parent, rank="n"):
		'''Add relationship in tree'''
		self.taxonomydb.add_link(child,parent,self.rank[rank])
		self.ids+=1
		return

	def add_node(self, description):
		'''Add node to tree
			extend databaseFunction add_node function using self.taxonomy
			dict to keep track of which nodes were already added
		'''
		self.taxid_base = self.taxonomydb.add_node(description)
		self.taxonomy[description] = self.taxid_base
		return self.taxid_base

	def add_rank(self, rank):
		'''Insert rank into database and store into rank dictionary'''
		rank_i = self.taxonomydb.add_rank(rank)
		self.rank[rank] = rank_i
		self.rank[rank_i] = rank
		return rank_i

	def get_id(self, desc):
		'''Check if node exists if so return id, if not create a new node and return id'''
		try:
			## Check if parent node exists
			i = self.nodeDict[desc]
			return i
		except KeyError:
			'''Node does not exist, add node to the database'''
			i = self.add_node(desc)
			self.nodeDict[desc] = i
			return i

	def _parse_new_links(self, parent,child,rank):
		'''Help function for parse_mod_file, gets existing node id or adds new node'''
		## add new nodes
		parent_i = self.get_id(parent)
		self.new_nodes.add(parent_i)
		child_i = self.get_id(child)
		self.new_nodes.add(child_i)
		## add link
		if self.verbose:
			print(parent,child,rank)
			print(parent_i,child_i,rank)
		self.new_links.add((parent_i,child_i,rank))

	def parse_mod_database(self, database):
		'''Retrieve all links to update from an existing database'''
		## Retrieve annotations of existing database
		modID = database.get_id(self.parent)
		mod_links = database.get_links( database.get_children(set([modID])), database)
		## Retrieve all nodes annotated in the modified database
		self.dbmod_annotation = database.get_nodes(col=1)
		self.new_links = set()
		self.new_nodes = set()
		### Get links from modified database
		mod_links = database.get_links(self.dbmod_annotation.keys())
		self.parent_link = self.taxonomydb.get_parent(self.taxonomydb.get_id(self.parent))

		self.old_nodes = self.taxonomydb.get_children(set([self.taxonomydb.get_id(self.parent)])) - set([self.taxonomydb.get_id(self.parent)])
		### Translate node ids between databases and add non existing nodes into current database
		for link in mod_links:
			link = list(link)
			if link[0] == link[1]: ### This should only occur if the root node of the mod database is used link replace with
				self.new_links.add(self.parent_link)
			else:
				parent,child,rank = self.dbmod_annotation[link[0]].strip(),self.dbmod_annotation[link[1]].strip(),link[2]
				#print([parent,child,rank])
				self._parse_new_links(parent,child,rank)
		## remove any overlapping nodes, they don't need new indexnumber but can keep the old ones
		self.old_nodes -= self.new_nodes
		### get links from current database
		self.existing_links = set(self.taxonomydb.get_links(self.old_nodes))
		### Check which links are overlapping and may need to be updated
		self.modified_links = self.existing_links - self.new_links
		'''Get all genomes annotated to new nodes in existing database'''
		self.mod_genomes = database.get_genomes()
		return self.modified_links

	def parse_mod_file(self, modfile):
		'''Parse mod_file information'''
		self.new_links = set()
		self.new_nodes = set()
		parent = self.taxonomydb.get_id(self.parent)
		self.parent_link = self.taxonomydb.get_parent(self.taxonomydb.get_id(self.parent))
		self.old_nodes = self.taxonomydb.get_children(set([parent])) - set([parent])
		with open(modfile, "r") as f:
			headers = f.readline().strip().split(self.sep)
			if len(headers) < 2 or len(headers)>3:
				raise InputError("The modification file must contain two or three columns defined by headers parent, child (and level) separated by separator '{sep}' (default \\t)".format(sep=self.sep))
			if headers[0].lower() == "child":  ## ncbi have child on the left, swap ids to follow these rules!
				swap = True
			elif headers[0].lower() == "parent":
				swap = False
			else:
				raise InputError("The modification file does not contain proper headers, must contain child and parent")
			for row in f:
				if row.strip() == "":
					continue
				line = row.strip().split(self.sep)
				if len(line) == 2:
					line +=["-"] ## Add no specified rank to line
				parent,child,rank = line
				try:
					rank = self.rank[rank]
				except:
					rank = self.add_rank(rank)
				if swap:
					## Make sure that the order between parent and child nodes are correct and consistent with the database
					child,parent = parent,child
				self._parse_new_links(parent,child,rank)
		## remove any overlapping nodes, they don't need new indexnumber but can keep the old ones
		self.old_nodes -= self.new_nodes
		### get links from current database
		self.existing_links = set(self.taxonomydb.get_links(self.old_nodes,swap=swap))
		### Check which links are overlapping and may need to be updated
		self.modified_links = self.existing_links & self.new_links
		return self.modified_links

	def update_annotations(self, genomeid2taxid):
		'''Function that adds annotation of genome ids to nodes'''
		if self.verbose: print("Update genome to taxid annotations using {genomeid2taxid}".format(genomeid2taxid=genomeid2taxid))
		update = {
			"set_column": "id",
			"where_column": "genome",
			"set_value": "",
			"where": ""
		}
		updated = 0
		added = 0
		with open(genomeid2taxid) as f:
			for row in f:
				genome,name = row.strip().split(self.sep)
				if self.verbose: print(genome,name)
				try:
					id = self.nodeDict[name.strip()]
				except KeyError:
					if self.verbose: print("# WARNING: there was no database entry for {name} annotation not updated for this entry!".format(name=name))
				else:
					## If no exception occured add genome
					update["set_value"] = id
					update["where"] = genome.strip()
					res = self.taxonomydb.update_genome(update)
					if self.taxonomydb.rowcount()!=0:
						if res:
							updated += 1
						else:
							added += 1
					#self.taxonomydb.commit()
		self.taxonomydb.commit()
		gid = self.taxonomydb.get_genomes()
		if self.verbose: print("{added} added and {updated} genome annotations were updated!".format(added=added, updated=updated))
		return

	def update_genomes(self):
		'''When a database is supplied as source for the update genome annotations from that database needs to be transfered to the taxonomydb'''
		update = {
			"set_column": "id",
			"where_column": "genome",
			"set_value": "",
			"where": ""
		}
		updated = 0
		added = 0
		'''Database has been updated, so the internal nodeDict needs to be updated'''
		self.nodeDict = self.taxonomydb.get_nodes()
		for genome in self.mod_genomes:
			'''Translate incoming database node id to taxonomydb node id'''
			id,genomeid = self.nodeDict[self.dbmod_annotation[self.mod_genomes[genome]].strip()],genome.strip()
			update["set_value"] = id
			update["where"] = genomeid
			res = self.taxonomydb.update_genome(update)
			if self.taxonomydb.rowcount()!=0:
				if res:
					updated += 1
				else:
					added += 1
		self.taxonomydb.commit()
		if self.verbose: print("{added} added and {updated} genome annotations were updated!".format(added=added, updated=updated))
		return

	def update_database(self):
		'''Update the database file'''
		if self.replace:
			self.taxonomydb.delete_links(self.modified_links)
			self.taxonomydb.delete_links(self.existing_links-set(self.parent_link))  ## delete existing links from old nodes, except parent
			self.taxonomydb.delete_nodes(self.old_nodes)
			if self.verbose: print("Deleted nodes {nodes}".format(nodes=self.old_nodes))
		added,nodes = self.taxonomydb.add_links(self.new_links)
		if added + len(nodes) + len(self.modified_links) > 0:
			if self.replace:
				print("Deleting {n} links and {n2} nodes that are no longer valid".format(n=len(self.modified_links & self.existing_links),n2=len(self.old_nodes)))
			print("Adding {n} new nodes".format(n=len(nodes)))
			print("Adding {n} updated and/or new links".format(n=added))

			''' Commit changes (only commit once both deletion and addition of new nodes and links are completed!)'''
			self.taxonomydb.commit()
			self.nodeDict = self.taxonomydb.get_nodes()
			if self.mod_genomes:
				print("Transfering genomeid2taxid annotation from incoming database")
				self.update_genomes()
		else:
			print("All updates already found in database, nothing has been changed!")
		return
