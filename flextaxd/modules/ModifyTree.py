#!/usr/bin/env python3

'''
Modify tree from database or source
'''

from .database.DatabaseConnection import ModifyFunctions
import logging,os
logger = logging.getLogger(__name__)
import math
#from .database.database import database

def progressBar(iterable, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
	"""
	Call in a loop to create terminal progress bar
	@params:
		iteration   - Required  : current iteration (Int)
		total       - Required  : total iterations (Int)
		prefix      - Optional  : prefix string (Str)
		suffix      - Optional  : suffix string (Str)
		decimals    - Optional  : positive number of decimals in percent complete (Int)
		length      - Optional  : character length of bar (Int)
		fill        - Optional  : bar fill character (Str)
		printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
	"""
	total = len(iterable)
	# Progress Bar Printing Function
	def printProgressBar (iteration):
		percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
		filledLength = int(length * iteration // total)
		bar = fill * filledLength + '-' * (length - filledLength)
		print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
	# Initial Call
	printProgressBar(0)
	# Update Progress Bar
	for i, item in enumerate(iterable):
		yield item
		printProgressBar(i + 1)
	# Print New Line on Complete
	print()

class InputError(Exception):
	"""Exception raised for errors in the input."""

	def __init__(self, message):
		#self.expression = expression
		self.message = message

class TreeError(Exception):
	"""Exception raised for errors in the Tree structure."""

	def __init__(self, message):
		#self.expression = expression
		self.message = message

class ModifyTree(object):
	"""docstring for ModifyTree."""
	def __init__(self, database=".taxonomydb", mod_database=False, mod_file=False, clean_database=False, purge_database=False,update_genomes=False, update_node_names=False,rename_node=False,separator="\t",verbose=False,parent=False,replace=False,**kwargs):
		super(ModifyTree, self).__init__()
		self.verbose = verbose
		logger.info("Modify Tree")
		self.sep = separator
		self.taxonomy = {}
		self.names = {}
		self.keep_rank = {} ## Place holder for parent rank

		self.parent=parent
		self.replace = replace
		self.ncbi_order = True
		self.mod_genomes =False

		self.fast_clean = True

		### Connect to or create database
		self.taxonomydb = ModifyFunctions(database,verbose=verbose)
		self.rank= self.taxonomydb.get_rank(col=2)
		## Save all nodes in the current database
		self.nodeDict = self.taxonomydb.get_nodes()
		self.clean = clean_database
		self.taxonomy_type = False
		self.do_not_delete_old = set()
		self.identical_nodes = set() #init needed. Needed in get_id in parse_modification in elif modfile
		try:
			if kwargs["taxonomy_type"] != "NCBI":  ## If taxonomy type in merge database is set to NCBI, keep ranks also above 9, default no.
				self.taxonomy_type = True
		except:
			pass
		if not self.clean:
			self.taxid_base = self.taxonomydb.get_taxid_base()

			'''Handle taxids'''
			logger.info("Taxid base updated!")
			self.taxid_base = self._taxfix(self.taxid_base)
			self.taxid_set = int(self.taxid_base)
			logger.debug("Taxid base: {taxidbase}".format(taxidbase=self.taxid_base))
		if mod_database:
			if not os.path.exists(mod_database):
				raise FileNotFoundError("The modification database was not found {path}".format(path=mod_database))
			self.moddb = ModifyFunctions(mod_database,verbose=verbose)
			self.identical_nodes = set(self.taxonomydb.get_nodes(col=2).keys()) & set(self.moddb.get_nodes(col=2).keys())
			self.dbmod_annotation = self.moddb.get_nodes(col=1)
			self.modsource = self.parse_modification(self.moddb,"database")
			
		elif mod_file:
			self.modsource = self.parse_modification(mod_file,"file")
		elif update_genomes or clean_database or update_node_names or rename_node or purge_database:
			pass
		else:
			raise InputError("No modification source could be found both mod_database and mod_file file are empty!")

	def _is_int(self,s):
		'''Small help function, test if value is int'''
		try:
			int(s)
			return True
		except ValueError:
			return False

	def _taxfix(self, x):
		'''Small function to fix taxid with NCBI, round up to nearest million in internal id'''
		return int(math.ceil(x / 1000000.0)) * 1000000

	def add_link(self, child=None, parent=None, rank="n"):
		'''Add relationship in tree'''
		if not child and parent:
			raise TreeError("A link requires both a child and a parent!")
		self.taxonomydb.add_link(child=child,parent=parent,rank=self.rank[rank])
		if not self.database.check_parent(child):
			raise TreeError("Node {child} has more than one parent!".format(child=child))
		self.ids+=1
		return

	def add_node(self, description):
		'''Add node to tree
			extend databaseFunction add_node function using self.taxonomy
			dict to keep track of which nodes were already added
		'''
		if self.taxid_set == self.taxid_base:
			self.taxid_base = self.taxonomydb.add_node(description,self.taxid_base)
			self.taxid_set = -1  #taxid base is not changing, make sure it´s not staying the same
		else:
			self.taxid_base = self.taxonomydb.add_node(description)
		self.taxonomy[description] = self.taxid_base
		return self.taxid_base

	def add_rank(self, rank):
		'''Insert rank into database and store into rank dictionary'''
		try:
			rank_i = self.rank[rank]
		except KeyError:
			rank_i = self.taxonomydb.add_rank(rank)
			self.rank[rank] = rank_i
			self.rank[rank_i] = rank
		#print(rank,rank_i)
		return rank_i

	def get_id(self, desc,ret=False,parent=False):
		'''Check if node exists if so return id, if not create a new node and return id'''
		try:
			## Check if parent node exists
			i = self.nodeDict[desc]
			if False: #desc == "Centipeda":
				print("taxonomydb parent {p}, {desc}".format(p=self.taxonomydb.get_parent(self.taxonomydb.get_id(desc)),desc=desc))
				print("moddb parent {p}, {i}".format(p=self.moddb.get_parent(self.moddb.get_id(desc)),i=i))
				print(self.taxonomydb.get_name(self.taxonomydb.get_parent(self.taxonomydb.get_id(desc))[0]))
				print(self.moddb.get_name(self.moddb.get_parent(self.moddb.get_id(desc))[0]))
			'''Due to identical names in different groups the parent of node always have to be checked (if seeming nessesary)'''
			'''This part of the code checks if the parent of the node ID is identical when existing, if not have to get unique ID'''
			#logger.info("Parent node already exists, check if duplicated entry, if so add _ to incoming node to allow separation {desc}".format(desc))
			#print(self.identical_nodes)
			#print(set([desc]))
			#print(set([desc]) & self.identical_nodes)
			if len(set([desc]) & self.identical_nodes) > 0 and parent:
				try:
					if self.taxonomydb.get_name(self.taxonomydb.get_parent(self.taxonomydb.get_id(desc))[0]) != self.moddb.get_name(self.moddb.get_parent(self.moddb.get_id(desc))[0]):
						#print("Old: ",desc, i)
						_i = i
						self.do_not_delete_old.add(i)
						desc += "_"
						try:
							i = self.nodeDict[desc]
						except KeyError:
							i = self.add_node(desc)
							logger.info("Identical nodes, add _ to make them unique: {Name}, new_taxid: {taxid}, dupid: {dupid}".format(Name=desc, taxid=i, dupid=_i))
							self.nodeDict[desc] = i

				except TypeError:
					i = self.add_node(desc)
					self.nodeDict[desc] = i

		except KeyError:
			#check if duplicate
			try:
				i = self.nodeDict[desc+"_"]
			except KeyError:
				'''Node does not exist, add node to the database'''
				i = self.add_node(desc)
				self.nodeDict[desc] = i
		if ret:
			return i,desc
		else:
			return i
			

	def _parse_new_links(self, parent=None,child=None,rank="no rank"):
		'''Help function for parse_mod_file, gets existing node id or adds new node'''
		## add new nodes
		if not child and parent:
			raise InputError("links requires both child and parent!")
		parent_i,parent = self.get_id(parent,ret=True)
		self.new_nodes.add(parent_i)
		child_i,child = self.get_id(child,ret=True,parent=True)
		self.new_nodes.add(child_i)
		'''If node has no level (no rank), check if parent database had a classified level. If so default add level'''
		try:
			level = self.parent_levels[int(child_i)]
			if True: #not self.taxonomy_type and level > 0: ## Higher ranks suggest wrong group (only 8 levels in bacteria), then skip
				level = False
			else:
				logger.info("Rank kept for {node}, {rank}".format(node=child, rank=level))
		except:
			level = False
		if not level:
			rank_i = self.add_rank(rank)
		else:
			rank_i = level
		#if parent == "CANJSG01" or parent == "Centipeda" or parent == "Selenomonadaceae":
		#	print("{parent}, {child}, {rank}, {parent_i}, {child_i},{rank_i}".format(parent=parent,child=child,rank=rank, parent_i=parent_i,child_i=child_i,rank_i=rank_i))
		self.new_links.add((parent_i,child_i,rank_i))
		return

	def database_mod(self,database,parent="root"):
		'''Handle database modification '''
		logger.debug(database)
		logger.info("Parse mod database... {db}".format(db=database))
		### Get translation dictionaries (from internal node index to description)
		self.dbmod_annotation = database.get_nodes(col=1)
		logger.debug(self.dbmod_annotation)
		self.dbmod_rank = database.get_rank(col=1)
		### Translate node ids between databases and add non existing nodes into current database
		for link in database.get_links(database.get_children(set([database.get_id(parent)]))):
			link = list(link)
			logger.debug(link)

			if link[0] == link[1]: ### This should only occur if the root node of the mod database is used link replace with
				if self.replace:
					if int(link[0]) != 1: ## Root to root
						self.new_links.add(self.taxonomydb.get_parent(self.taxonomydb.get_id(self.parent)))
						logger.debug("New root link: {}".format(self.new_links))
			else:
				try:
					parent,child,rank = self.dbmod_annotation[link[0]].strip(),self.dbmod_annotation[link[1]].strip(),self.dbmod_rank[link[2]]
					if parent == "Centipeda":
						logger.debug("New link: [{p}, {c}, {r}], {link}".format(p=parent,c=child,r=rank,link=link))
					if child == "Centipeda":
						logger.debug("New link: [{p}, {c}, {r}], {link}".format(p=parent,c=child,r=rank,link=link))
					logger.debug("New link: [{p}, {c}, {r}]".format(p=parent,c=child,r=rank))
					self._parse_new_links(parent=parent,child=child,rank=rank)

				except FileNotFoundError:
					logger.warning("Warning something is wrong! Check debug option")
					logger.debug(self.dbmod_annotation)
					logger.debug(self.dbmod_annotation[link[0]])
					logger.debug(self.dbmod_annotation[link[1]])
					exit()
		return database.get_genomes()

	def file_mod(self,modfile):
		'''Handle file modification'''
		logger.debug(modfile)
		logger.info("Parse modification file...")
		with open(modfile, "r") as f:
			headers = f.readline().strip().split(self.sep)
			logger.debug(headers)
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
					line +=["no rank"] ## Add no specified rank to line
				parent,child,rank = line
				try:
					rank_i = self.rank[rank]
				except:
					rank_i = self.add_rank(rank)
				if swap:
					## Make sure that the order between parent and child nodes are correct and consistent with the database
					child,parent = parent,child
				self._parse_new_links(parent=parent,child=child,rank=rank)
		return

	def translate(self,node,mod=False):
		'''from node id to node name'''
		n1,n2,r = node
		if mod:
			translate_dict = self.dbmod_annotation
		else:
			translate_dict = self.nodeDict
		n1=translate_dict[n1]
		n2=translate_dict[n2]
		return (n1,n2,r)

	def parse_modification(self, input,modtype="database"):
		'''Retrieve all links to update from an existing database'''
		## Retrieve all nodes annotated in the modified database
		logger.debug("Parent: {parent}".format(parent=self.parent))
		self.new_links = set()
		self.new_nodes = set()
		self.non_overlapping_old_links = set()
		self.existing_links = set()
		# ### Get the connecting link between the two databases
		self.parent_link = self.taxonomydb.get_parent(self.taxonomydb.get_id(self.parent))
		if not self.parent_link:
			raise InputError("The selected parent node ({parent}) could not be found in the source database!".format(parent=self.parent))
		self.existing_nodes = self.taxonomydb.get_children(set([self.taxonomydb.get_id(self.parent)])) ## - set([self.taxonomydb.get_id(self.parent)] )
		logger.info("{n} children to {parent}".format(n=len(self.existing_nodes),parent=self.parent))
		if len(self.existing_nodes) > 0:
			self.existing_links = set(self.taxonomydb.get_links(self.existing_nodes))
		parentlinks = set(self.taxonomydb.get_links([self.taxonomydb.get_id(self.parent)]))
		if modtype == "database":
			modparent = self.moddb.get_parent(self.moddb.get_id(self.parent))
			logger.info("{n} existing links to {parent} ({parentlinks}) ({modparent})".format(n=len(self.existing_links),parent=self.parent,parentlinks=parentlinks,modparent=modparent))
		self.parent_levels = self.keep_levels(self.existing_links | parentlinks | set(self.taxonomydb.get_links((self.existing_nodes & self.new_nodes))))
		if modtype == "database":
			self.mod_genomes = self.database_mod(input,self.parent)
		elif modtype == "file":
			self.file_mod(input)
		else:
			raise InputError("Wrong modification input database or file must be supplied")
		### get links from current database
		logger.info("Old links keep from duplicated taxonomy nodes {nodes}".format(nodes=self.do_not_delete_old))
		self.old_nodes = self.existing_nodes - (self.new_nodes)
		logger.info("nodes:")
		logger.info("old: {old}".format(old=len(self.old_nodes)))
		logger.info("new: {new}".format(new=len(self.new_nodes)))
		logger.info("ovl: {ovl}".format(ovl=len(self.existing_nodes & self.new_nodes)))
		if self.replace and len(self.existing_nodes) > 0:  ## remove nodes connected to old nodes that are not replaced
			if len((self.existing_nodes & self.new_nodes)) > 0:
				self.non_overlapping_old_links = set(self.taxonomydb.get_links((self.existing_nodes & self.new_nodes) - set([self.taxonomydb.get_id(self.parent)])))  ## Remove all links related to new nodes
				#logger.info("Existing nodes: {enodes}, New nodes: {nnodes}, Old nodes: {onodes}".format(enodes=len(self.existing_nodes),nnodes=len(self.new_nodes), onodes=len(self.old_nodes)))
				#logger.info("Existing nodes: {enodes}, Old nodes: {onodes}".format(enodes=self.existing_nodes, onodes=self.old_nodes))
				#logger.info("save links: {fill}".format(fill=self.non_overlapping_old_links))
		self.overlapping_links = self.existing_links & self.new_links## (links existing in both networks)
		self.old_links = self.existing_links - self.new_links
		logger.info("links:")
		logger.info("old: {old}".format(old=len(self.old_links)))
		logger.info("new: {new}".format(new=len(self.new_nodes)))
		logger.info("ovl: {ovl}".format(ovl=len(self.overlapping_links)))
		if self.replace:
			logger.debug("rm: {rm}".format(rm=len(self.non_overlapping_old_links)))
			logger.debug(self.non_overlapping_old_links)
		'''Get all genomes annotated to new nodes in existing database'''
		return True

	def update_annotations(self, genomeid2taxid, reference=False):
		'''Function that adds annotation of genome ids to nodes'''
		logger.info("Update genome to taxid annotations using {genomeid2taxid}".format(genomeid2taxid=genomeid2taxid))
		_ref = reference
		update = {
			"set_column": "id",
			"where_column": "genome",
			"set_value": "",
			"where": ""
		}
		#if _ref:
		#	update[""]
		updated = 0
		added = 0
		with open(genomeid2taxid) as f:
			for row in f:
				try:
					if len(row.strip().split(self.sep)) > 2:
						genome,name,reference = row.strip().split(self.sep)
					else:
						genome,name = row.strip().split(self.sep)
				except ValueError:
					genome,name = row.strip().split("    ")
				logger.debug("genome: {genome}, name: {name}".format(genome=genome,name=name))
				try:
					if not self._is_int(name.strip()):
						id = self.nodeDict[name.strip()]
					else:
						id = int(name.strip())
						## The input was already an index not a name, output warning or assume valid index?
						logger.debug("# WARNING: Input was an index not a name, make sure indexes match the current database!")

				except KeyError:
					logger.debug("# WARNING: there was no database entry for {name} annotation not updated for this entry!".format(name=name))
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
		self.taxonomydb.commit()
		gid = self.taxonomydb.get_genomes()
		logger.info("{added} added and {updated} genome annotations were updated!".format(added=added, updated=updated))
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
		notadded = 0
		'''Database has been updated, so the internal nodeDict needs to be updated'''
		if len(self.mod_genomes) >= 100:
			allupdates = []
		self.nodeDict = self.taxonomydb.get_nodes()
		for genome in self.mod_genomes:
			'''Translate incoming database node id to taxonomydb node id'''
			inc_taxid = self.dbmod_annotation[self.mod_genomes[genome]].strip()
			try:
				id,genomeid = self.nodeDict[inc_taxid],genome.strip()
			except KeyError:  ## taxid does not exist in receiving database, skip genome
				#logger.debug("Genome not added {genome}".format(genome=genome))
				notadded +=1
				continue
			if len(self.mod_genomes) < 100:
				update["set_value"] = id
				update["where"] = genomeid
				res = self.taxonomydb.update_genome(update)
				if self.taxonomydb.rowcount()!=0:
					if res:
						updated += 1
					else:
						added += 1
			else:
				allupdates.append("({id},'{genomeid}')".format(id=id,genomeid=genomeid))
		if len(self.mod_genomes) >= 100:
			update["set_value"] = "id"
			update["where_value"] = "genome"
			update["data"] = allupdates
			self.taxonomydb.multi_update(update,"genomes")
			logger.info("Genomes transfered to increase speed no statistics were calculated")
		else:
			if notadded > 0:
				logger.info("{notadded} genomes not added, taxonomy id does not exist in the receiving database".format(notadded=notadded))
			logger.info("{added} added and {updated} genome annotations were updated!".format(added=added, updated=updated))
		self.taxonomydb.commit()
		return

	def update_node_names(self, refdict):
		'''Function that adds annotation of genome ids to nodes'''
		logger.info("Update node names in the database from {refdict}".format(refdict=refdict))
		update = {
			"set_column": "name",
			"where_column": "name",
			"set_value": "",
			"where": ""
		}
		updated = 0
		added = 0
		with open(refdict) as f:
			for row in f:
				try:
					old_name,name = row.strip().split(self.sep)
				except ValueError:
					old_name,name = row.strip().split("    ")
				logger.debug("old_name: {old_name}, name: {name}".format(old_name=old_name,name=name))
				## If no exception occured add old_name
				update["set_value"] = name
				update["where"] = old_name.strip()
				res = self.taxonomydb.update_table(update,table="nodes")
				if self.taxonomydb.rowcount()!=0:
					if res:
						updated += 1
		self.taxonomydb.commit()
		gid = self.taxonomydb.get_genomes()
		logger.info("{updated} annotations were updated!".format(added=added, updated=updated))
		return

	def rename_node(self, data, table):
		'''Function that renames a node in the database'''
		try:
			print('[WARNING] Naivly attempting rename. This function needs validation-checks and warnings.')
			self.taxonomydb.update(data, table)
			self.taxonomydb.commit()
		except:
			print('Could not rename node. Make sure that it exists in the database')
		return

	def clean_database(self, ncbi=False):
		'''Function that removes all node and node paths without annotation'''
		logger.info("Fetch annotated nodes")
		self.annotated_nodes = set(self.taxonomydb.get_genomes(cols="genome,id").keys())
		an = len(self.annotated_nodes)
		logger.info("Annotated nodes: {an}".format(an=an))
		if an == 0 and not ncbi:
			raise InputError("Database has no annotations, the whole database would be cleaned")
		logger.info("Get all links in database")
		self.all_links = set(self.taxonomydb.get_links())
		logger.info("Get all nodes in database")
		self.all_nodes = set(self.taxonomydb.get_nodes(col=1).keys())
		'''Add parents to all nodes that may not have annotations'''
		logger.info("Retrieve all parents of annotated nodes")
		self.annotated_nodes |= self.taxonomydb.get_parents(self.annotated_nodes,find_all=True,ncbi=ncbi)
		logger.info("Parents added: {an}".format(an=len(self.annotated_nodes)-an))
		if ncbi:
			logger.info("Keep main nodes of the NCBI taxonomy (parents on level 3 and above)")
			self.keep = set(self.taxonomydb.get_children([1],maxdepth=3))  #set([self.nodeDict[node] for node in self.taxonomydb.get_children([1],maxdepth=1)])
			logger.info("Adding root levels {nlev}".format(nlev=len(self.keep-self.annotated_nodes)))
			self.annotated_nodes |= self.keep
		'''Get all links related to an annotated node and its parents'''
		self.annotated_links = set(self.taxonomydb.get_links(self.annotated_nodes,only_parents=True))
		self.clean_links = self.all_links - self.annotated_links
		self.clean_nodes = self.all_nodes - self.annotated_nodes
		logger.info("Links to remove {nlinks}".format(nlinks=len(self.clean_links)))
		logger.debug("Links remaining {nlinks}".format(nlinks=len(self.annotated_links)))
		logger.info("Nodes to remove {nnodes}".format(nnodes=len(self.clean_nodes)))
		logger.debug("Nodes remaining {nnodes}".format(nnodes=len(self.annotated_nodes)))
		logger.info("Clean annotations related to removed nodes")
		if len(self.clean_links) < len(self.all_links) and len(self.clean_links) > 0:
			logger.info("Cleaning {nlinks} links".format(nlinks=len(self.clean_links)))
			if self.fast_clean:
				self.taxonomydb.fast_delete_links(self.clean_links)
			else:
				logger.info("Cleaning tree (this will take a long time for large trees)")
				self.taxonomydb.delete_links(self.clean_links)
		if len(self.clean_nodes) < len(self.all_nodes) and len(self.clean_nodes) > 0:
			self.taxonomydb.delete_nodes(self.clean_nodes)
			if not ncbi:
				self.taxonomydb.delete_genomes(self.clean_nodes)
		if self.taxonomydb.validate_tree():
			self.taxonomydb.commit()

			logger.info("Vacuum database")
			self.taxonomydb.query("vacuum")
		logger.info("Database is cleaned!")

	def purge_database(self, genomes_list = None, force_genome_delete = False):
		''' This function takes as input a list of genomes that are missing in the genomes directory specified as --genomes_path in flextaxd-create.
            It gets the nodes for these genomes and then traverses the tables of the FlexTaxD-database, removing entries (e.g. parent-relations)
            that are exclusive to such nodes. It tests so that any parent node is not a parent for a node that is not deleted.
            Optionally, missing genomes that have no distinct nodes in the tree can be removed (force_genome_delete)
                This occurs when using GTDB taxonomy (full GTDB-database) but only using the representative genomes.
                The nodes of such datasets have multiple children, of which one is the representative genome and the other children can be removed.
                When this flag is applied, the function will remove missing entries from all tables wherever possible ("truly missing" datasets)
                in addition to naively removing all missing genomes (datasets that are missing due to exhaustive import of taxonomy information)
        '''
		logger.info('Entered purge_database')
		if not type(genomes_list) == set:
			genomes_list = set(genomes_list) # Convert input list to hashmap for lookup speed
		if genomes_list:
			logger.info('Have '+str(len(genomes_list))+' missing genomes, looking them up in the database...')
			# Get nodes and genomes from database
			logger.info('Fetching id,genome from database')
			genomes_nodes = self.taxonomydb.get_genomes(table='genomes',cols='id,genome') # get_gnomes returns dict of "genome"->"id"
			#print(next(iter(genomes_nodes)) )
			#/
			# Determine which nodes to keep (We have to do it this way. If using GTDB, one node ID may carry multiple genomes. Thus we must decide which nodes to keep as opposed which to delete.
			nodes_keep = set()
			#/ Speed implementation of match function, reverse dictionary, use union on keys vs genomes, fetch nodes that belong to genomes @Davod

			#genomes_nodes = self.taxonomydb.get_genomes(table='genomes',cols='genome,id')
			#print(next(iter(genomes_nodes)) )
			if True:
				remove_genome_keys = set(genomes_nodes) - genomes_list## Genomes missing and in database

			num_rows_before = self.taxonomydb.num_rows('genomes')
			self.taxonomydb.delete_genomes(nodes_keep,genomes=genomes_list,match_genome_only=True)

			
			self.clean_database(ncbi=True)
			num_rows_after = self.taxonomydb.num_rows('genomes')
			num_rows_deleted = num_rows_before - num_rows_after
			logger.info('Removed '+str(num_rows_deleted)+' node genomes')
			#print(list(nodes_keep)[1])
			# Replaced function
			#for genome,node in genomes_nodes.items():
			#	if not genome in genomes_list:  ## if not in missing
			#		nodes_keep.add(node)
			#/

			#/ End speed implementation

			return
			# Determine which nodes to delete (leaf-nodes [e.g. those that can hold a genome] in database) minus those that we decided to keep
			nodes_delete = set(genomes_nodes.values()).difference(nodes_keep3)

			logger.info('Number of genome-nodes from table "genomes": '+str(len(genomes_nodes)))
			logger.info('Number of genome-nodes to keep: '+str(len(nodes_keep3)))
			logger.info('Number of genome-nodes to delete: '+str(len(nodes_delete)))
			#/
			# Get parents for leaf-nodes (to-keep and to-delete). Determine which parents are delete-only
			## NOTE: for some reason the input child is also returned here. Its not intuitive, but it does not matter.
			parents_keep = set(self.taxonomydb.get_parents(nodes_keep,find_all=True))#,only_parents=True))
			parents_delete = set(self.taxonomydb.get_parents(nodes_delete,find_all=True))#,only_parents=True))
			parents_deleteOnly = parents_delete.difference(parents_keep)
			#/
			# Perform node deletion (leaf-nodes and exclusive parent-nodes)
			logger.info('Number of nodes to be deleted: '+str(len(parents_deleteOnly)))
			#@ delete links from table TREE (parent in _list_ and child in _list_)
			logger.info("Attempting delete of nodes in table 'tree'")
			self.taxonomydb.ambigious_delete_links(parents_deleteOnly)
			#@/
			#@ delete nodes from table NODES
			logger.info("Attempting delete of nodes in table 'nodes'")
			self.taxonomydb.delete_nodes(parents_deleteOnly)
			#@/
			#@ delete nodes from table GENOMES
			logger.info("Attempting delete of nodes in table 'genomes'")
			num_rows_before = self.taxonomydb.num_rows('genomes')
			self.taxonomydb.delete_genomes(parents_deleteOnly,genomes=list(remove_genome_keys),match_genome_only=True)
			num_rows_after = self.taxonomydb.num_rows('genomes')
			num_rows_deleted = num_rows_before - num_rows_after
			logger.info('Removed '+str(num_rows_deleted)+' node genomes')


			if num_rows_deleted != len(genomes_list) and not force_genome_delete:
				logger.warning('WARNING: Some nodes were not removed as they did not have exclusive nodes in the tree. This occurs for instance when building from GTDB taxonomy (that has taxonomy for their full database) but downloading only their representative set of genomes')
				logger.warning('To force the removal of these genomes, apply the flag --purge_database_force')
			if force_genome_delete: # if force delete specified, then remove all missing genomes
				logger.info("Force-deleting missing genomes from table 'genomes' N="+str(len(remove_genome_keys)))
				self.taxonomydb.delete_genomes('',genomes=remove_genome_keys,match_genome_only=True)
				logger.info('Removed '+str(len(genomes_list)-num_rows_deleted)+' genomes by force')
			#/

	def keep_levels(self, links):
		parent_levels = {}
		for c,p,l in links:
			if self.rank[l] != 1:
				parent_levels[p] = l
		return parent_levels

	def update_database(self):
		'''Update the database file'''
		if self.replace:
			logger.info("Clean up genomes annotated to child nodes from  {parent}".format(parent=self.parent))
			nodes = self.taxonomydb.get_children(set([self.taxonomydb.get_id(self.parent)])) | set([self.taxonomydb.get_id(self.parent)] )
			logger.debug(nodes)
			self.taxonomydb.delete_genomes(nodes)
			self.taxonomydb.query("vacuum") ## Actually remove the data from database
			if len(self.non_overlapping_old_links) + len(self.old_nodes) > 0:
				logger.info("Replace tree, deleting all nodes downstream of selected parent!")
			if len(self.old_links | self.non_overlapping_old_links-set(self.parent_link)) > 0:
				logger.debug("Delete links no longer valid!")
				'''Add function to keep taxonomy level from previous database (belongs to parent)'''
				self.taxonomydb.delete_links((self.old_links | self.non_overlapping_old_links)-set(self.parent_link))
			if len(self.old_nodes):
				logger.debug("Delete nodes!")
				self.taxonomydb.delete_nodes(self.old_nodes)
			self.taxonomydb.query("vacuum") ## Actually remove the data from database
		logger.debug("New links: [{links}]".format(links=self.new_links))
		links,nodes = self.taxonomydb.add_links(self.new_links)
		if len(links) + len(nodes) + len(self.non_overlapping_old_links) > 0:
			if self.replace and len(self.old_nodes) > 0: logger.info("Deleted {n} links and {n2} nodes that are no longer valid".format(n=len(self.old_links | self.non_overlapping_old_links-set(self.parent_link)),n2=len(self.old_nodes)))
			if len(self.new_nodes) > 1: logger.info("Adding {n} new nodes".format(n=len(nodes)))
			if len(links) > 1: logger.info("Adding {n} new links".format(n=len(links)))

			''' Commit changes (only commit once both deletion and addition of new nodes and links are completed!)'''
			self.taxonomydb.commit()
			self.nodeDict = self.taxonomydb.get_nodes()
			if self.mod_genomes:
				logger.info("Transfering genomeid2taxid annotation from incoming database")
				self.update_genomes()
				self.taxonomydb.query("vacuum")
		else:
			logger.info("All updates already found in database, nothing has been changed!")
		if self.taxid_set > 0:
			logger.info("Taxid base: {taxidbase}".format(taxidbase = self.taxid_base))
		else:
			logger.debug("Taxid base: {taxidbase}".format(taxidbase = self.taxid_base))
		logging.getLogger().setLevel(logging.INFO)
		self.taxonomydb.query("vacuum")
		logging.info("Validate modified database!")
		self.taxonomydb.validate_tree()
		return
