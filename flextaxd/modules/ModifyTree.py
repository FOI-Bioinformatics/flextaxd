#!/usr/bin/env python3

'''
Modify tree from database or source
'''

from .database.DatabaseConnection import ModifyFunctions, MAX_BIND_PARAMS
from . import ranks, functions
import logging,os,datetime
logger = logging.getLogger(__name__)
import math,time
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
		self.mod_genomes =False

		self.fast_clean = True
		
		## Hidden option specific request, replace name of root node, with root of file.
		try:
			self.replace_parent = kwargs["replace_parent"]
		except:
			logger.debug("kwarg replace_parent not present")
			self.replace_parent = False
		### Connect to or create database
		self.taxonomydb = ModifyFunctions(database,verbose=verbose)
		self.rank= self.taxonomydb.get_rank(col=2)
		## Save all nodes in the current database, only needed when nodes are looked up by name
		if mod_database or mod_file or update_node_names or rename_node or update_genomes:
			self.nodeDict = self.taxonomydb.get_nodes()
		else:
			self.nodeDict = {}
			logger.debug("Node names are not used by this modification, nodeDict is not loaded")
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
		self.modfile = mod_file	
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

	def add_link(self, child=None, parent=None, rank=None):
		'''Add relationship in tree. The rank is a rank name: it is added to the rank table when it is new, and
		   a link that nothing names carries no rank level at all (issue #142) - it used to fall back on the
		   rank "n", which is the "no rank" level of an NCBI derived database.'''
		if not child and parent:
			raise TreeError("A link requires both a child and a parent!")
		self.taxonomydb.add_link(child=child,parent=parent,rank=self.add_rank(rank))
		if self.taxonomydb.check_parent():
			raise TreeError("Node {child} has more than one parent!".format(child=child))
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
		'''Insert rank into database and store into rank dictionary. A rank that names nothing (no rank
		   column, an empty field, the text None) gets no rank level: the link stays without one (issue #142)'''
		if ranks.unnamed(rank):
			return None
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
			

	def _parse_new_links(self, parent=None,child=None,rank=None):
		'''Help function for parse_mod_file, gets existing node id or adds new node. Without a rank name the
		   link is added without a rank level, it is no longer put on "no rank" (issue #142)'''
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
			if not self.taxonomy_type and level > 0: ## Higher ranks suggest wrong group (only 8 levels in bacteria), then skip
				level = False
			else:
				logger.info("Rank kept for {node}, {rank}".format(node=child, rank=level))
		except:
			level = False
		if not level:
			rank_i = self.add_rank(rank)
		else:
			rank_i = level
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
		untranslatable = 0
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
					logger.debug("New link: [{p}, {c}, {r}]".format(p=parent,c=child,r=rank))
					self._parse_new_links(parent=parent,child=child,rank=rank)

				except KeyError as error:
					'''a link that names a node or a rank the source database does not describe is a broken
					   source database: name what is missing and skip that link, the bare KeyError traceback
					   this used to die with said nothing about the node (issue #153)'''
					logger.warning("Link {link} of {db} names {missing}, which that database has no description or rank for, the link is skipped".format(
									link=link,db=database.database,missing=error))
					untranslatable += 1
		if untranslatable:
			logger.warning("{n} links of {db} could not be translated, they name nodes or ranks that database does not describe".format(
							n=untranslatable,db=database.database))
		return database.get_genomes()

	def file_mod(self,modfile,replace_parent=False):
		'''Handle file modification'''
		logger.debug(modfile)
		logger.info("Parse modification file...")
		'''The gz sniffing opener of the builders and the taxonomy readers: a modification file may be
		   compressed, and plain open() would read the compressed bytes as text (issue #127)'''
		with functions.zopen(modfile, "r") as f:
			headers = f.readline().strip().split(self.sep)
			if replace_parent:
				new_root = f.readline().strip().split(self.sep)
				change_root = new_root[0]
				new_root[0] = self.parent
				parent,child,rank = new_root
				self._parse_new_links(parent=parent,child=child,rank=rank)
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
					line +=[None] ## the line names no rank: the level stays without one, it is not the "no rank" name NCBI uses (issue #142)
				parent,child,rank = line
				if replace_parent:
					if parent == change_root:
						parent = self.parent
				try:
					rank_i = self.rank[rank]
				except:
					rank_i = self.add_rank(rank)
				if swap:
					## Make sure that the order between parent and child nodes are correct and consistent with the database
					child,parent = parent,child
				'''the root convention is the only link that points a node at itself: a row linking a
				   node to itself hides that node's subtree from the walks of the validation, and the
				   Brucella modification file of issue #162 started with exactly one such row'''
				if str(parent).strip() == str(child).strip():
					logger.warning("#Skipped: the modification file links {node} to itself, only the root links to itself".format(node=child))
					continue
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
		logger.info(self.parent_link)
		if not self.parent_link:
			raise InputError("The selected parent node ({parent}) could not be found in the source database!".format(parent=self.parent))
		
		self.existing_nodes = self.taxonomydb.get_children(set([self.taxonomydb.get_id(self.parent)])) ## - set([self.taxonomydb.get_id(self.parent)] )
		if modtype == "file":
			self.existing_nodes	-= set([self.taxonomydb.get_id(self.parent)])
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
			self.file_mod(input,self.replace_parent)
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
		if modtype == "file":
			logger.info("Modification type is file...")
			if self.replace_parent:  ## Update the name of the parent to match with parent of "incoming" file
				logger.info("Replace parent option initiated,  change the name of parent in the database to allow merge")
				with functions.zopen(self.modfile, "r") as f:
					headers = f.readline().strip().split(self.sep)
					branchroot = f.readline().strip().split(self.sep)
					new_name = branchroot[0]
					#logger.info("New root name: {nn}".format(nn=new_name))
				logger.info("New name of parent, changing from {oldp}, to {newp}".format(oldp=self.parent, newp=new_name))
				self._name_update(self.parent,new_name)
				self.taxonomydb.commit()  # Make sure change is done to active database
				time.sleep(1) ## Sleep one second to make change happening in the database
				self.parent = new_name
				logger.info("Update complete, new self.parent = {sp}".format(sp=self.parent))
				logger.info("Check node id:")
				logger.info(set([self.taxonomydb.get_id(self.parent)]))
		return True

	def update_annotations(self, genomeid2taxid, reference=False):
		'''Function that adds annotation of genome ids to nodes
			The file lists genome, node (name or index) and optionally the reference
			(refseq/genbank/nt) of every genome, which is the format of the genome annotation dump
			(WriteTaxonomy.dump_genome_annotations), so that a dump can be edited and read back
			without losing the labels. A reference given as argument labels every genome and
			overrides that column, a file without the column is labelled with the given reference.'''
		logger.info("Update genome to taxid annotations using {genomeid2taxid}".format(genomeid2taxid=genomeid2taxid))
		_ref = reference
		update = {
			"set_column": "id",
			"where_column": "genome",
			"set_value": "",
			"where": ""
		}
		reference_update = {
			"set_column": "reference",
			"where_column": "genome",
			"set_value": "",
			"where": ""
		}
		updated = 0
		added = 0
		skipped = []
		'''The rows are collected before anything is written. Two things come out of that: the taxids the
		   file names by number are asked of the database once instead of trusted (issue #73), and the
		   annotations can be written in batches the way the merge route does it, one statement per batch
		   instead of one to three per genome, which is about 800000 statements for the 394932 genomes of
		   the GTDB case (issue #126).'''
		rows = []
		'''the same gz sniffing opener the build route uses for this file (parse_genomeid2taxid): a
		   "--genomeid2taxid file.gz" was readable during a build and binary garbage through the
		   modify route (issue #127)'''
		with functions.zopen(genomeid2taxid, "r") as f:
			for row in f:
				'''The reference belongs to the row it was read with: the argument of the function
				   labels every genome, otherwise the third column of the file does'''
				genome_reference = _ref
				try:
					if len(row.strip().split(self.sep)) > 2:
						genome,name,column_reference = row.strip().split(self.sep)
						if not _ref:
							genome_reference = column_reference
					else:
						genome,name = row.strip().split(self.sep)
				except ValueError:
					genome,name = row.strip().split("    ")
				logger.debug("genome: {genome}, name: {name}".format(genome=genome,name=name))
				'''the node column names a node, and a node may be named by digits alone: the Brucella
				   tree of issue #161 is full of names like "21158" and "1330", and forcing a digit
				   column to be a taxid skipped 25 of its 91 genomes in a run that had just annotated
				   all of them through the build route. The name is tried first; only when no node
				   carries it does a pure number become a taxid, which the database is asked about
				   below (issue #73)'''
				node_name = name.strip()
				if node_name in self.nodeDict:
					id = self.nodeDict[node_name]
					if self._is_int(node_name) and int(node_name) in self.nodeDict and int(node_name) != id:
						logger.warning("{name} is both the name of node {found} and the taxid of another node: the name won, the annotation goes to {found}".format(name=node_name,found=id))
					rows.append((id,genome.strip(),genome_reference))
				elif self._is_int(node_name):
					## a number naming no node is taken as a taxid, which the database is asked about below
					rows.append((int(node_name),genome.strip(),genome_reference))
				else:
					logger.debug("# WARNING: there was no database entry for {name} annotation not updated for this entry!".format(name=name))
					skipped.append((genome.strip(),"no database entry for node {name}".format(name=node_name)))

		unknown_taxids = self._unknown_taxids([id for id,_,_ in rows if isinstance(id,int)])
		if unknown_taxids:
			'''a taxid the file names that the database has no node for: the genome keeps the annotation
			   it had and is listed, instead of being pointed at a node that does not exist and vanishing
			   from the dump and the classifier build (issue #73)'''
			logger.warning("{n} taxids named by {file} have no node in the database, the first are: {first}".format(
				n=len(unknown_taxids),file=genomeid2taxid,
				first=", ".join(str(taxid) for taxid in sorted(unknown_taxids)[:5])))
			for id,genome,_ in rows:
				if id in unknown_taxids:
					skipped.append((genome,"taxid {id} has no node in the database".format(id=id)))
			rows = [row for row in rows if row[0] not in unknown_taxids]

		if len(rows) >= 100:
			updated,added = self._batch_annotations(rows)
		else:
			for id,genome,genome_reference in rows:
				update["set_value"] = id
				update["where"] = genome
				res = self.taxonomydb.update_genome(update)
				if self.taxonomydb.rowcount()!=0:
					if res:
						updated += 1
					else:
						added += 1
				if genome_reference:
					reference_update["set_value"] = genome_reference
					reference_update["where"] = genome
					self.taxonomydb.update_genome(reference_update)
		self.taxonomydb.commit()
		if skipped:
			skipped_file = "{file}.skipped".format(file=genomeid2taxid)
			with open(skipped_file,"w") as of:
				for entry,reason in skipped:
					logger.warning("Skipped genome {entry}: {reason}".format(entry=entry,reason=reason))
					print("{entry}\t{reason}".format(entry=entry,reason=reason), end="\n", file=of)
			self.taxonomydb.save_report("skipped_import",
										genomeid2taxid,
										len(skipped),
										skipped_file)
			logger.warning("{n} genome annotations skipped, listed in {file}".format(n=len(skipped),file=skipped_file))
		gid = self.taxonomydb.get_genomes()
		logger.info("{added} added and {updated} genome annotations were updated!".format(added=added, updated=updated))
		return

	def _unknown_taxids(self, taxids):
		'''Which of the taxids an annotation file names by number the database has no node for. Asked in
				chunks of placeholders, the limit of the connection, so that a file of half a million rows
				costs a few hundred index lookups instead of half a million queries.

		------
		Parameters:
			list - the taxids read as numbers, may repeat
		------
		Returns:
			set - the taxids with no node in the nodes table
		'''
		unique = sorted(set(taxids))
		if not unique:
			return set()
		known = set()
		for start in range(0, len(unique), MAX_BIND_PARAMS):
			chunk = unique[start:start+MAX_BIND_PARAMS]
			rows = self.taxonomydb.select("SELECT id FROM nodes WHERE id IN ({})".format(
				",".join(["?"]*len(chunk))), params=tuple(chunk))
			known.update(int(id) for id, in rows)
		return set(unique) - known

	def _batch_annotations(self, rows, batch=20000):
		'''Write the collected annotations through the batched update the merge route already uses: one
				statement per batch of literal rows instead of one to three statements per genome (issue
				#126). A row whose genome is not in the table yet adds it, the rest update the node and the
				reference of the row that is there, an empty reference keeping the label that is stored
				(multi_update, issue #116).

		------
		Parameters:
			list - rows as (node id, genome, reference)
			int - how many rows a single statement carries, the statements of the merge route were
				measured at 14 MB for 400000 rows, which is far past what this route needs
		------
		Returns:
			tuple - (updated, added)
		'''
		update = {
			"set_column": "id",
			"where_column": "genome",
			"columns": ["id","genome","reference"]
		}
		updated = 0
		added = 0
		for start in range(0, len(rows), batch):
			chunk = rows[start:start+batch]
			update["data"] = ["({id},'{genome}',{reference})".format(
				id=id,
				genome=str(genome).replace("'","''"),
				reference="NULL" if not reference else "'{ref}'".format(ref=str(reference).replace("'","''")))
				for id,genome,reference in chunk]
			num_before = self.taxonomydb.num_rows('genomes')
			if not self.taxonomydb.multi_update(update,"genomes"):
				logger.error("{n} genome annotations could not be written, the database keeps the annotations it had".format(n=len(chunk)))
				continue
			num_added = self.taxonomydb.num_rows('genomes') - num_before
			added += num_added
			updated += len(chunk) - num_added
		logger.info("{n} genome annotations written in batches of {batch}".format(n=len(rows),batch=batch))
		return updated, added

	def update_genomes(self):
		'''When a database is supplied as source for the update genome annotations from that database needs to be transfered to the taxonomydb,
		   including the reference (refseq/genbank/nt) label every incoming genome carries'''
		update = {
			"set_column": "id",
			"where_column": "genome",
			"set_value": "",
			"where": ""
		}
		reference_update = {
			"set_column": "reference",
			"where_column": "genome",
			"set_value": "",
			"where": ""
		}
		updated = 0
		added = 0
		notadded = 0
		notadded_genomes = []
		'''The reference of the incoming genomes is read from the source database, unlabelled
		   source genomes keep the reference the receiving database already has'''
		genome_references = {}
		if getattr(self,"moddb",False):
			for genome, reference in self.moddb.select("SELECT genome,reference FROM genomes"):
				genome_references[genome] = reference
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
				notadded_genomes.append("{genome}\ttaxid {taxid} is not in the receiving database".format(genome=genome,taxid=inc_taxid))
				continue
			reference = genome_references.get(genome) or genome_references.get(genomeid)
			if len(self.mod_genomes) < 100:
				update["set_value"] = id
				update["where"] = genomeid
				res = self.taxonomydb.update_genome(update)
				if self.taxonomydb.rowcount()!=0:
					if res:
						updated += 1
					else:
						added += 1
				if reference:
					reference_update["set_value"] = reference
					reference_update["where"] = genomeid
					self.taxonomydb.update_genome(reference_update)
			else:
				reference_value = "NULL" if not reference else "'{ref}'".format(ref=str(reference).replace("'","''"))
				allupdates.append("({id},'{genomeid}',{reference})".format(id=id,genomeid=str(genomeid).replace("'","''"),reference=reference_value))
		if len(self.mod_genomes) >= 100:
			update["columns"] = ["id","genome","reference"]
			update["data"] = allupdates
			num_rows_before = self.taxonomydb.num_rows('genomes')
			if self.taxonomydb.multi_update(update,"genomes"):
				added = self.taxonomydb.num_rows('genomes') - num_rows_before
				updated = len(allupdates) - added
				logger.info("{added} added and {updated} genome annotations were transfered in one batch".format(added=added,updated=updated))
			else:
				logger.error("{n} genome annotations could not be transfered, the receiving database keeps the annotations it had".format(n=len(allupdates)))
		else:
			logger.info("{added} added and {updated} genome annotations were updated!".format(added=added, updated=updated))
		if notadded > 0:
			logger.info("{notadded} genomes not added, taxonomy id does not exist in the receiving database".format(notadded=notadded))
			'''The genomes that could not be transfered are listed in a file next to the receiving
			   database: the database itself keeps the number, not the individual genomes'''
			notadded_file = "{database}.notadded".format(database=self.taxonomydb.database)
			with open(notadded_file,"w") as of:
				for entry in notadded_genomes:
					print(entry, end="\n", file=of)
			source = self.moddb.database if getattr(self,"moddb",False) else "modification"
			self.taxonomydb.save_report("skipped_import",source,len(notadded_genomes),notadded_file)
			logger.warning("{n} genomes were not added, they are listed in {file}".format(n=notadded,file=notadded_file))
		self.taxonomydb.commit()
		return

	def _name_update(self, old_name,new_name,updated=0):
		'''Funciton that updates the name'''
		## If no exception occured add old_name
		update = {
			"set_column": "name",
			"where_column": "name",
			"set_value": "",
			"where": ""
		}
		update["set_value"] = new_name
		update["where"] = old_name.strip()
		res = self.taxonomydb.update_table(update,table="nodes")
		if self.taxonomydb.rowcount()!=0:
			if res:
				updated += 1
		return updated

	def update_node_names(self, refdict):
		'''Function that adds annotation of genome ids to nodes'''
		logger.info("Update node names in the database from {refdict}".format(refdict=refdict))
		updated = 0
		added = 0
		with functions.zopen(refdict, "r") as f:
			for row in f:
				try:
					old_name,name = row.strip().split(self.sep)
				except ValueError:
					old_name,name = row.strip().split("    ")
				logger.debug("old_name: {old_name}, name: {name}".format(old_name=old_name,name=name))
				
				updated = self._name_update(old_name,name,updated)
				
		self.taxonomydb.commit()
		gid = self.taxonomydb.get_genomes()
		logger.info("{updated} annotations were updated!".format(added=added, updated=updated))
		return gid

	#def rename_node(self, data, table):
	#	'''Function that renames a node in the database'''
	#	try:
	#		print('[WARNING] Naivly attempting rename. This function needs validation-checks and warnings.')
	#		self.taxonomydb.update(data, table)
	#		self.taxonomydb.commit()
	#	except:
	#		print('Could not rename node. Make sure that it exists in the database')
	#	return

	def clean_database(self, ncbi=False):
		'''Function that removes all node and node paths without annotation.
		    The nodes that survive are the annotated nodes, all their parents and, for NCBI
		    databases, the top levels of the tree. Which rows to keep is decided inside the
		    database, so the size of the database does not decide the amount of data read'''
		logger.info("Fetch annotated nodes")
		self.annotated_nodes = set(self.taxonomydb.get_genome_nodes())
		an = len(self.annotated_nodes)
		logger.info("Annotated nodes: {an}".format(an=an))
		if an == 0 and not ncbi:
			raise InputError("Database has no annotations, the whole database would be cleaned")
		## An NCBI clean of an unannotated database is allowed (used to slim a database before
		## adding the taxonomy that is actually wanted), but the size of the result is reported
		'''Add parents to all nodes that may not have annotations'''
		logger.info("Retrieve all parents of annotated nodes")
		parents = self.taxonomydb.get_ancestor_nodes(self.annotated_nodes, max_depth=50)
		if parents is False:
			logger.error("Parent nodes could not be fetched, cleaning is aborted!")
			return
		self.annotated_nodes |= set(parents)
		logger.info("Parents added: {an}".format(an=len(set(parents))-an))
		if ncbi:
			logger.info("Keep main nodes of the NCBI taxonomy (parents on level 3 and above)")
			children = self.taxonomydb.get_descendant_nodes([1], max_depth=3)
			if children is False:
				logger.error("Root levels of the taxonomy could not be fetched, cleaning is aborted!")
				return
			self.keep = set(children) | set([1])  ## the root node itself is kept, its children link to it
			logger.info("Adding root levels {nlev}".format(nlev=len(self.keep-self.annotated_nodes)))
			self.annotated_nodes |= self.keep
		num_nodes = self.taxonomydb.num_rows('nodes')
		num_links = self.taxonomydb.num_rows('tree')
		if len(self.annotated_nodes) == 0:
			logger.info("No nodes would survive the cleaning, nothing is removed")
			return
		if len(self.annotated_nodes) >= num_nodes:
			logger.info("All {nnodes} nodes are annotated or parents of annotations, nothing to remove".format(nnodes=num_nodes))
			return
		logger.info("Nodes to keep {nnodes} of {all}".format(nnodes=len(self.annotated_nodes),all=num_nodes))
		if an == 0:
			logger.warning("Database has no genome annotations, this clean keeps only the top levels of the taxonomy: {keep} of {all} nodes survive".format(keep=len(self.annotated_nodes),all=num_nodes))
		keep_table = self.taxonomydb._temp_value_table("ftx_temp_keep_nodes", self.annotated_nodes)
		if not keep_table:
			logger.error("The list of nodes to keep could not be built, cleaning is aborted!")
			return
		'''Remove all links leading to nodes that are not kept, then the nodes themselves'''
		logger.info("Cleaning {nlinks} links".format(nlinks=num_links - len(self.annotated_nodes)))
		links_removed = self.taxonomydb.delete_rows_not_in('tree','child',keep_table,hold=True)
		if links_removed is False:
			self.taxonomydb.drop_temp_tables([keep_table])
			logger.error("Links could not be removed, database changes are not committed and cleaning is aborted!")
			return
		nodes_removed = self.taxonomydb.delete_rows_not_in('nodes','id',keep_table,hold=True)
		if nodes_removed is False:
			self.taxonomydb.drop_temp_tables([keep_table])
			logger.error("Nodes could not be removed, database changes are not committed and cleaning is aborted!")
			return
		if not ncbi:
			genomes_removed = self.taxonomydb.delete_rows_not_in('genomes','id',keep_table,hold=True)
			if genomes_removed is False:
				self.taxonomydb.drop_temp_tables([keep_table])
				logger.error("Genome annotations could not be removed, cleaning is aborted!")
				return
		self.taxonomydb.drop_temp_tables([keep_table])
		self.taxonomydb.commit()
		if self.taxonomydb.validate_tree():
			self.taxonomydb.commit()

			logger.info("Vacuum database")
			self.taxonomydb.query("vacuum")
		logger.info("Database is cleaned!")

	def purge_database(self, genomes_list = None, genome_path = False):
		''' This function takes as input a list of genomes that are missing in the genomes directory
		    given as --purge_database. All annotations of those genomes are removed, matched on the
		    genome name, which also covers genomes that share their node with a genome that is present:
		    with GTDB taxonomy (the full GTDB-database) and only the representative dataset on disk,
		    one species node carries the representative and all other genomes of that species.
		    Nodes and links that end up with no annotated genome below them are removed by the
		    database clean that runs afterwards, which is the same function as --clean_database.
		    The removed option --purge_database_force needs no replacement: removing missing genomes
		    by name is what the purge always does now, and removing the nodes of the genomes that
		    have no annotated genome downstream is what the clean always does.
		    The disk report of --purge_database (genome_path) describes the state before the purge:
		    every genome it listed as annotated without a file has just been removed, so its numbers
		    are updated to the state after it. The file is compacted at the end: the clean only
		    vacuums when it removed nodes itself, while a GTDB purge removes genome rows below nodes
		    that keep their representative and therefore no nodes at all.
        '''
		logger.info('Entered purge_database')
		if not type(genomes_list) == set:
			genomes_list = set(genomes_list) # Convert input list to hashmap for lookup speed
		if genomes_list:
			logger.info('Have '+str(len(genomes_list))+' missing genomes, removing their annotations...')
			num_rows_before = self.taxonomydb.num_rows('genomes')
			genomes_removed = self.taxonomydb.delete_genomes(set(),genomes=genomes_list,match_genome_only=True)
			num_rows_after = self.taxonomydb.num_rows('genomes')
			num_rows_deleted = num_rows_before - num_rows_after
			logger.info('Removed '+str(num_rows_deleted)+' node genomes')
			if not genomes_removed:
				logger.error("Genome annotations could not be removed, database is left unchanged and cleaning is skipped!")
				return
			if num_rows_deleted < len(genomes_list):
				logger.warning(str(len(genomes_list)-num_rows_deleted)+' missing genomes had no genome annotation in the database')
			if genome_path:
				'''The scan that produced this list stored its findings for that directory, the purge
				   just removed all of them except the ones that had no annotation to remove'''
				missing_after = len(genomes_list) - num_rows_deleted
				self.taxonomydb.set_meta("genomes_path",genome_path)
				self.taxonomydb.set_meta("report_date",datetime.datetime.now().isoformat(timespec="seconds"))
				self.taxonomydb.set_meta("genomes_missing_on_disk",str(missing_after))
				self.taxonomydb.save_report("missing_on_disk",genome_path,missing_after)
			if num_rows_after == 0:
				logger.warning("The purge removed every genome annotation in the database, the taxonomy is left untouched")
				self._compact_database()
				return
			try:
				self.clean_database(ncbi=True)
			except InputError as e:
				logger.warning("Cleaning after the purge was refused: {err}".format(err=str(e)))
			self._compact_database()
			return

	def _compact_database(self):
		'''Give the pages the purge freed back to the file: the vacuum of clean_database only runs
			when that function had nodes to remove, which a purge of a GTDB style database has not'''
		logger.info("Vacuum database after the purge")
		'''A vacuum cannot run inside an open transaction, and the clean can return early without
		   committing when it finds nothing to remove'''
		self.taxonomydb.commit()
		res = self.taxonomydb.query("vacuum")
		if isinstance(res,Exception):
			logger.warning("The database file could not be compacted: {err}".format(err=str(res)))
			
	def keep_levels(self, links):
		'''The levels the branch being modified already has, so that a node which is re-added can keep the rank
		   of the node above it. A level that is not a level - "no rank", or a link carrying no rank at all -
			must not be handed down: the test used to compare the rank NAME with the id 1, so it never fired
			and every node inherited whatever its parent carried (issue #140).'''
		parent_levels = {}
		for c,p,l in links:
			if not ranks.is_no_rank(self.rank.get(l)):
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
