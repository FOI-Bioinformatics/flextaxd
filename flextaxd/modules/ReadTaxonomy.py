#!/usr/bin/env python3 -c

'''
Read SILVA taxonomy dmp files (nodes or names) and holds a dictionary
'''

from .database.DatabaseConnection import DatabaseFunctions
from . import functions
from . import ranks
import logging
logger = logging.getLogger(__name__)

class InputError(Exception):
	"""Exception raised for errors in the input."""
	def __init__(self, message):
		#self.expression = expression
		self.message = message

class ReadTaxonomy(object):
	"""docstring for ReadTaxonomy."""
	'''The importers that add nodes by name work against the names the database already holds: an import into
	   a database that already has the tree would otherwise add a second node for every name and annotate the
	   genomes of the second run to that copy (issue #145). The NCBI reader works with taxids and turns this
	   off, it would only pay for a dictionary of millions of names.'''
	match_nodes = True

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
		self.mp = 0			## an input file whose levels are codes (d, p, g) rather than rank names
		self.skipped = []
		self.new_nodes = 0          ## nodes this run actually inserted, the ones already in the database are matched (issue #145)
		self.unknown_levels = set()   ## level codes of an input file that name no rank, reported once (issue #141)

		if self.match_nodes:
			'''the names already in the database, so that a name met again keeps the node it has (issue #145)'''
			self.taxonomy = self.database.get_nodes(col=2)
			logger.debug("Matching the {n} node names already in the database".format(n=len(self.taxonomy)))

		## Add base node
		self.root = self.add_node("root")
		self.add_rank("no rank")
		## Add basic link
		self.add_link(child=1, parent=1,rank="no rank")

	def set_mp(self,opt):
		self.mp = opt

	def write_skipped(self):
		'''Write and print rows/genomes skipped during parsing, so the user knows what was not added to the database'''
		if not self.skipped:
			return False
		if self.taxonomy_file:
			skipped_file = "{file}.skipped".format(file=self.taxonomy_file)
		else:
			skipped_file = "FlexTaxD.skipped"
		with open(skipped_file,"w") as of:
			for entry,reason in self.skipped:
				logger.warning("Skipped {entry}: {reason}".format(entry=entry,reason=reason))
				print("{entry}\t{reason}".format(entry=entry,reason=reason), end="\n", file=of)
		'''Keep only the number in the database, the rows themselves stay in the file written above'''
		self.database.save_report("skipped_import",
									self.taxonomy_file or "FlexTaxD",
									len(self.skipped),
									skipped_file)
		logger.warning("{n} rows/genomes skipped, listed in {file}".format(n=len(self.skipped),file=skipped_file))
		return skipped_file

	def zopen(self, path,*args, **kwargs):
		'''Redefine open to handle zipped files automatically, a mode containing "b" keeps the
		    returned data binary. The implementation is shared with the database builders'''
		return functions.zopen(path,*args,**kwargs)

	#@staticmethod
	def parse_taxonomy(self,treefile=False):
		'''Parse taxonomy information'''
		self.read_nodes(treefile=treefile)

	def set_separator(self,sep):
		'''Change separator from default (\t)'''
		self.sep = sep

	def add_rank(self, rank,mp=False):
		'''Insert the rank level of a rank name into the database and remember its id.

		A name that says nothing - no rank column, an empty field, the text None - gets no rank level: the
		link of that level stays without one, which the dmp export writes as an empty rank field (issue #134)
		and which is not the same thing as the name "no rank" that NCBI gives its unnamed taxids (issue #142).
		'''
		name = self.rank_name(rank,mp=mp)
		if ranks.unnamed(name):
			return None
		try:
			return self.rank[rank]
		except KeyError:
			pass
		rank_i = self.database.add_rank(name)
		self.rank[rank] = rank_i
		return rank_i

	def rank_name(self,rank,mp=False):
		'''The rank name a caller means: an MP/GTDB file carries level codes (g, s, sk) where the others
		   carry the rank name itself. A code outside that vocabulary is reported once, not on every row of
		   the file, and the levels it names are added without a rank (issue #141)'''
		if mp or self.mp:
			name = self.levelDict.get(rank)
			if name is None and rank is not None:
				if rank not in self.unknown_levels:
					self.unknown_levels.add(rank)
					logger.warning("Level {level} is not an MP level code: the nodes of that level are added without a rank".format(level=rank))
			return name
		return rank

	def add_link(self, child=None, parent=None,rank=None):
		'''Add relationship in tree

		The rank is a rank name, and it is added to the rank table when it is new: an input that brings a rank
		in keeps it (issue #142). A link without a rank name carries no rank level at all, and takes no rank id
		by default - rank_i 1 used to do that, which is the "no rank" level of every NCBI derived database.
		'''
		if not child and parent:
			raise InputError("A link requires both a child and a parent!")
		self.database.add_link(child=child,parent=parent,rank=self.add_rank(rank))
		self.ids+=1

	def rank_summary(self):
		'''Say which rank levels the database ended up with: the NCBI defined ones, the custom ones a taxonomy
		   file brought in, and the number of levels that carry no rank'''
		return ranks.summary(self.database, log=logger)

	def add_node(self, description,id=False):
		'''Add node to tree
			extend databaseFunction add_node function using self.taxonomy
			dict to keep track of which nodes were already added

			A name that is already known - added earlier in this run, or already in the database when the
			importer matches names - keeps the node it has: this project identifies nodes by name, and an
			import that met an existing name used to add a second node for it and put the annotations of the
			second run on that copy (issue #145). A file that brings its own taxid for a name that exists
			under another one is reported, the name decides.
		'''
		logger.debug("Node desc: {desc} id: {id}".format(desc=description,id=id))
		if description in self.taxonomy:
			found = self.taxonomy[description]
			if id and id != found:
				logger.warning("Node {node} is already in the database as {found}, the taxid {given} of the input is not used".format(node=description,found=found,given=id))
			return found
		self.taxid_base = self.database.add_node(description,id)
		self.taxonomy[description] = self.taxid_base
		self.new_nodes += 1
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
			file_parents = {}   ## node id and name for every id the file puts in a parent column
			file_children = set()   ## node ids the file links to something, to find the ids that are only parents
			for tree_row in _treefile:
				data = tree_row.strip().split(self.sep)
				logger.debug(data)
				if len(data) > 2:
					rank = data.pop().strip()
					if rank.strip() != "":
						self.add_rank(rank)
				if len(data) < 2:
					'''A row must hold both a child and a parent to define a link'''
					continue
				if swap:
					data[0],data[1] = data[1],data[0]
				'''Add nodes into database'''
				if data[0] == "":
					'''Check for empty rows'''
					continue
				'''a row whose parent and child are the same node is the root convention, not a level of
				   the file: the Brucella file of issue #162 opened with "Brucella melitensis	Brucella
				   melitensis	species" meaning "this is the top of my tree", and the row wrote a
				   self-link on a normal node instead; the top of the file is attached to the root after
				   the loop instead'''
				if data[0].strip() == data[1].strip():
					logger.warning("#Skipped: the taxonomy file links {node} to itself, the top of an imported tree is attached to the root".format(node=data[0]))
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
				file_parents[self.taxonomy[data[1].strip()]] = data[1].strip()
				file_children.add(self.taxonomy[data[0].strip()])
				self.length +=1
			'''the nodes the file puts above everything else: a file that names its own top node instead
			   of "root" ("Brucella melitensis" in issue #162, whose author tried the root self-link to
			   say it) would otherwise leave the whole imported tree hanging off nothing, invisible to
			   every walk of the validation. They get the database root as parent and no rank level, and
			   the import says so by name'''
			for node_id, node_name in file_parents.items():
				if node_id in file_children or node_id == 1:
					continue
				'''the root link may already be there: a second file imported into the same database
				   names its top after the node the first one brought in, the name matches the existing
				   id (#145), and the attach must not write the same root link a second time'''
				if self.database.select("SELECT COUNT(*) FROM tree WHERE parent = 1 AND child = ?;",
										params=(node_id,))[0][0] > 0:
					logger.info("{node} is the top of this file too and the root already links it, the two files join there".format(node=node_name))
					continue
				self.add_link(child=node_id, parent=1, rank=None)
				logger.warning("#Attached: {node} is the top of the imported tree and nothing links it, the database root becomes its parent".format(node=node_name))
			self.database.commit()

	def parse_genomeid2taxid(self,genomeid2taxid,reference=False):
		'''Parse file that annotates genome_id´s to nodes in the tree
			The file lists genome, node (name or taxid) and optionally the reference
			(refseq/genbank/nt) of that genome. The reference belongs to the row it was read with:
			a reference given as argument labels every genome and overrides that column, a row
			without one stays unlabelled.'''
		nodeDict = self.database.get_nodes()
		_ref = reference
		with self.zopen(genomeid2taxid,"rt") as f:
			# headers = f.readline().strip().split("\t") # This line was commented away since this genomeid2taxid is not documented to have a header
			for row in f:
				if row.strip() != "": ## If there are trailing empty lines in the file
					columns = row.strip().split("\t")
					if len(columns) < 2:
						columns = row.strip().split("    ") ## fallback shared with ModifyTree.update_annotations
					if len(columns) < 2:
						logger.warning("# WARNING: row without genome and taxid column skipped: {row}".format(row=row.strip()))
						self.skipped.append((row.strip(),"row has no genome and taxid columns"))
						continue
					genomeid,taxid = columns[0],columns[1]
					## the argument wins over the third column, and extra columns are ignored
					row_reference = _ref if _ref else (columns[2] if len(columns) > 2 else False)
					try:
						self.database.add_genome(genome=genomeid.strip(),_id=nodeDict[taxid.strip()],reference=row_reference)
					except KeyError:
						logger.warning("# WARNING: {taxid} not found in the database".format(taxid=taxid))
						self.skipped.append((genomeid.strip(),"taxid {taxid} not found in the database".format(taxid=taxid)))
			self.database.commit()
		self.write_skipped()
		return
