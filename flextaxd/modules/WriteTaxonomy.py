#!/usr/bin/env python3 -c

'''
Read NCBI taxonomy dmp files (nodes or names) and holds a dictionary
'''

from .database.DatabaseConnection import DatabaseFunctions
import logging
logger = logging.getLogger(__name__)

class DumpFormatError(Exception):
	'''The dump that was asked for cannot be written: the columns a program reads and the minimal dump
	   exclude each other, and writing the minimal file for a program leaves it with taxids and no names
	   (issue #136)'''
	def __init__(self, message):
		super(DumpFormatError, self).__init__(message)
		self.message = message

class WriteTaxonomy(object):
	"""docstring for WriteTaxonomy."""
	def __init__(self, path, database=".taxonomydb",separator="\t|\t",minimal=False,prefix="names,nodes",desc=False,dbprogram=None,dump_genomes=False,**kwargs):
		super(WriteTaxonomy, self).__init__()
		self.database = DatabaseFunctions(database)
		logging.debug("Write settings: ")
		self.path = path.rstrip("/")+"/"
		logging.debug("Output path: {outdir}".format(outdir=self.path))
		self.separator = separator
		logging.debug("Output separator: '{separator}'".format(separator=self.separator))
		self.prefix = self.check_prefix(prefix)
		logging.debug("Prefix: nodes:{nodes} names:{names} ".format(nodes=self.prefix[1],names=self.prefix[0]))
		self.dump_descriptions = desc
		logging.debug("Add descriptions: {desc}".format(desc=self.dump_descriptions))
		### Allows a minimal output file with only nessesary fields default is NCBI id | name | empty | scientific name
		self.updated = False
		self.minimal = minimal
		'''both settings are kept: the writers of the dmp files refuse the combination of a minimal layout and
		   a program format, while the other dumps (the genome lists) do not depend on the columns'''
		self.dbprogram  = dbprogram
		if self.minimal:
			if separator != "\t|\t":
				self.separator = separator
			else:
				self.separator = "\t"
		if self.dbprogram: logging.debug("Output format for program {program}".format(program=self.dbprogram))
		self.link_order = False ## Default print is NCBI structure with child in the first column
		logging.debug("NCBI structure (child first): {parent}".format(parent=self.link_order))

	def dump_taxid_map(self):
		logging.info("Dump taxid maps")
		with open('{}{}.map'.format(self.path,"seqid2taxid"),"w") as outputfile,open('{}{}.txt'.format(self.path,"prelim_map"),"w") as prelim:
			genomes = self.get_all('genomes', 'genome,id')
			for genome in genomes:
				print(*genome, sep="\t", end="\n", file=outputfile)
				genome=tuple(["TAXID"]) + genome
				print(*genome, sep="\t", end="\n", file=prelim)
		return

	def dump_genomes(self):
		'''Write the list of annotated genomes to a file'''
		with open('{}{}.dmp'.format(self.path,"genomes"),"w") as outputfile:
			genomes = self.get_all('genomes', 'genome,reference', sort="reference")
			for genome in genomes:
				print(*genome, sep="\t", end="\n", file=outputfile)
		return

	def dump_genome_annotations(self, sort="reference"):
		'''Dump all genomes, including their taxonomy reference (will work as input file for genomeid2taxid)'''
		select = "genome,name,reference"
		QUERY = "SELECT {select} FROM genomes JOIN nodes ON nodes.id=genomes.id".format(select=select)
		if sort:
			QUERY += " ORDER BY {col} DESC".format(col=sort)
		logging.debug(QUERY)
		with open('{}{}.dmp'.format(self.path,"genomes"),"w") as outputfile:
			genomes = self.database.query(QUERY).fetchall()
			for genome in genomes:
				print(*genome, sep="\t", end="\n", file=outputfile)
		return

	def set_separator(self,sep):
		self.separator=sep
		return self.separator

	def set_order(self,order):
		'''Set parent column (if parent or child is first in order)'''
		logging.debug("Changing order child first: {parent}".format(parent=self.link_order))
		self.link_order = order

	def set_minimal(self):
		'''Switch the writer to the minimal layout: tab separated, without the trailing columns. Those columns
		   are what kraken2 and krakenuniq read for the name of a taxid, so the two cannot be combined
		   (issue #136)'''
		logging.debug("Set minimal output to True!")
		if self.dbprogram:
			raise DumpFormatError("the minimal dump cannot be written for --dbprogram {program}: it drops the trailing columns {program} reads, leave either the minimal layout or the program out".format(program=self.dbprogram))
		self.dbprogram = None
		if self.separator == "\t|\t":
			self.separator = "\t"
		self.minimal=True

	def check_prefix(self,prefix):
		'''The prefix must hold exactly two names: the names prefix and the nodes prefix'''
		parts = prefix.split(",")
		if len(parts) != 2:
			raise ValueError("The dump prefix requires two comma separated names (names,nodes), got: {prefix}".format(prefix=prefix))
		return parts

	def set_prefix(self,prefix):
		self.prefix=self.check_prefix(prefix)
		logging.debug("Update output prefix nodes:{nodes} names:{names} ".format(nodes=self.prefix[1],names=self.prefix[0]))
		return self.prefix

	def get_all(self, table, select="*", sort=False):
		QUERY = "SELECT {select} FROM {table}".format(select=select, table=table)
		if sort:
			QUERY += " ORDER BY {col} DESC".format(col=sort)
		logging.debug(QUERY)
		return self.database.query(QUERY).fetchall()

	def get_links(self, table, select="child,parent,rank"):
		'''The links of a table with the name of their rank level beside them.

			The rank is taken with a LEFT JOIN, so every link of the tree is returned: a database holds far
			more nodes than it has rank levels (only the official ranks have a row, and a custom taxonomy
			deliberately leaves the levels below species without one), and a link that is not returned never
			reaches the exported files, which takes its whole subtree away from the program that reads them
			(issue #134). A link without a rank row comes back with the rank as None, which the writers turn
			into an empty rank field.'''
		QUERY = "SELECT {select} FROM {table} LEFT JOIN rank on rank.rank_i = {table}.rank_i".format(select=select, table=table)
		logging.debug(QUERY)
		return self.database.query(QUERY).fetchall()

	def unique_indexes(self):
		'''Check duplicated indexes and give them unique IDs before print'''
		QUERY = "SELECT child FROM tree GROUP BY child HAVING count(parent) > 1"  ## Thanks to andrewjmc@github for this suggestion
		child_w_dpi = self.database.query(QUERY).fetchall()  ## Fetch all conflicting links and give them unique index before printing
		child_w_dpi = [x[0] for x in child_w_dpi]
		lmax = 10000000
		if len(child_w_dpi) > 0:
			self.nodeDict = self.database.get_nodes(col=1)
			lmax = max(self.nodeDict)
		return lmax,child_w_dpi

	def check_program_layout(self):
		'''The trailing columns names.dmp and nodes.dmp write for kraken2, krakenuniq and bracken are exactly
		   the columns the minimal layout drops. Measured with kraken2 2.17.2: the same 8 sequences get indexed
		   either way, but the database built from the minimal dump lists every taxid without a name, and the
		   run says nothing - so the combination is refused (issue #136).'''
		if self.minimal and self.dbprogram:
			raise DumpFormatError("--dump_mini cannot be combined with --dbprogram {program}: the minimal dump drops the columns {program} reads, use --dump --dbprogram {program} or leave the program out".format(program=self.dbprogram))

	def nodes(self):
		'''Write database tree to nodes.dmp'''
		self.check_program_layout()
		logging.info('Write tree to: {}{}.dmp'.format(self.path,self.prefix[1]))
		child_w_dpi = []
		checklist = {}
		written = 0
		unranked = 0
		with open('{}{}.dmp'.format(self.path,self.prefix[1]),"w") as outputfile:
			## Retrieve all links that exists in the database
			if self.dump_descriptions:
				self.nodeDict = self.database.get_nodes()
				print("child\tparent\trank", sep=self.separator, end="\n", file=outputfile)
			else:
				lmax,child_w_dpi = self.unique_indexes()
			links = self.get_links('tree','child,parent,rank')
			for link in links:
				link = list(link)
				if link[2] is None:
					'''A level with no rank row keeps an empty rank field: only the official ranks have a
						row in the database and the levels a custom taxonomy adds below species deliberately
						have none, and the program reading the file numbers those itself (kraken2 gives them
						S1, S2 ...). Printing the value as it comes back would name the rank "None".'''
					link[2] = ""
					unranked += 1
				'''If child with duplicate index, make sure link has unique index, match with node'''
				if link[0] in child_w_dpi:
					try:
						checklist[link[0]] += 1
						name = self.nodeDict[link[0]]
						'''Update node index to uniqe'''
						link[0] = lmax+100+len(checklist.keys())+checklist[link[0]]
						self.nodeDict[link[0]] = name
						self.updated = True
					except KeyError:
						checklist[link[0]] = 0
						## Do nothing
				if self.link_order:
					link[0],link[1] = link[1],link[0]
				if self.dump_descriptions:
					link[0],link[1] = self.nodeDict[link[0]],self.nodeDict[link[1]]
				if self.dbprogram in ["bracken"]:
					link = list(link)+["-"]
				if self.dbprogram == "kraken2":
					link = list(link)+["",""] ## Make sure to add enough extra columns so that kraken2 does not trim away nessesary columns
				if not self.minimal:
					link = list(link)+[""]
				print(*link, sep=self.separator, end="\n", file=outputfile)
				written += 1
		'''Every link of the tree must reach the file: a node that is missing from it takes its whole subtree
			away from the program that reads the dump, and nothing else in a run says so (issue #134)'''
		counted = self.database.query("SELECT COUNT(*) FROM tree")
		if isinstance(counted, Exception):
			logger.warning("The number of tree links written could not be checked: {err}".format(err=counted))
			return
		links = counted.fetchone()[0]
		if written != links:
			logger.warning("{w} of {n} tree links were written to {file}: {m} links did not reach the file".format(
							w=written,n=links,file=self.path+self.prefix[1]+".dmp",m=links-written))
		else:
			logger.info("{w} tree links written to {file}, {u} of them without a rank level".format(
						w=written,u=unranked,file=self.path+self.prefix[1]+".dmp"))

	def names(self):
		'''Write node annotations to names.dmp'''
		self.check_program_layout()
		logging.info('Write annotations to: {}{}.dmp'.format(self.path,self.prefix[0]))
		end = "\n"
		if self.dbprogram in ["krakenuniq","kraken2"]:
			end = "\t|\n"
		with open('{}{}.dmp'.format(self.path,self.prefix[0]),"w") as outputfile:
			## Retrieve all nodes that exists in the database
			if self.updated:
				nodes = self.nodeDict.items()
			else:
				nodes = self.get_all('nodes', 'id,name')
			empty = ""
			for node in nodes:
				if not self.minimal:
					empty = ""
					if self.dbprogram == "bracken":
						empty = "-"
					node = list(node) + [empty,"scientific name"]
				print(*node, sep=self.separator, end=end, file=outputfile)
