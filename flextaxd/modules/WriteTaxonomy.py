#!/usr/bin/env python3 -c

'''
Read NCBI taxonomy dmp files (nodes or names) and holds a dictionary
'''

from .database.DatabaseConnection import DatabaseFunctions
import logging
logger = logging.getLogger(__name__)

class WriteTaxonomy(object):
	"""docstring for WriteTaxonomy."""
	def __init__(self, path, database=".taxonomydb",separator="\t|\t",minimal=False,prefix="names,nodes",desc=False,dbprogram=None):
		super(WriteTaxonomy, self).__init__()
		self.database = DatabaseFunctions(database)
		logging.debug("Write settings: ")
		self.path = path.rstrip("/")+"/"
		logging.debug("Output path: {outdir}".format(outdir=self.path))
		self.separator = separator
		logging.debug("Output separator: '{separator}'".format(separator=self.separator))
		self.prefix = prefix.split(",")
		logging.debug("Prefix: nodes:{nodes} names:{names} ".format(nodes=self.prefix[1],names=self.prefix[0]))
		self.dump_descriptions = desc
		logging.debug("Add descriptions: {desc}".format(desc=self.dump_descriptions))
		### Allows a minimal output file with only nessesary fields default is NCBI id | name | empty | scientific name
		self.minimal = minimal
		if self.minimal:
			if dbprogram:
				logger.warning("# WARNING: dbprogram cannot be used in combination with dump_mini, parameter ignored! (use --dump)")
			self.dbprogram = None
			if separator != "\t|\t":
				self.separator = separator
			else:
				self.separator = "\t"
		else:
			self.dbprogram  = dbprogram
		if self.dbprogram: logging.debug("Output format for program {program}".format(program=self.dbprogram))
		self.link_order = False ## Default print is NCBI structure with child in the first column
		logging.debug("NCBI structure (child first): {parent}".format(parent=self.link_order))


	def set_separator(self,sep):
		self.separator=sep
		return self.separator

	def set_order(self,order):
		'''Set parent column (if parent or child is first in order)'''
		logging.debug("Changing order child first: {parent}".format(parent=self.link_order))
		self.link_order = order

	def set_minimal(self):
		logging.debug("Set minimal output to True!")
		self.dbprogram = None
		if self.separator != "\t|\t":
			self.separator = separator
		else:
			self.separator = "\t"
		self.minimal=True

	def set_prefix(self,prefix):
		self.prefix=prefix.split(",")
		logging.debug("Update output prefix nodes:{nodes} names:{names} ".format(nodes=self.prefix[1],names=self.prefix[0]))
		return self.prefix

	def get_all(self, table, select="*"):
		QUERY = "SELECT {select} FROM {table}".format(select=select, table=table)
		logging.debug(QUERY)
		return self.database.query(QUERY).fetchall()

	def get_links(self, table, select="child,parent,rank"):
		QUERY = "SELECT {select} FROM {table} JOIN (rank) on rank.rank_i = tree.rank_i".format(select=select, table=table)
		logging.debug(QUERY)
		return self.database.query(QUERY).fetchall()

	def nodes(self):
		'''Write database tree to nodes.dmp'''
		logging.info('Write tree to: {}{}.dmp'.format(self.path,self.prefix[1]))
		with open('{}{}.dmp'.format(self.path,self.prefix[1]),"w") as outputfile:
			## Retrieve all links that exists in the database
			if self.dump_descriptions:
				self.nodeDict = self.database.get_nodes()
				print("child\tparent\trank", sep=self.separator, end="\n", file=outputfile)

			links = self.get_links('tree','child,parent,rank')
			for link in links:
				link = list(link)
				if self.link_order:
					link[0],link[1] = link[1],link[0]
				if self.dump_descriptions:
					link[0],link[1] = self.nodeDict[link[0]],self.nodeDict[link[1]]
				if self.dbprogram in ["bracken","kraken2"]:
					link = list(link)+["-"]
				# if self.dbprogram == "kraken2":
				# 	link = list(link)+["",""] ## Make sure to add enough extra columns so that kraken2 does not trim away nessesary columns
				if not self.minimal:
					link = list(link)+[""]
				print(*link, sep=self.separator, end="\n", file=outputfile)

	def names(self):
		'''Write node annotations to names.dmp'''
		logging.info('Write annotations to: {}{}.dmp'.format(self.path,self.prefix[0]))
		end = "\n"
		if self.dbprogram in ["krakenuniq","kraken2"]:
			end = "\t|\n"
		with open('{}{}.dmp'.format(self.path,self.prefix[0]),"w") as outputfile:
			## Retrieve all nodes that exists in the database
			nodes = self.get_all('nodes', 'id,name')
			empty = ""
			for node in nodes:
				if not self.minimal:
					empty = ""
					if self.dbprogram == "bracken":
						empty = "-"
					node = list(node) + [empty,"scientific name"]
				print(*node, sep=self.separator, end=end, file=outputfile)
