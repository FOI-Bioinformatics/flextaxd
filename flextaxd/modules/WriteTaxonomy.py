#!/usr/bin/env python3 -c

'''
Read NCBI taxonomy dmp files (nodes or names) and holds a dictionary
'''

from .database.DatabaseConnection import DatabaseFunctions

class WriteTaxonomy(object):
	"""docstring for WriteTaxonomy."""
	def __init__(self, path, database=".taxonomydb",verbose=False,separator="\t|\t",minimal=False,prefix="names,nodes",desc=False,dbprogram=False):
		super(WriteTaxonomy, self).__init__()
		self.verbose = verbose
		self.database = DatabaseFunctions(database)
		self.path = path.rstrip("/")+"/"
		self.separator = separator
		self.prefix = prefix.split(",")
		self.dump_descriptions = desc
		### Allows a minimal output file with only nessesary fields default is NCBI id | name | empty | scientific name
		self.minimal = minimal
		if self.minimal:
			if dbprogram:
				print("# WARNING: dbprogram cannot be used in combination with dump_mini, parameter ignored! (use --dump)")
			self.dbprogram = False
			if separator != "\t|\t":
				self.separator = separator
			else:
				self.separator = "\t"
		else:
			self.dbprogram  = dbprogram
		self.parent = False ## Default print is NCBI structure with child in the first column


	def set_separator(self,sep):
		self.separator=sep
		return self.separator

	def set_order(self,ord):
		'''Set parent column (if parent or child is first in order)'''
		self.parent = ord

	def set_minimal(self):
		if self.verbose: print("Set minimal output!")
		self.minimal=True

	def set_prefix(self,prefix):
		self.prefix=prefix.split(",")
		return self.prefix

	def get_all(self, table, select="*"):
		QUERY = "SELECT {select} FROM {table}".format(select=select, table=table)
		if self.verbose: print(QUERY)
		return self.database.query(QUERY).fetchall()

	def get_links(self, table, select="child,parent,rank"):
		QUERY = "SELECT {select} FROM {table} JOIN (rank) on rank.rank_i = tree.rank_i".format(select=select, table=table)
		if self.verbose: print(QUERY)
		return self.database.query(QUERY).fetchall()

	def nodes(self):
		'''Write database tree to nodes.dmp'''
		if self.verbose: print('Write tree to: {}{}.dmp'.format(self.path,self.prefix[1]))
		with open('{}{}.dmp'.format(self.path,self.prefix[1]),"w") as outputfile:
			## Retrieve all links that exists in the database
			if self.dump_descriptions:
				self.nodeDict = self.database.get_nodes()
				print("child\tparent\trank", sep=self.separator, end="\n", file=outputfile)

			links = self.get_links('tree','child,parent,rank')
			for link in links:
				link = list(link)
				if self.parent:
					link[0],link[1] = link[1],link[0]
				if self.dump_descriptions:
					link[0],link[1] = self.nodeDict[link[0]],self.nodeDict[link[1]]
				if self.dbprogram == "kraken2":
					link = list(link)+["",""]  ## Make sure to add enough extra columns so that kraken2 does not trim away nessesary columns
				if not self.minimal:
					link = list(link)+[""]
				print(*link, sep=self.separator, end="\n", file=outputfile)

	def names(self):
		'''Write node annotations to names.dmp'''
		if self.verbose: print('Write annotations to: {}{}.dmp'.format(self.path,self.prefix[0]))
		end = "\n"
		if self.dbprogram == "ganon" or self.dbprogram == "krakenuniq":
			end = "\t|\n"
		with open('{}{}.dmp'.format(self.path,self.prefix[0]),"w") as outputfile:
			## Retrieve all nodes that exists in the database
			nodes = self.get_all('nodes', 'id,name')
			for node in nodes:
				if not self.minimal:
					node = list(node) + ["","scientific name"]
				if self.dbprogram == "kraken2":
					node += ["",""] ## Make sure to add enough extra columns so that kraken2 does not trim away nessesary columns
				print(*node, sep=self.separator, end=end, file=outputfile)
