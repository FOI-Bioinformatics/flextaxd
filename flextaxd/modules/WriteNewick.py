#!/usr/bin/env python3 -c

'''
Read NCBI taxonomy dmp files (nodes or names) and holds a dictionary
'''

from .database.DatabaseConnection import DatabaseConnection
import logging
logger = logging.getLogger(__name__)

class WriteNewick(object):
	"""docstring for WriteNewick."""
	def __init__(self, database=".taxonomydb",verbose=False):
		super(WriteNewick, self).__init__()
		self.verbose = verbose
		self.database = self.set_database(self.open_database(database))
		logger.info("Load WriteNewick")

	def get(self, table, select="*"):
		QUERY = "SELECT {select} FROM {table}".format(select=select, table=table)
		return self.database.query(QUERY).fetchall()

	def open_database(self, database):
		'''Create connection to the database'''
		return DatabaseConnection(database)

	def set_database(self,database_obj):
		'''Set a database object as WriteNewick database'''
		self.database = database_obj

	def print_newick_tree(self):
		'''Print a newick formatted tree from the database of a species'''
		logger.info("Print Tree")
		tree = self.get("tree")
		nodes = self.get("nodes")
