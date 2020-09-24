import sys
import os
import sqlite3
import logging
logger = logging.getLogger(__name__)

class ConnectionError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class NameError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class TreeError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class DatabaseConnection(object):
	"""docstring for DatabaseConnection"""
	def __init__(self, database, verbose=False):
		super().__init__()
		self.verbose = verbose
		self.database = database
		BASE_DIR = os.path.dirname(os.path.abspath(__file__))  ## Retrieve path
		if not os.path.exists(self.database):
			if self.verbose:
				logger.debug("python {path}/CreateDatabase.py {database}".format(path=BASE_DIR,database=self.database))
			os.system("python {path}/CreateDatabase.py {database}".format(path=BASE_DIR,database=self.database))
		try: ## If database connection already exists
			self.conn
		except AttributeError:
			logger.debug("Connecting to {database}".format(database=self.database))
			self.conn = self.connect(self.database)
			self.cursor = self.create_cursor(self.conn)

	def __str__(self):
		return "Object of class DatabaseConnection, connected to {database}".format(database=self.database)

	def __repr__(self):
		return "DatabaseConnection()"

	def set_verbose(self,val):
		self.verbose = val

	def connect(self,database):
		'''Create database connection'''
		try:
			self.conn = sqlite3.connect(database)
			logger.info("{database} opened successfully.".format(database=database))
			return self.conn
		except Exception as e:
			logger.debug.write(str(e))
		raise ConnectionError("Count not connect to the database {database} see above message for details!".format(database=database))

	def create_cursor(self,conn):
		'''Create a db cursor'''
		try:
			self.cursor = conn.cursor()
			logger.debug("cursor created.")
			return self.cursor
		except Exception as e:
			sys.stderr.write(str(e))
		raise ConnectionError("Count not create db cursor to {database}".format(database=database))

	def commit(self):
		self.conn.commit()

	def query(self,query,insert_val = False, cursor=False,error=False):
		'''    The query function controls  the error handling of sqlite3 execute'''
		res = False
		if not cursor:
			cursor = self.cursor
		try:
			if insert_val:
				res = cursor.execute(query,insert_val)
				if error:
					return res
				return cursor.lastrowid
			else:
				return cursor.execute(query)
		except Exception as e:
			if "UNIQUE constraint failed" not in str(e):
				## UNIQUE constraint is an accepted error as it keeps multiple edges from being added
				logger.warning("Error in DatabaseConnection query")
				logger.debug(query)
				logger.debug("Insert val: ", insert_val)
				sys.stderr.write(str(e))
			return(e)

	def insert(self,data,table):
		'''Insert function
				data is a dictionary with keys matching
				table columns with respective value to be inserted
		'''
		INSERT_QUERY = '''
			INSERT INTO {table}({columns})
			VALUES ({values})
		'''
		columns = data.keys()
		values = tuple([data[col] for col in columns])
		insertStr = INSERT_QUERY.format(
				columns=",".join(columns),
				values=','.join(["?" for x in values]),
				table=table
		)
		return self.query(insertStr,insert_val=values)

	def update(self,data,table):
		'''Update function requires table column which column to identify row with and value to replace'''
		UPDATE_QUERY = '''
			UPDATE {table}
			SET {set_column} = ?
			WHERE {where_column} = ?
		'''.format(table=table,set_column=data["set_column"],where_column=data["where_column"])
		udata = tuple([data["set_value"],data["where"]])
		logger.debug("{q} {d}".format(q=UPDATE_QUERY,d=udata))
		res = self.query(UPDATE_QUERY,udata,error=True)
		if self.rowcount() == 0:
			res = self.insert({data["set_column"]:data["set_value"],data["where_column"]:data["where"]}, table)
			return False
		return True

	def delete(self,nodes,table):
		'''Deleting a node should make sure all related data is also removed'''
		DELETE_QUERY = '''
				DELETE FROM {table}
					WHERE id in({values})
		'''.format(table=table,
					values=','.join(["?" for x in nodes])
		)
		return self.query(DELETE_QUERY,insert_val=values)

	def rowcount(self):
		return self.cursor.rowcount

class DatabaseFunctions(DatabaseConnection):
	"""docstring for DatabaseFunctions."""
	def __init__(self, database, verbose=False):
		super().__init__(database, verbose)
		logger.debug("Load DatabaseFunctions")

	'''Validate tree function'''
	def validate_tree(self):
		'''This function validates the tree structure in the databases
			1. All nodes must only have one parent (check_parent)
			2. All edges must be attatched to the tree (match tree links and edges)
			3. All nodes must be attatched to the tree (match tree nodes and nodes)
			4. All nodes from links must have a node description (match the number of nodes from links with number of annotated nodes)
		'''
		logger.info("Get database nodes")
		self.nodes = self.get_nodes(col=1)					## get all nodes
		logger.info("Get database edges")
		self.edges = self.get_links()					## Get all links
		logger.info("Get children")
		self.tree_connections = self.get_children([1])  ## Get all chilren from root
		logger.info("Get tree edges")
		self.tree_links = self.get_links(self.tree_connections) ## get all links from rooted tree
		logger.info("Get tree nodes")
		self.tree_nodes = set()
		for edge in self.edges:
			self.tree_nodes.add(edge[0])
			self.tree_nodes.add(edge[1])
		logger.info("Validate parents")
		for n in self.nodes.keys():
			self.check_parent(n)
		stats = """Tree statistics
					Nodes: {nodes}
					Links: {links}
					Tree: n({tnodes}), l({tlinks})
					LinkNodes: {link_nodes}
					Parent_ok: True
					""".format(nodes = len(self.nodes),
								links=len(self.edges),
								tnodes=len(self.tree_connections),
								tlinks = len(self.tree_links),
								link_nodes = len(self.tree_nodes)
		)
		logger.info(stats)
		nodeset = set(self.nodes.keys())
		if len(nodeset - self.tree_connections) != 0:
			raise TreeError("The number of nodes and the number nodes under root does not match!")
		if len(set(self.edges) - set(self.tree_links)) != 0:
			raise TreeError("The number of edges and the number edges under root does not match!")
		if len(self.tree_nodes - nodeset) != 0:
			raise TreeError("The number of annotated nodes does not match with the number of nodes connected to edges!")
		if len(self.tree_connections - nodeset) != 0:
			raise TreeError("There are nodes in the database not connected to the tree")
		logger.info("Validation OK!")

	def check_parent(self,name):
		'''check tree structure parent
			make sure no tree node has multiple parents
		'''
		#QUERY = '''SELECT parent,child,rank FROM tree LEFT JOIN rank on (tree.rank_i = rank.rank_i) WHERE child = "{node}"'''.format(node=name)
		try:
			QUERY = '''SELECT parent,child,rank_i FROM tree WHERE child = {node}'''.format(node=name)
			logger.debug(QUERY)
			res = self.query(QUERY).fetchall()
			if len(res) == 1:
				return True
			elif len(res) == 0:
				logger.info(QUERY)
				logger.info(res)
				raise TreeError("Node: {node} has no parent".format(node=name))
			if res[0] == (1,1,1):
				return True
			logger.info(QUERY)
			logger.info(res)
		except AttributeError:
			logger.info(QUERY)
			logger.info(res)
		raise TreeError("Node: {node} has more than one parent!".format(node=name))

	'''Get functions of class'''
	def get_all(self, database=False, table=False):
		'''Get full table from table'''
		QUERY = '''SELECT * FROM {table}'''.format(table)
		return self.query(QUERY).fetchall()

	def get_taxid_base(self):
		'''Fetch the next incremental node from the current database'''
		QUERY = "SELECT MAX(id) AS max_id FROM nodes"
		res = self.query(QUERY).fetchone()[0]
		if res == None:
			res = 1
		else:
			res += 1
		return res

	def get_genomes(self, database=False,limit=0,table="genomes",cols="id,genome"):
		'''Get the list of genomes in the database'''
		## This is a many to many relation, so all genomes has to be put in a set for each taxonomy id
		genomeDict = {}
		QUERY = '''SELECT {cols} FROM {table}'''.format(cols=cols,table=table)
		if limit > 0:
			QUERY += " LIMIT {limit}".format(limit=limit)
		for id,genome in self.query(QUERY).fetchall():
			genomeDict[genome] = id
		return genomeDict

	def get_nodes(self, database=False,col=False):
		'''Retrieve the whole node info table of the database to decrease the number of database calls!'''
		nodeDict = {}
		QUERY = '''SELECT id,name FROM nodes'''
		logger.debug(QUERY)
		if not database:
			database = self.database
		for node in self.query(QUERY).fetchall():
			nodeDict[node[0]] = node[1]
			if col == 1:
				continue
			nodeDict[node[1]] = node[0]
		return nodeDict

	def get_links(self, nodes=False,database=False,swap=False,only_parents=False):
		'''This function returns all links in the given database'''
		order = ["parent","child"]
		if swap:
			order[1],order[0] = order[0],order[1]

		if only_parents:
			QUERY = '''SELECT {order},rank_i FROM tree WHERE child in ({nodes})'''.format(nodes=",".join(map(str,nodes)),order=",".join(order))
		elif nodes:
			QUERY = '''SELECT {order},rank_i FROM tree WHERE parent in ({nodes}) OR child in ({nodes})'''.format(nodes=",".join(map(str,nodes)),order=",".join(order))
		else:
			QUERY = '''SELECT {order},rank_i FROM tree'''.format(order=",".join(order))
		logger.debug(QUERY)
		if not database:
			database = self.database
		links = self.query(QUERY).fetchall()
		return links

	'''Add functions of class'''
	def add_node(self, description, id=False, table="nodes"):
		'''Add node to tree'''
		if description.strip() == "":  ## for some reason an empty node has been given, do not add empty nodes
			return False
		info = { "name": description }
		### If ID is supplied skip autoincrement and add specific ID
		if id:
			info["id"] = id
		taxid_base = self.insert(info, table=table)
		#logger.debug("node added [description {}, taxid base {}]: ".format(info, taxid_base))
		return taxid_base

	def add_rank(self, rank ,id=False):
		'''Add node to tree'''
		if rank == None:
			raise AttributeError("Rank cannot be None")
		info = { "rank": rank }
		### If ID is supplied skip autoincrement and add specific ID
		#logger.debug("rank added: {info}".format(info=info))
		if id:
			info["id"] = id
		rank_id = self.insert(info, table="rank")
		#logger.debug("rank added: {info}".format(info=info))
		return rank_id

	def add_link(self, child=None, parent=None, rank=1, table="tree"):
		'''Add relationship in tree'''
		info = {
			"child": child,
			"parent": parent,
			"rank_i": rank
		}
		logger.debug("link added:  child {}, parent {}, rank {} ".format(child,parent,rank))
		return self.insert(info, table="tree")

	def add_genome(self, genome, _id=False):
		'''Add genome annotation to nodes'''
		info = {
			"genome": genome
		}
		if _id:
			info["id"] = _id
		return self.insert(info, table="genomes")

	def add_links(self,links, table="tree",hold=False):
		'''Add links from a list to tree'''
		added_links = []
		nodes = set()
		for parent,child,rank in links:
			#logger.debug("{} {} rank: {} link added!".format(parent,child,rank))
			res = self.add_link(child,parent,rank,table=table)
			### Check if the link already exist in the database, this overlap may occur when a large new branch is added
			if "UNIQUE constraint failed" not in str(res):
				added_links.append([parent,child,rank])
				nodes.add(parent)
				nodes.add(child)
		## Commit changes
		if not hold:
			self.commit()
		return added_links,nodes

	def add_nodes(self,nodes, table="tree",hold=False):
		'''Add nodes from a list of nodes'''
		added_nodes = 0
		for node in nodes:
			res = self.add_node(node)
			### Check if the link already exist in the database, this overlap may occur when a large new branch is added
			if "UNIQUE constraint failed" not in str(res):
				added_nodes +=1
		## Commit changes
		if not hold:
			self.commit()
		return added_nodes

	'''Delete functions of class'''
	def delete_links(self,links, table="tree",hold=False):
		'''This function deletes all links given in links'''
		QUERY = "DELETE FROM {table} WHERE parent = {parent} AND child = {child}"
		logger.debug("Deleting {nlinks} links!".format(nlinks=len(links)))
		logger.debug(QUERY.format(table=table,parent="",child=""))
		for parent,child,rank in links:
			logger.debug("{}-{} rank: {} deleted!".format(parent,child,rank))
			res = self.query(QUERY.format(table=table, parent=parent, child=child))
		## Commit changes
		if not hold:
			logger.debug("Commit changes!")
			self.commit()

	def delete_nodes(self, nodes, table="nodes",hold=False):
		'''This function deletes all nodes given in nodes'''
		QUERY = "DELETE FROM {table} WHERE id = {node}"
		logger.debug("Deleting {nnodes} nodes!".format(nnodes=len(nodes)))
		logger.debug(QUERY.format(table=table,node=""))
		for node in nodes:
			logger.debug("Delete node {node}".format(node=node))
			res = self.query(QUERY.format(table=table, node=node))
		## Commit changes
		if not hold:
			logger.debug("Commit changes!")
			self.commit()

	def num_rows(self,table):
		'''Return the number of rows in a table'''
		QUERY = '''SELECT Count(*) FROM {table}'''
		logger.debug(QUERY)
		return self.query(QUERY.format(table=table)).fetchall()[0][0]


class ModifyFunctions(DatabaseFunctions):
	"""ModifyFunctions adds a few important functions only nessesary when modifying a database"""
	def __init__(self, database, verbose=False):
		super().__init__(database, verbose)
		logger.debug("Load ModifyFunctions")

	def get_rank(self,col=1):
		'''Get rank index from database'''
		QUERY = "SELECT rank_i,rank FROM rank"
		logger.debug(QUERY)
		rankDict = {}
		for rank in self.query(QUERY).fetchall():
			rankDict[rank[0]] = rank[1]
			if col == 1:
				continue
			rankDict[rank[1]] = rank[0]
		return rankDict

	def get_children(self,parents,children=set(),level=0,maxdepth=50):
		'''Get all children from a parent'''
		QUERY = '''SELECT child FROM tree WHERE parent in({nodes})'''.format(nodes=",".join(map(str,list(parents))))
		logger.debug(QUERY)
		res = self.query(QUERY).fetchall()
		if (len(res) + len(children)) != 0:
			children = set([child_i[0] for child_i in res])
			if level < maxdepth:
				children |= self.get_children(parents=children,level=level+1,maxdepth=maxdepth)
		return children

	def get_parent(self,name):
		'''return parent'''
		#QUERY = '''SELECT parent,child,rank FROM tree LEFT JOIN rank on (tree.rank_i = rank.rank_i) WHERE child = "{node}"'''.format(node=name)
		QUERY = '''SELECT parent,child,rank_i FROM tree WHERE child = "{node}"'''.format(node=name)
		logger.debug(QUERY)
		res = self.query(QUERY).fetchone()
		return res

	def get_parents(self,name,parents=set(),depth=0):
		'''Get all parents until root'''
		ld = depth+1
		if isinstance(name, int):
			name = [name]
		QUERY = '''SELECT parent,child FROM tree WHERE child in ({node})'''.format(node=",".join(map(str,name)))
		logger.debug(QUERY)
		res = self.query(QUERY).fetchall()
		try:
			res = res[0]
		except IndexError:
			logger.warning("WARNING: parent could not be found for node {res} \n{query}".format(query=QUERY, res=name))
			return set()
		if res[0] == res[1]:  ## Special for root
			return set([res[0]])
		else:
			## Add all parents
			try:
				parents |= self.get_parents([int(res[0])],depth=ld)
			except RecursionError:
				print(res, name)
		## Add current node
		parents |= set([int(res[0])])
		return parents

	def get_id(self,name):
		'''get node id from name'''
		QUERY = '''SELECT id FROM nodes WHERE name = "{node}"'''.format(node=name)
		try:
			res = self.query(QUERY).fetchone()[0]
		except TypeError:
			logger.debug(QUERY)
			raise NameError("Name not found in the database! {name}".format(name=name))
		return res

	def update_genome(self,data):
		'''Add genome annotation to nodes'''
		return self.update(data, table="genomes")
