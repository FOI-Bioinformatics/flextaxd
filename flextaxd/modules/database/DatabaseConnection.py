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
		'''Set the verbose level of the DatabaseConnection object

		------
		Returns
			boolean - response value'''
		self.verbose = val
		return self.verbose

	def connect(self,database):
		'''Create database connection

		------
		Returns
			connection object (sqlite3)
		'''
		try:
			self.conn = sqlite3.connect(database)
			logger.info("{database} opened successfully.".format(database=database))
			return self.conn
		except Exception as e:
			sys.stderr.write(str(e))
		raise ConnectionError("Count not connect to the database {database} see above message for details!".format(database=database))

	def create_cursor(self,conn):
		'''Create a db cursor

		------
		Returns
			cursor - return database pointer
		'''
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
		'''The query function is a wrapper around sqlite3 execute to form responses related to the request

		Returns
			int - 	if insert: lastrowid
			list -	if request: sqlite3 execute result
			str - 	if Exception: returs error message
		'''
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
				logger.warning(query)
				if insert_val:
					logger.debug("Insert val: {vals}".format(vals=insert_val))
				sys.stderr.write(str(e)+"\n")
			return(e)

	def insert(self,data,table):
		'''Insert function
				data is a dictionary with keys matching
				table columns with respective value to be inserted
		------
		Returns
			see query responses
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
		'''Update function requires table column which column to identify row with and value to replace

		------
		Returns
			boolean - 	True: if rowcount of update is not zero
						False if column is not set (value is inserted)
		'''
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

	def multi_update(self,data,table):
		'''Update function requires table column which column to identify row with and value to replace
		Parameters
			list - data in form of a list with updates
			table - which table to update
		------
		Returns
			boolean - 	True: if rowcount of update is not zero
						False if column is not set (value is inserted)
		'''
		UPDATE_QUERY = '''
			INSERT INTO {table} ({set_column},{where_column})
			VALUES {data}
			ON CONFLICT({where_column}) DO UPDATE SET
			{where_column} = {where_column},
			{set_column}   = {set_column};
		'''.format(table=table,set_column=data["set_column"],where_column=data["where_column"],data=",".join(data["data"]))
		#udata = tuple([data["set_value"],data["where"]])
		logger.debug("{q}".format(q=UPDATE_QUERY))
		#exit()
		res = self.query(UPDATE_QUERY,error=True)
		#if self.rowcount() == 0:
		#	res = self.insert({data["set_column"]:data["set_value"],data["where_column"]:data["where"]}, table)
		#	return False
		return True

	def delete(self,nodes,table):
		'''Deleting a node should make sure all related data is also removed

		------
		Returns
			see query responses
		'''
		DELETE_QUERY = '''
				DELETE FROM {table}
					WHERE id in({values})
		'''.format(table=table,
					values=','.join(["?" for x in nodes])
		)
		return self.query(DELETE_QUERY,insert_val=values)

	def rowcount(self):
		'''Get the number row number the cursor is currently at

		------
		Returns
			int - cursor rowcount'''
		return self.cursor.rowcount

class DatabaseFunctions(DatabaseConnection):
	"""DatabaseFunctions class defines additional functions to the DatabaseConnection class

	"""
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

		------
		Returns
			True: if all edges has one and only one parent
		'''
		self.taxonomy = self.get_nodes()

		logger.info("Get all database nodes")
		self.nodes = self.get_nodes(col=1)					## get all nodes
		logger.info("Get all database edges")
		self.edges = self.get_links()					## Get all links
		logger.info("Get all children from root node")
		self.tree_connections = self.get_children([1])  ## Get all chilren from root
		logger.info("Get tree edges from children")
		self.tree_links = self.get_links(self.tree_connections) ## get all links from rooted tree
		logger.info("Get nodes from tree edges")
		self.tree_nodes = set()
		for edge in self.edges:
			self.tree_nodes.add(edge[0])
			self.tree_nodes.add(edge[1])
		logger.info("Validate parents")
		res = self.check_parent()
		extra = ""
		if len(res) > 0:
			failed_nodes = [self.nodes[x[0]] for x in list(res)]
			logger.debug(failed_nodes)
			if len(failed_nodes) > 10:
				logger.info(failed_nodes[:10])
				failed_nodes = len(failed_nodes)
				extra = "total failed nodes: "
			raise TreeError("There are nodes with multiple parents {extra}{nodes}".format(nodes=failed_nodes,extra=extra))
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
			lset = nodeset - self.tree_connections
			if len(lset) > 10:
				lset = len(lset)
				logger.debug([self.taxonomy[x] for x in lset])
			logger.info("{}".format(lset))
			raise TreeError("The number of nodes and the number nodes under root does not match!")
		if len(set(self.edges) - set(self.tree_links)) != 0:
			lset = set(self.edges) - set(self.tree_links)
			if len(lset) > 10:
				lset = len(lset)
				logger.debug([self.taxonomy[x] for x in lset])
			logger.info("{}".format(lset))
			raise TreeError("The number of edges and the number edges under root does not match!")
		if len(self.tree_nodes - nodeset) != 0:
			lset = self.tree_nodes - nodeset
			if len(lset) > 10:
				lset = len(lset)
				logger.debug([self.taxonomy[x] for x in lset])
			logger.info("{}".format(lset))
			raise TreeError("The number of annotated nodes does not match with the number of nodes connected to edges!")
		if len(self.tree_connections - nodeset) != 0:
			lset = self.tree_connections - nodeset
			if len(lset) > 10:
				lset = len(lset)
				logger.debug([self.taxonomy[x] for x in lset])
			logger.info("{}".format(lset))
			raise TreeError("There are nodes in the database not connected to the tree")
		logger.info("Validation OK!")
		return True

	def statistics(self):
		'''Print statistics of the database

		------
		Returns:
			true
		'''
		logger.info("Get all database nodes")
		self.nodes = self.get_nodes(col=1)					## get all nodes
		logger.info("Get all database edges")
		self.edges = self.get_links()					## Get all links
		logger.info("Get annotated genomes")
		self.genomes = self.get_genomes()

		stats = """Tree statistics
					Nodes: {nodes}
					Links: {links}
					Genomes: {ngenomes}
					""".format(nodes = len(self.nodes),
								links=len(self.edges),
								ngenomes=len(self.genomes)
		)
		print(stats)
		return True

	def check_parent(self):
		'''check tree structure parent
			make sure no tree node has multiple parents
		------
		Returns
			True: if node has one and only one parent
		'''
		try:
			QUERY = "SELECT child FROM tree GROUP BY child HAVING count(parent) > 1"  ## Thanks to andrewjmc@github for this suggestion
			logger.debug(QUERY)
			child_w_dp = self.query(QUERY).fetchall()
			QUERY = "SELECT parent,rank_i FROM tree WHERE child in ({children})".format(children=",".join(map(str,list(*child_w_dp))))
			p_ranks = self.query(QUERY).fetchall()
			if len(p_ranks) != len(set(p_ranks)):
				logger.error("Nodes with two parents have identical ranks, Fatal Error")
			elif len(child_w_dp) > 0:
				logger.info("Found: {n} nodes with multiple parents, however in different lineages, OK.".format(n=len(child_w_dp)))
				return []
			return child_w_dp
		except AttributeError:
			logger.info("AttributeError occured")
			logger.info(QUERY)
			logger.info(child_w_dp)
		raise TreeError("Node: {node} has more than one parent!".format(node=name))

	'''Get functions of class'''
	def get_all(self, database=False, table=False):
		'''Get full table from table

		------
		Returns
			list - results from table
		'''
		if not table:
			raise Exception("Table was not supplied to DatabaseFuction get_all()")
		QUERY = '''SELECT * FROM {table}'''.format(table)
		return self.query(QUERY).fetchall()

	def get_taxid_base(self):
		'''Fetch the next incremental node from the current database

		------
		Returns
			int - returns highest index number from nodes table
		'''
		QUERY = "SELECT MAX(id) AS max_id FROM nodes"
		res = self.query(QUERY).fetchone()[0]
		if res == None:
			res = 1
		else:
			res += 1
		return res

	def get_genomes(self, database=False,limit=0,table="genomes",cols="id,genome"):
		'''Get the list of genomes in the database

		------
		Returns
			dict - dictionary with id to genome (or selected cols)
		'''
		## This is a many to many relation, so all genomes has to be put in a set for each taxonomy id
		genomeDict = {}
		QUERY = '''SELECT {cols} FROM {table}'''.format(cols=cols,table=table)
		if limit > 0:
			QUERY += " LIMIT {limit}".format(limit=limit)
		for id,genome in self.query(QUERY).fetchall():
			genomeDict[genome] = id
		return genomeDict

	def get_node(self, name, database=False):
		'''Retrieve one node by name

		------
		Returns
			list - list with node id
		'''
		QUERY = '''SELECT id FROM nodes where name = "{name}"'''.format(name=name)
		if not database:
			database = self.database
		return [n[0] for n in self.query(QUERY).fetchall()]

	def get_nodes(self, database=False,col=False):
		'''Retrieve the whole node info table of the database to decrease the number of database calls!

		------
		Returns
			dict - dictionary with node id to name and name to node id
		'''
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

	def get_links(self, nodes=False,database=False,swap=False,only_parents=False,order=False,simple=False):
		'''This function returns all links in the given database

		------
		Returns
			list - list of child to parent (or vice versa) with rank_i
		'''
		order = ["parent","child"]
		if not simple:
			order += ["rank_i"]
		logger.debug(nodes)
		if swap:
			order[1],order[0] = order[0],order[1]

		if only_parents and nodes:
			QUERY = '''SELECT {order} FROM tree WHERE child in ({nodes})'''.format(nodes=",".join(map(str,nodes)),order=",".join(order))
		elif nodes and order:
			QUERY = '''SELECT {order} FROM tree WHERE child in ({nodes}) ORDER BY parent ASC'''.format(nodes=",".join(map(str,nodes)),order=",".join(order))
		elif nodes:
			QUERY = '''SELECT {order} FROM tree WHERE parent in ({nodes}) OR child in ({nodes})'''.format(nodes=",".join(map(str,nodes)),order=",".join(order))
		else:
			QUERY = '''SELECT {order} FROM tree'''.format(order=",".join(order))
		logger.debug(QUERY)
		if not database:
			database = self.database
		links = self.query(QUERY).fetchall()
		return links

	'''Add functions of class'''
	def add_node(self, description, id=False, table="nodes"):
		'''Add node to tree

		Returns
		------
		int - lastrowid/taxid_base
		'''
		if description.strip() == "":  ## for some reason an empty node has been given, do not add empty nodes
			return False
		info = { "name": description }
		### If ID is supplied skip autoincrement and add specific ID
		if id:
			info["id"] = id
		taxid_base = self.insert(info, table=table)
		logger.debug("node added [description {}, taxid base {}]: ".format(info, taxid_base))
		return taxid_base

	def add_rank(self, rank ,id=False):
		'''Add node to tree

		Returns
		------
		int - id of added rank level
		'''
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
		'''Add relationship in tree

		Returns
		------
			see insert responses
		'''
		info = {
			"child": child,
			"parent": parent,
			"rank_i": rank
		}
		logger.debug("link added:  child {}, parent {}, rank {} ".format(child,parent,rank))
		return self.insert(info, table="tree")

	def add_genome(self, genome, _id=False):
		'''Add genome annotation to nodes

		Returns
		------
			see insert responses
		'''
		info = {
			"genome": genome
		}
		if _id:
			info["id"] = _id
		return self.insert(info, table="genomes")

	def add_links(self,links, table="tree",hold=False):
		'''Add links from a list to tree

		Returns
		------
			list 	- tree links added
			set		- set of unique nodes in added links
		'''
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
		'''Add nodes from a list of nodes

		Returns
		------
			int - number of nodes added ()
		'''
		added_nodes = 0
		for node in nodes:
			res = self.add_node(node)
			### If the node already exist in the database
			if "UNIQUE constraint failed" not in str(res):
				added_nodes +=1
		## Commit changes
		if not hold:
			self.commit()
		return added_nodes

	'''Delete functions of class'''
	def delete_links(self,links, table="tree",hold=False):
		'''This function deletes all links given in links

		Returns
		------
			boolean
		'''
		logger.info("Slow clean")
		QUERY = "DELETE FROM {table} WHERE parent = {parent} AND child = {child} AND rank_i = {rank_i}"
		logger.debug("Deleting {nlinks} links!".format(nlinks=len(links)))
		logger.debug(QUERY.format(table=table,parent="",child="",rank_i=""))
		logger.info("{links}".format(links=len(links)))
		for parent,child,rank in links:
			logger.debug("{}-{} rank: {} deleted!".format(parent,child,rank))
			res = self.query(QUERY.format(table=table, parent=parent, child=child, rank_i=rank))
		## Commit changes
		if not hold:
			logger.debug("Commit changes!")
			self.commit()
		return True

	def fast_delete_links(self,links,table="tree",hold=False):
		'''This function is used when general clean function is executed
		Returns
		------
			boolean
		'''
		logger.info("Fast clean")
		parents = set()
		children = set()
		QUERY = "DELETE FROM {table} WHERE parent in ({parents}) AND child in ({children})"
		logger.debug("Deleting {nlinks} links!".format(nlinks=len(links)))
		for parent,child,rank in links:
			parents |= set([parent])
			children |= set([child])
			#logger.debug("{}-{} rank: {} deleted!".format(parent,child,rank))
		logger.debug(QUERY.format(table=table,parents=len(parents),children=len(children)))
		res = self.query(QUERY.format(table=table, parents=",".join(list(map(str,parents))), children=",".join(list(map(str,children)))))
		## Commit changes
		if not hold:
			logger.debug("Commit changes!")
			self.commit()
		return True

	def delete_nodes(self, nodes, table="nodes",hold=False):
		'''This function deletes all nodes given in nodes

		Returns
		------
			boolean
		'''
		QUERY = "DELETE FROM {table} WHERE id in ({nodes})"
		logger.debug("Deleting {nnodes} nodes!".format(nnodes=len(nodes)))
		logger.debug(QUERY.format(table=table,nodes=""))
		res = self.query(QUERY.format(table=table, nodes=",".join(list(map(str,nodes)))))
		## Commit changes
		if not hold:
			logger.debug("Commit changes!")
			self.commit()
		return True

	def delete_genomes(self,nodes,table="genomes",hold=False):
		'''This function clean up genomes from given nodes

		Returns
		------
			boolean
		'''
		QUERY = "DELETE FROM {table} WHERE id in({nodes})"
		logger.info("Deleting {nnodes} annotations!".format(nnodes=len(nodes)))
		logger.debug(QUERY.format(table=table,nodes=""))
		#for node in nodes:
			#logger.info("Delete genomes from: {node}".format(node=node))
		res = self.query(QUERY.format(table=table, nodes=",".join(list(map(str,nodes)))))
		## Commit changes
		if not hold:
			logger.debug("Commit changes!")
			self.commit()
		return True

	def num_rows(self,table):
		'''Find out the number of rows in a table

		Returns
		------
			int - the number of rows in a table'''
		QUERY = '''SELECT Count(*) FROM {table}'''
		logger.debug(QUERY)
		return self.query(QUERY.format(table=table)).fetchall()[0][0]


class ModifyFunctions(DatabaseFunctions):
	"""ModifyFunctions adds nessesary functions when modifying a database"""
	def __init__(self, database, verbose=False):
		super().__init__(database, verbose)
		logger.debug("Load ModifyFunctions")

	def get_rank(self,col=1):
		'''Get rank index from database

		Returns
		------
			dict - rank to rank_index and vice versa
		'''
		QUERY = "SELECT rank_i,rank FROM rank"
		logger.debug(QUERY)
		rankDict = {}
		for rank in self.query(QUERY).fetchall():
			rankDict[rank[0]] = rank[1]
			if col == 1:
				continue
			rankDict[rank[1]] = rank[0]
		return rankDict

	def get_children(self,parents,children=set(),level=0,maxdepth=50,selected=False):
		'''Get all children from a parent

		Returns
		------
			set - unique list of children from a decending tree
		'''
		if level >= maxdepth:
			return children
		QUERY = '''SELECT child FROM tree WHERE parent in({nodes})'''.format(nodes=",".join(map(str,list(parents))))
		if selected:
			QUERY = '''SELECT child FROM tree WHERE parent in({nodes}) and rank_i = {rank}'''.format(nodes=",".join(map(str,list(parents))), rank=selected)
		#logger.debug(QUERY)
		#print(QUERY)
		res = self.query(QUERY).fetchall()
		if (len(res) + len(children)) != 0:
			children = set([child_i[0] for child_i in res])
			if level < maxdepth:
				children |= self.get_children(parents=children,level=level+1,maxdepth=maxdepth)
		return children

	def get_parent(self,name,all=False):
		'''Get parent from node id parent

		Returns
		------
			list - parent link and rank
		'''
		#QUERY = '''SELECT parent,child,rank FROM tree LEFT JOIN rank on (tree.rank_i = rank.rank_i) WHERE child = "{node}"'''.format(node=name)
		QUERY = '''SELECT parent,child,rank_i FROM tree WHERE child = "{node}"'''.format(node=name)
		logger.debug(QUERY)
		if all:
			res = self.query(QUERY).fetchall()
		else:
			res = self.query(QUERY).fetchone()
		return res

	def parse_parents(self,parents):
		pret = set()
		for pc in parents:
			pret |= set([int(pc[0])])
			pret |= set([int(pc[1])])
		return pret

	def get_parents(self,name,parents=set(),depth=0,find_all=False,simple=False,ncbi=False):
		'''Get all parents until root

		Returns
		------
			list - all parents of a node
		'''
		ld = depth+1
		lmax = 12
		if ncbi and depth == 0:
			logger.info("NCBI mode on, minimum level to check 50 (since euk has red many levels) Will take time!")
		lmax = 50
		if simple:
			simple=""
		else:
			simple=",rank_i"
		try:
			name = set(list(map(int,name))) ## Make sure all names are int
		except TypeError:
			name = set(list(map(int,[x[0] for x in name]))) ## Make sure all names are int
		if isinstance(name, int):
			name = [name]
		if len(name & set([1])) > 0 and not ncbi:
			return parents
		QUERY = '''SELECT parent,child{simple} FROM tree WHERE child in ({node})'''.format(simple=simple, node=",".join(map(str,name)))
		#if find_all: logger.info(QUERY)
		res = self.query(QUERY).fetchall()
		if not find_all:
			try:
				res = int(res[0])
			except IndexError:
				logger.warning("WARNING: parent could not be found for node {res} \n{query}".format(query=QUERY, res=name))
				return set()
			except TypeError:
				#print(res)
				logger.error("ERROR: cannot fetch parents when parents are duplicated!")
				raise TreeError("cannot fetch parents when parents are duplicated!")
			if res[0] == res[1]:  ## Special for root
				return set([res[0]])
			else:
				## Add all parents
				try:
					#if ncbi:
					#	logger.info("Loop depth: {ld}".format(ld=ld))
					parents |= self.get_parents([int(res[0])],depth=ld,ncbi=ncbi)

				except RecursionError:
					print(res, name)
			## Add current node
			parents |= set([int(res[0])])
		else:
			parents = self.parse_parents(res) #set([(int(pc[0]),int(pc[-1])) for pc in res])
			if len(parents-name) == 1 and ld > lmax:  ## There are no more parents to fetch
				return parents
			else:
				## Add all parents
				if ld > lmax:
					return parents
				try:
					#if ncbi:
					#	logger.info("Loop depth: {ld}".format(ld=ld))
					parents |= self.get_parents(parents,find_all=find_all,depth=ld,ncbi=ncbi)
				except RecursionError:
					logger.info("RecursionError, node requested")
					print(res, name)
		return parents

	def get_id(self,name):
		'''get node id from name

		Returns
		------
			int - node id from node name
		'''
		QUERY = '''SELECT id FROM nodes WHERE name = "{node}" COLLATE NOCASE'''.format(node=name)
		try:
			res = self.query(QUERY).fetchone()[0]
		except TypeError:
			logger.debug(QUERY)
			raise NameError("Name not found in the database! {name}".format(name=name))
		return res

	def update_genome(self,data):
		'''Add genome annotation to nodes

		Returns
		------
			see update responses
		'''
		return self.update(data, table="genomes")
