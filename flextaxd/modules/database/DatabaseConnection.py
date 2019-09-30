import sys
import os
import sqlite3

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

class DatabaseConnection(object):
	"""docstring for DatabaseConnection"""
	def __init__(self, database, verbose=False):
		super().__init__()
		self.verbose = verbose
		self.database = database
		BASE_DIR = os.path.dirname(os.path.abspath(__file__))  ## Retrieve path
		if not os.path.exists(self.database):
			if self.verbose:
				print("python CreateDatabase.py {database}".format(database=self.database))
			os.system("python {path}/CreateDatabase.py {database}".format(path=BASE_DIR,database=self.database))
		self.conn = self.connect(self.database)
		self.cursor = self.create_cursor(self.conn)

	def set_verbose(self,val):
		self.verbose = val

	def connect(self,database):
		'''Create database connection'''
		try:
			conn = sqlite3.connect(database)
			if self.verbose:
				print("{database} opened successfully.".format(database=database))

			return conn
		except Exception as e:
			sys.stderr.write(str(e))
		raise ConnectionError("Count not connect to the database {database} see above message for details!".format(database=database))
		return None

	def create_cursor(self,conn):
		'''Create a db cursor'''
		try:
			cursor = conn.cursor()
			if self.verbose:
				print("cursor created.")
			return cursor
		except Exception as e:
			sys.stderr.write(str(e))
		return None

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
				print("Error in DatabaseConnection query")
				if self.verbose: print(query)
				#if self.verbose: print("Insert val", insert_val)
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
		if self.verbose: print(UPDATE_QUERY,udata)
		#print(UPDATE_QUERY, udata)
		res = self.query(UPDATE_QUERY,udata,error=True)
		'''If the genome does not exist, it should instead be inserted in the database'''

		if self.rowcount() == 0:
			res = self.insert({data["set_column"]:data["set_value"],data["where_column"]:data["where"]}, table)
			return False
		return True

	# def delete(self,nodes,table):
	# 	'''Deleting a node should make sure all related data is also removed'''
	# 	DELETE_QUERY = '''
	# 			DELETE FROM {table}
	# 				WHERE id in({values})
	# 	'''.format(table=table,
	# 				values=','.join(["?" for x in nodes])
	# 	)
	# 	return self.query(DELETE_QUERY,insert_val=values)

	def rowcount(self):
		return self.cursor.rowcount

class DatabaseFunctions(DatabaseConnection):
	"""docstring for DatabaseFunctions."""
	def __init__(self, database, verbose=False):
		super().__init__(database, verbose)
		if self.verbose: print("Load DatabaseFunctions")

	'''Get functions of class'''
	def get_taxid_base(self):
		'''Fetch the next incremental node from the current database'''
		QUERY = "SELECT MAX(id) AS max_id FROM nodes"
		res = self.query(QUERY).fetchone()[0]
		if res == None:
			res = 1
		else:
			res += 1
		return res

	def get_genomes(self, database=False,limit=0):
		'''Get the list of genomes in the database'''
		## This is a many to many relation, so all genomes has to be put in a set for each taxonomy id
		genomeDict = {}
		QUERY = '''SELECT id,genome FROM {table}'''.format(table="genomes")
		if limit > 0:
			QUERY += " LIMIT {limit}".format(limit=limit)
		for id,genome in self.query(QUERY).fetchall():
			genomeDict[genome] = id
		return genomeDict

	def get_all(self, database=False, table=False):
		'''Get full table from table'''
		QUERY = '''SELECT * FROM {table}'''.format(table)
		return self.query(QUERY).fetchall()

	def get_nodes(self, database=False,col=False):
		'''Retrieve the whole node info table of the database to decrease the number of database calls!'''
		nodeDict = {}
		QUERY = '''SELECT id,name FROM nodes'''
		if not database:
			database = self.database
		for node in self.query(QUERY).fetchall():
			nodeDict[node[0]] = node[1]
			if col == 1:
				continue
			nodeDict[node[1]] = node[0]
		return nodeDict

	def get_links(self, nodes,database=False,swap=False):
		'''This function returns all links in the given database'''
		order = ["parent","child"]
		if swap:
			order[1],order[0] = order[0],order[1]
		QUERY = '''SELECT {order},rank_i FROM tree WHERE parent in ({nodes}) OR child in ({nodes})'''.format(nodes=",".join(map(str,nodes)),order=",".join(order))
		if not database:
			database = self.database
		links = self.query(QUERY).fetchall()
		return links

	'''Add functions of class'''
	def add_node(self, description, id=False):
		'''Add node to tree'''
		if description.strip() == "":  ## for some reason an empty node has been given, do not add empty nodes
			return False
		info = { "name": description }
		### If ID is supplied skip autoincrement and add specific ID
		if id:
			info["id"] = id
		taxid_base = self.insert(info, table="nodes")
		return taxid_base

	def add_rank(self, rank,id=False):
		'''Add node to tree'''
		info = { "rank": rank }
		### If ID is supplied skip autoincrement and add specific ID
		if id:
			info["id"] = id
		rank_id = self.insert(info, table="rank")
		return rank_id

	def add_link(self, child, parent, rank=1, table="tree"):
		'''Add relationship in tree'''
		info = {
			"child": child,
			"parent": parent,
			"rank_i": rank
		}
		if self.verbose: print("link added: ",child,parent,rank)
		return self.insert(info, table="tree")

	def add_genome(self, id, genome):
		'''Add genome annotation to nodes'''
		info = {
			"id": id,
			"genome": genome
		}
		return self.insert(info, table="genomes")

	def add_links(self,links, table="tree",hold=False):
		'''Add links from a list to tree'''
		added_links = 0
		nodes = set()
		for parent,child,rank in links:
			#print(parent,child,rank)
			res = self.add_link(child,parent,rank,table=table)
			### Check if the link already exist in the database, this overlap may occur when a large new branch is added
			if "UNIQUE constraint failed" not in str(res):
				added_links +=1
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
		if self.verbose: print(QUERY.format(table=table,parent="",child=""), links)
		for parent,child,rank in links:
			res = self.query(QUERY.format(table=table, parent=parent, child=child))

		## Commit changes
		if not hold:
			self.commit()

	def delete_nodes(self, nodes, table="nodes",hold=False):
		'''This function deletes all nodes given in nodes'''
		QUERY = "DELETE FROM {table} WHERE id = {node}"
		if self.verbose: print(QUERY.format(table=table,node=""), nodes)
		for node in nodes:
			res = self.query(QUERY.format(table=table, node=node))
		## Commit changes
		if not hold:
			self.commit()

	def num_rows(self,table):
		'''Return the number of rows in a table'''
		QUERY = '''SELECT Count(*) FROM {table}'''
		return self.query(QUERY.format(table=table)).fetchall()[0][0]


class ModifyFunctions(DatabaseFunctions):
	"""ModifyFunctions adds a few important functions only nessesary when modifying a database"""
	def __init__(self, database, verbose=False):
		super().__init__(database, verbose)
		if self.verbose: print("Load ModifyFunctions")

	def get_rank(self,col=1):
		'''Get rank index from database'''
		QUERY = "SELECT rank_i,rank FROM rank"
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
		res = self.query(QUERY).fetchall()
		if (len(res) + len(children)) != 0 and level < maxdepth:
			children = set([child_i[0] for child_i in res])
			children |= self.get_children(parents=children,level=level+1)
		return children

	def get_parent(self,name):
		'''return parent'''
		#QUERY = '''SELECT parent,child,rank FROM tree LEFT JOIN rank on (tree.rank_i = rank.rank_i) WHERE child = "{node}"'''.format(node=name)
		QUERY = '''SELECT parent,child,rank_i FROM tree WHERE child = "{node}"'''.format(node=name)
		res = self.query(QUERY).fetchone()
		return res

	def get_id(self,name):
		'''get node id from name'''
		QUERY = '''SELECT id FROM nodes WHERE name = "{node}"'''.format(node=name)
		try:
			res = self.query(QUERY).fetchone()[0]
		except TypeError:
			print(QUERY)
			raise NameError("Name not found in the database! {name}".format(name=name))
		return res


	def update_genome(self,data):
		'''Add genome annotation to nodes'''
		return self.update(data, table="genomes")
