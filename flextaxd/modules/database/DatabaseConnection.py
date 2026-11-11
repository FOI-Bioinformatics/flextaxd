import sys
import os
import sqlite3
import subprocess
import datetime
import logging
logger = logging.getLogger(__name__)

MAX_BIND_PARAMS = 900 ## Maximum number of "?" placeholders per query, keep below the lowest
                      ## compile time limit of sqlite (999), the limit varies between builds

def flextaxd_version():
	'''The version of flextaxd that is running, recorded in the meta table of every database.
			Imported inside the function: the modules package is also used as a top level package,
			where the installed flextaxd package may not be reachable'''
	try:
		from flextaxd import __version__
		return __version__
	except ImportError:
		return "unknown"

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
		self.meta_ready = False
		BASE_DIR = os.path.dirname(os.path.abspath(__file__))  ## Retrieve path
		created = not os.path.exists(self.database)
		if created:
			cmd = [sys.executable, os.path.join(BASE_DIR,"CreateDatabase.py"), self.database]
			if self.verbose:
				logger.debug(" ".join(cmd))
			if subprocess.call(cmd) != 0:
				raise ConnectionError("Could not create the database {database}, see message above for details!".format(database=self.database))
		try: ## If database connection already exists
			self.conn
		except AttributeError:
			logger.debug("Connecting to {database}".format(database=self.database))
			self.conn = self.connect(self.database)
			self.cursor = self.create_cursor(self.conn)
		## Databases written by older flextaxd versions have no meta table
		self.ensure_meta_table()
		## and no unique index covering the links that carry no rank
		self.ensure_tree_link_index()
		if created and self.get_meta("built_with") is None:
			self.set_meta("built_with",flextaxd_version())
			self.set_meta("built",datetime.datetime.now().isoformat(timespec="seconds"))

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
		raise ConnectionError("Count not create db cursor to {database}".format(database=self.database))

	def ensure_meta_table(self):
		'''Create the meta and report tables of databases written before flextaxd recorded build
				information, does nothing to databases that already have them

		------
		Returns:
			true
		'''
		self.query("""CREATE TABLE IF NOT EXISTS meta (
						key VARCHAR(255) PRIMARY KEY,
						value TEXT
					);""")
		self.query("""CREATE TABLE IF NOT EXISTS report (
						kind VARCHAR(31) NOT NULL,
						source VARCHAR(255) NOT NULL,
						entry TEXT NOT NULL,
						UNIQUE (kind,source,entry)
					);""")
		self.conn.commit()
		self.meta_ready = True
		return True

	def ensure_tree_link_index(self):
		'''Give a database written before the link index existed the unique index over the links of the
				tree, first removing the links that an older version wrote twice. The unique constraint of the
				tree table cannot cover a link that carries no rank: NULL never equals NULL in a sqlite unique
				constraint, so importing a taxonomy into a database that already held those links added the
				identical rank-less link again on every run (issue #147). Databases created now get the index
				from CreateDatabase and are marked there, so this runs once per older database.

		------
		Returns:
			boolean - whether the index is in place
		'''
		if self.get_meta("tree_link_index") is not None:
			return True
		## a link is the same link when parent, child and rank match, with no rank standing for -1, the
		## same stand-in the delete of a link uses
		find_copies = ("rowid NOT IN (SELECT MIN(rowid) FROM tree"
						" GROUP BY parent,child,ifnull(rank_i,-1))")
		copies = self.select("SELECT COUNT(*) FROM tree WHERE {where};".format(where=find_copies))
		if not copies:
			## no tree table at all: a database that nothing has been imported into yet
			return True
		num_duplicates = int(copies[0][0])
		if num_duplicates:
			self.query("DELETE FROM tree WHERE {where};".format(where=find_copies))
			logger.warning("{num} tree links of {database} were written twice by an older flextaxd,"
							" they are removed now that the link is unique (issue #147)".format(
								num=num_duplicates, database=self.database))
		res = self.query("""CREATE UNIQUE INDEX IF NOT EXISTS tree_link
							ON tree (parent,child,ifnull(rank_i,-1));""")
		if isinstance(res, Exception):
			logger.warning("Could not add the unique tree link index to {database}: {err}".format(
				database=self.database, err=str(res)))
			return False
		self.set_meta("tree_link_index","1")
		return True

	def set_meta(self,key,value):
		'''Store a single value in the meta table, replacing any previous value of that key

		------
		Returns:
			str - the value stored
		'''
		self.query("INSERT OR REPLACE INTO meta (key,value) VALUES (?,?);",insert_val=(key,value))
		self.conn.commit()
		return value

	def get_meta(self,key):
		'''Read a single value from the meta table

		------
		Returns:
			str - the value of the key, None when the key is not set
		'''
		res = self.select("SELECT value FROM meta WHERE key = ?;",params=(key,))
		if res:
			return res[0][0]
		return None

	def get_build_info(self):
		'''The flextaxd build information stored in the meta table

		------
		Returns:
			dict - built_with, built, modified_with and modified keys, values None when unset
		'''
		return {info: self.get_meta(info) for info in ["built_with","built","modified_with","modified"]}

	def save_report(self,kind,source,count,detail=False):
		'''Store what one step found for one input as a number only: the individual genomes and rows
				are never kept in the database, they are either written to a file by the step (the
				"<input>.skipped" files) or can be recomputed from the genome folder with --genomes_path.
				Findings of the same kind from the same input of an earlier run are replaced, a count of
				zero removes the stored finding.

		------
		Parameters:
			str - kind of finding, missing_on_disk, unannotated_file, unannotated_sequence or skipped_import
			str - the input the findings came from, a genome path or an input file
			int - how many were found
			str - where the individual entries can be read, when the step wrote them to a file
		------
		Returns:
			int - the number stored
		'''
		count = int(count)
		self.query("DELETE FROM report WHERE kind = ? AND source = ?;",insert_val=(kind,source))
		if count:
			self.query("INSERT INTO report (kind,source,entry) VALUES (?,?,?);",
						insert_val=(kind,source,self.report_count_entry(count,detail)))
		self.conn.commit()
		return count

	def report_count_entry(self,count,detail=False):
		'''The single entry that holds a finding count, marked so that report_lines can tell it apart
				from the individual entries that databases written before this change contain'''
		if detail:
			return "total\t{count}\t{detail}".format(count=count,detail=detail)
		return "total\t{count}".format(count=count)

	def report_findings(self,kind):
		'''The findings of one kind as they are stored: one entry with the number per input

		------
		Returns:
			tuple - (total, list of "input: number, listed in file" lines, list of individual
					entries of databases written by an earlier flextaxd version)
		'''
		total = 0
		lines = []
		legacy = []
		for source,entry in self.get_report(kind):
			if entry.startswith("total\t"):
				fields = entry.split("\t")
				number = int(fields[1])
				total += number
				line = "{source}: {number}".format(source=source,number=number)
				if len(fields) > 2:
					line += ", listed in {file}".format(file=fields[2])
				lines.append(line)
			else:
				total += 1
				legacy.append("{source}\t{entry}".format(source=source,entry=entry))
		return total, lines, legacy

	def get_report(self,kind,source=False):
		'''Read the stored findings of one kind

		------
		Returns:
			list - tuples of (source, entry)
		'''
		if source:
			return self.select("SELECT source,entry FROM report WHERE kind = ? AND source = ? ORDER BY entry;",params=(kind,source))
		return self.select("SELECT source,entry FROM report WHERE kind = ? ORDER BY entry;",params=(kind,))

	def report_kinds(self):
		'''How many findings of each kind the database holds

		------
		Returns:
			dict - kind to number of entries
		'''
		return {row[0]: row[1] for row in self.select("SELECT kind,COUNT(*) FROM report GROUP BY kind;")}

	def commit(self):
		'''Commit the transaction and record which flextaxd version changed the database last'''
		self.conn.commit()
		if getattr(self,"meta_ready",False):
			self.set_meta("modified_with",flextaxd_version())
			self.set_meta("modified",datetime.datetime.now().isoformat(timespec="seconds"))

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

	def select(self,query,params=False):
		'''Run a SELECT through query and return the rows, an empty list if the query failed
				(query returns the exception object instead of a cursor on failure)

		Returns
			list - rows returned by the query
		'''
		if params:
			res = self.query(query,insert_val=tuple(params),error=True)
		else:
			res = self.query(query)
		if isinstance(res,Exception):
			return []
		return res.fetchall()

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
			optional key "columns" in data - every column of the value rows in the order they
				are written, default the two update columns (set_column, where_column). Rows that
				already exist get every column but the identifying one from the incoming value,
				an incoming NULL keeps the value that is already stored.
		------
		Returns
			boolean - 	True: the statement ran
						False: the statement failed (the rows were not written)
		'''
		columns = data.get("columns",[data["set_column"],data["where_column"]])
		UPDATE_QUERY = '''
			INSERT INTO {table} ({columns})
			VALUES {data}
			ON CONFLICT({where_column}) DO UPDATE SET
			{assignments};
		'''.format(table=table,
					columns=",".join(columns),
					data=",".join(data["data"]),
					where_column=data["where_column"],
					assignments=", ".join("{col} = COALESCE(excluded.{col}, {table}.{col})".format(col=col,table=table) for col in columns if col != data["where_column"]))
		#udata = tuple([data["set_value"],data["where"]])
		logger.debug("{q}".format(q=UPDATE_QUERY))
		#exit()
		res = self.query(UPDATE_QUERY,error=True)
		if isinstance(res,Exception):
			logger.error("The batched update of the {table} table failed: {err}".format(table=table,err=str(res)))
			return False
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
		return self.query(DELETE_QUERY,insert_val=tuple(nodes))

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
	def _table_count(self, statement):
		'''A count run by the validation: 0 when the tables are not there yet, the way validation of a
				database that nothing has been imported into used to pass on empty lists.

		------
		Returns:
			int
		'''
		rows = self.select(statement)
		if not rows:
			return 0
		return int(rows[0][0])

	def _node_names(self, node_ids):
		'''The names of a handful of node ids for a message, the id itself for an id the database has no
				node row for (which is one of the things validation reports).

		------
		Returns:
			list of str
		'''
		ids = [int(node_id) for node_id in node_ids]
		if not ids:
			return []
		names = {int(node_id): name for node_id, name in self.select(
			"SELECT id,name FROM nodes WHERE id IN ({})".format(",".join(["?"]*len(ids))), params=tuple(ids))}
		return [str(names.get(node_id, node_id)) for node_id in ids]

	def _root_walk_table(self, name="ftx_temp_under_root"):
		'''Walk the tree once into a temporary table: the ids that can be reached as children from the root
				node, which is the set the python walk of the children of node 1 built before. The checks of
				the validation all look at that set, so it is walked once and holds its id as primary key for
				the checks to probe; walking it again for every check cost a quarter of a minute each on the
				2.7 million links of a taxdump the size of the NCBI one (issue #107). sqlite ends a recursive
				with clause at the first row it has seen before, so a database that holds a link back to an
				ancestor stops here instead of running on for ever. A database whose tree table is not there
				yet ends with an empty set, which is what the validation read from an empty tree before.

		------
		Returns:
			str - the name of the temporary table
		'''
		self.query("DROP TABLE IF EXISTS {name}".format(name=name))
		self.query("CREATE TEMP TABLE {name} (id integer PRIMARY KEY)".format(name=name))
		res = self.query(
			"INSERT OR IGNORE INTO {name} (id)"
			" WITH RECURSIVE under_root(id) AS ("
			" SELECT tree.child FROM tree WHERE tree.parent = 1"
			" UNION"
			" SELECT tree.child FROM tree JOIN under_root ON tree.parent = under_root.id)"
			" SELECT id FROM under_root;".format(name=name))
		if isinstance(res, Exception):
			logger.debug("The tree under root could not be walked: {err}".format(err=str(res)))
		return name

	def _validate_ids(self, message, id_select):
		'''One check of the validation: how many ids the select finds, and a few of them named in the
				message. Raises the TreeError of that check, the wording the set comparison in python used
				to give.

		------
		Returns:
			boolean - False when nothing was found
		'''
		num = self._table_count("SELECT COUNT(*) FROM ({ids});".format(ids=id_select))
		if not num:
			return False
		examples = self.select("SELECT id FROM ({ids}) LIMIT 10;".format(ids=id_select))
		found = ", ".join(self._node_names([row[0] for row in examples]))
		logger.info("{num} ids fail here, for example: {found}".format(num=num, found=found))
		raise TreeError("{message} total: {num}".format(message=message, num=num))

	def validate_tree(self):
		'''This function validates the tree structure in the databases
			1. All nodes must only have one parent (check_parent)
			2. All edges must be attatched to the tree (match tree links and edges)
			3. All nodes must be attatched to the tree (match tree nodes and nodes)
			4. All nodes from links must have a node description (match the number of nodes from links with number of annotated nodes)

			Every check but the names in the message runs as SQL against the two tables and gives a number
			and a few ids: the nodes and the links were read into python containers before, and every node
			under root was written into the text of one query, which took 221 seconds on a database of 2.7
			million nodes and never finished the purge of a real one (issue #107).

		------
		Returns
			True: if all edges has one and only one parent
		'''
		logger.info("Count nodes, links and the tree under root")
		under_root = self._root_walk_table()
		num_nodes = self._table_count("SELECT COUNT(*) FROM nodes;")
		num_links = self._table_count("SELECT COUNT(*) FROM tree;")
		num_under = self._table_count("SELECT COUNT(*) FROM {u};".format(u=under_root))
		num_under_links = self._table_count(
			"SELECT COUNT(*) FROM tree WHERE EXISTS (SELECT 1 FROM {u} WHERE {u}.id = tree.parent);".format(u=under_root))
		num_endpoints = self._table_count(
			"SELECT COUNT(*) FROM (SELECT child AS id FROM tree UNION SELECT parent AS id FROM tree);")

		logger.info("Validate parents")
		res = self.check_parent()
		extra = ""
		if len(res) > 0:
			failed_nodes = self._node_names([x[0] for x in res])
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
					""".format(nodes = num_nodes,
								links=num_links,
								tnodes=num_under,
								tlinks = num_under_links,
								link_nodes = num_endpoints
		)
		logger.info(stats)

		try:
			## nodes that no link under root reaches: they never reach a dump or a classifier database
			self._validate_ids("The number of nodes and the number nodes under root does not match!",
				"SELECT nodes.id FROM nodes WHERE NOT EXISTS"
				" (SELECT 1 FROM {u} WHERE {u}.id = nodes.id)".format(u=under_root))
			## links whose parent is not under root, named by the node they cut off
			self._validate_ids("The number of edges and the number edges under root does not match!",
				"SELECT tree.child AS id FROM tree WHERE NOT EXISTS"
				" (SELECT 1 FROM {u} WHERE {u}.id = tree.parent)".format(u=under_root))
			## the two ends of a link that have no node description of their own
			self._validate_ids("The number of annotated nodes does not match with the number of nodes connected to edges!",
				"SELECT ends.id FROM (SELECT child AS id FROM tree UNION SELECT parent AS id FROM tree) AS ends"
				" WHERE NOT EXISTS (SELECT 1 FROM nodes WHERE nodes.id = ends.id)")
			## a node under root that has no description: the last check of the original order
			self._validate_ids("There are nodes in the database not connected to the tree",
				"SELECT {u}.id FROM {u} WHERE NOT EXISTS"
				" (SELECT 1 FROM nodes WHERE nodes.id = {u}.id)".format(u=under_root))
		finally:
			self.drop_temp_tables([under_root])
		logger.info("Validation OK!")
		return True

	def statistics(self):
		'''Print statistics and the stored build report of the database

		------
		Returns:
			true
		'''
		logger.info("Counting database rows")
		stats = """Tree statistics
					Nodes: {nodes}
					Links: {links}
					Genomes: {ngenomes}
					""".format(nodes = self.num_rows('nodes'),
								links=self.num_rows('tree'),
								ngenomes=self.num_rows('genomes')
		)
		print(stats)
		for line in self.report_lines():
			print(line)
		return True

	def report_lines(self,limit=20):
		'''The build report stored in the database: which flextaxd versions built and changed it, and
				what the last run that used a genome path or an input file found missing

		------
		Parameters:
			int - entries shown for each kind of finding, None shows all of them
		------
		Returns:
			list - lines to print
		'''
		build = self.get_build_info()
		lines = ["Database built with flextaxd: {v}".format(v=build["built_with"] or "unknown, the database was built by an earlier flextaxd version"),
				"Database last changed with flextaxd: {v}".format(v=build["modified_with"] or "unknown")]
		if build["built"]:
			lines.append("Database created: {t}".format(t=build["built"]))
		if build["modified"]:
			lines.append("Database last changed: {t}".format(t=build["modified"]))
		report_date = self.get_meta("report_date")
		if report_date:
			lines.append("Disk report from {d}, genome path {p}".format(d=report_date,p=self.get_meta("genomes_path") or "unknown"))
		for kind,description in [("missing_on_disk","Genome annotations without sequence file"),
						("sequence_annotation","Sequence level annotations (packed file, no genome file expected)"),
						("unannotated_file","Sequence files without genome annotation"),
						("unannotated_sequence","Library sequences without taxid"),
						("skipped_import","Rows skipped by the taxonomy import")]:
			total, findings, legacy = self.report_findings(kind)
			if not total:
				continue
			lines.append("{d}: {n}".format(d=description,n=total))
			lines.extend(["\t{line}".format(line=line) for line in findings])
			if legacy:
				## databases written by an earlier flextaxd version still hold the entries themselves
				show = legacy if limit is None else legacy[:limit]
				lines.extend(["\t{e}".format(e=entry) for entry in show])
				if limit is not None and len(legacy) > limit:
					lines.append("\tand {n} more entries stored by an earlier flextaxd version".format(n=len(legacy)-limit))
		return lines

	def check_parent(self):
		'''check tree structure parent
			make sure no tree node has multiple parents
		------
		Returns
			True: if node has one and only one parent
		'''
		nrank = 0
		child_w_dp = []
		try:
			QUERY = "SELECT child FROM tree GROUP BY child HAVING count(parent) > 1"  ## Thanks to andrewjmc@github for this suggestion
			logger.debug(QUERY)
			res = self.query(QUERY)
			if isinstance(res,Exception):
				raise TreeError("The tree structure could not be checked: {err}".format(err=res))
			child_w_dp = res.fetchall()
			logger.debug(child_w_dp)
			## Bugfix  only allowed one parent to be different, now multiple is ok as long as ranks are different
			for child in map(str,[x[0] for x in child_w_dp]):
				QUERY = "SELECT parent,rank_i FROM tree WHERE child = ?"
				p_ranks = self.select(QUERY,params=(child,))
				if len(p_ranks) != len(set(p_ranks)):
					logger.error("Nodes with two parents have identical ranks, Fatal Error")
				else:
					nrank += 1
			if nrank > 0:
				logger.info("Found: {n} nodes with multiple parents, however in different lineages, OK.".format(n=nrank))
				return []
			return child_w_dp
		except AttributeError as ae:
			logger.info("AttributeError occured")
			logger.info(QUERY)
			logger.info(child_w_dp)
			logger.debug(ae)
		except TypeError as ae:
			logger.info("TypeError occured")
			logger.info(QUERY)
			logger.info(child_w_dp)
			logger.debug(ae)
		raise TreeError("Following Nodes: {node} has more than one parent!".format(node=child_w_dp))

	'''Get functions of class'''
	def get_all(self, database=False, table=False,sort=False):
		'''Get full table from table

		------
		Returns
			list - results from table
		'''
		if not table:
			raise Exception("Table was not supplied to DatabaseFuction get_all()")
		QUERY = '''SELECT * FROM {table}'''.format(table=table)
		if sort:
			QUERY += " ORDER BY {sort} ASC".format(sort=sort)
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
		for node_i,genome in self.query(QUERY).fetchall():
			genomeDict[genome] = node_i
		return genomeDict

	def get_sequence_annotations(self, database=False, reference="nt"):
		'''The genome names of the annotations that were made on a sequence instead of on a genome: the
				packed multi sequence route (ReadTaxonomyNCBI.parse_nt_file, the third column of an
				accession2taxid file through ReadTaxonomy.parse_genomeid2taxid) stores the source of such
				a row as its reference, and by the decision of 2026-11-05 those files are never split up -
				the sequence stays inside the packed file and its taxid travels in the map that goes with
				library.fna. A row of this kind therefore has no sequence file of its own in a genome
				folder and must not be counted or deleted as one that is missing (issue #131).

		------
		Parameters:
			str - the reference value that marks a sequence level annotation, "nt" as the importers write it
		------
		Returns:
			set - the genome names annotated on sequence level
		'''
		rows = self.select('''SELECT genome FROM genomes WHERE reference = ?''', params=(reference,))
		return set(row[0] for row in rows)

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

	def get_nodes(self, database=False,col=5):
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
			if col == 1 or col == 5:
				nodeDict[node[0]] = node[1]
			if col == 2 or col == 5:
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
		'''Add rank level to the rank table, or return the level that already carries that name

		Returns
		------
		int - id of the rank level
		'''
		if rank == None:
			raise AttributeError("Rank cannot be None")
		'''A rank is named, not numbered: every caller asks for a rank by name and works with the id it gets
		   back, so one row per name is what the tree needs. Inserting unconditionally handed every import
		   against the same database its own copy of every rank it named, which left one level under several
		   ids (issue #138). There is no unique constraint on the column: the databases already on disk hold
		   duplicates, and an index would refuse them. Such a level keeps its lowest id, which is the one new
		   links hang on; the old duplicates are left where they are.'''
		levels = self.select("SELECT rank_i FROM rank WHERE rank = ? ORDER BY rank_i",params=(rank,))
		if levels:
			logger.debug("Rank level already present: {rank} [{rank_i}]".format(rank=rank,rank_i=levels[0][0]))
			return levels[0][0]
		info = { "rank": rank }
		### If ID is supplied skip autoincrement and add specific ID
		'''the column of the rank table is rank_i, there is no id: writing "id" made sqlite answer "no such
		   column: id" and query handed that exception back as the rank id of the new level (issue #139)'''
		if id:
			info["rank_i"] = id
		rank_id = self.insert(info, table="rank")
		#logger.debug("rank added: {info}".format(info=info))
		return rank_id

	def add_link(self, child=None, parent=None, rank=None, table="tree"):
		'''Add relationship in tree

		The rank is the id of a rank level, or None when nothing names this level of the tree: there is no
		default rank id. rank_i 1 was, and it was the "no rank" level of every database derived from an NCBI
		taxdump (the first rank row of NCBI_base.db is rank_i 1, "no rank"), so a level that had no rank
		silently ended up named with that (issue #142).

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

	def add_genome(self, genome, _id=False,reference=False):
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
		if reference:
			info["reference"] = reference
		logger.debug(info)
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
		logger.info("Slow clean, deleting {links} links from {table}".format(links=len(links),table=table))
		QUERY = "DELETE FROM {table} WHERE parent = ? AND child = ? AND ifnull(rank_i,-1) = ?"
		num_deleted = 0
		for parent,child,rank_i in links:
			logger.debug("{}-{} rank: {} deleted!".format(parent,child,rank_i))
			res = self.query(QUERY.format(table=table), insert_val=(parent,child,-1 if rank_i is None else rank_i),error=True)
			if isinstance(res,Exception):
				logger.error("Error deleting link {}-{} from {table}: {err}".format(parent,child,table=table,err=str(res)))
				return False
			num_deleted += res.rowcount
		logger.info("{num} links deleted!".format(num=num_deleted))
		## Commit changes
		if not hold:
			logger.debug("Commit changes!")
			self.commit()
		return True

	def fast_delete_links(self,links,table="tree",hold=False):
		'''This function is used when general clean function is executed.
		    Links are removed as exact (parent,child,rank_i) combinations through a temporary table,
		    the query stays the same size no matter how many links are removed

		Returns
		------
			boolean
		'''
		links = [tuple(link)[:3] for link in links]
		logger.info("Fast clean, deleting {nlinks} links from {table}".format(nlinks=len(links),table=table))
		if not links:
			return True
		temp_table = self._temp_link_table("ftx_temp_{table}_links".format(table=table),links)
		if not temp_table:
			return False
		QUERY = ("DELETE FROM {table} WHERE exists (select 1 from {temp}"
				" where {temp}.parent = {table}.parent"
				" and {temp}.child = {table}.child"
				" and ifnull({temp}.rank_i,-1) = ifnull({table}.rank_i,-1))").format(table=table,temp=temp_table)
		logger.debug(QUERY)
		try:
			res = self.cursor.execute(QUERY)
			num_deleted = res.rowcount
		except Exception as e:
			logger.error("Error deleting links from {table}: {err}".format(table=table,err=str(e)))
			self.drop_temp_tables([temp_table])
			return False
		self.drop_temp_tables([temp_table])
		logger.info("{num} links deleted!".format(num=num_deleted))
		## Commit changes
		if not hold:
			logger.debug("Commit changes!")
			self.commit()
		return True

	def ambigious_delete_links(self,nodes,hold=False):
		'''This funtion takes as input a list of nodes. It will search the "tree" table for rows that return on "parent in list_of_nodes" and "child in list_of_nodes"
		Returns
		------
			boolean
		'''
		logger.info("Entered ambigious_delete_links")
		QUERY = "DELETE FROM {table} WHERE parent in ({parents}) AND child in ({children})"
		logger.info("Attempting delete where parent&child in: ")
		logger.debug( str(nodes))
		res = self.query(QUERY.format(table='tree', parents=",".join(list(map(str,nodes))), children=",".join(list(map(str,nodes)))))
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
		nodes = list(nodes)
		logger.info("Deleting {nnodes} nodes from {table}!".format(nnodes=len(nodes),table=table))
		if not nodes:
			return True
		if len(nodes) > MAX_BIND_PARAMS:
			## Too many id values to pass as bind variables in one query, use a temporary table
			temp_table = self._temp_value_table("ftx_temp_{table}_ids".format(table=table),nodes)
			if not temp_table:
				return False
			QUERY = "DELETE FROM {table} WHERE id in (select val from {temp})".format(table=table,temp=temp_table)
			logger.debug(QUERY)
			try:
				res = self.cursor.execute(QUERY)
				num_deleted = res.rowcount
			except Exception as e:
				logger.error("Error deleting nodes from {table}: {err}".format(table=table,err=str(e)))
				self.drop_temp_tables([temp_table])
				return False
			self.drop_temp_tables([temp_table])
		else:
			QUERY = "DELETE FROM {table} WHERE id in ({nodes})".format(table=table, nodes=",".join(["?"]*len(nodes)))
			logger.debug(QUERY)
			res = self.query(QUERY, insert_val=tuple(nodes),error=True)
			if isinstance(res,Exception):
				logger.error("Error deleting nodes from {table}: {err}".format(table=table,err=str(res)))
				return False
			num_deleted = res.rowcount
		logger.info("{num} nodes deleted!".format(num=num_deleted))
		## Commit changes
		if not hold:
			logger.debug("Commit changes!")
			self.commit()
		return True

	def _temp_value_table(self,name,values):
		'''Fill a temporary table with a large list of values, used to sidestep the
		    sqlite limit on the number of bind variables in a single query

		Returns
		------
			str - 	name of the temporary table
			boolean - False if the table could not be created or filled
		'''
		self.query("DROP TABLE IF EXISTS {name}".format(name=name))
		try:
			self.cursor.execute("CREATE TEMP TABLE {name} (val primary key)".format(name=name))
			self.cursor.executemany("INSERT OR IGNORE INTO {name} (val) VALUES (?)".format(name=name),[(value,) for value in values])
			return name
		except Exception as e:
			logger.warning("Could not fill temporary table {name}: {err}".format(name=name,err=str(e)))
			return False

	def _temp_link_table(self,name,links):
		'''Fill a temporary table with links (parent,child,rank_i), used to remove many links
		    with a single query, without building a list of values in the query text

		Returns
		------
			str - 	name of the temporary table
			boolean - False if the table could not be created or filled
		'''
		self.query("DROP TABLE IF EXISTS {name}".format(name=name))
		try:
			self.cursor.execute("CREATE TEMP TABLE {name} (parent integer, child integer, rank_i integer, primary key (parent,child,rank_i))".format(name=name))
			self.cursor.executemany("INSERT OR IGNORE INTO {name} (parent,child,rank_i) VALUES (?,?,?)".format(name=name),[tuple(link)[:3] for link in links])
			return name
		except Exception as e:
			logger.warning("Could not fill temporary table {name}: {err}".format(name=name,err=str(e)))
			return False

	def drop_temp_tables(self,names):
		'''Remove temporary value tables created by _temp_value_table

		Returns
		------
			boolean
		'''
		for name in names:
			self.query("DROP TABLE IF EXISTS {name}".format(name=name))
		return True

	def delete_rows_not_in(self, table, column, temp_table, hold=False):
		'''Delete all rows of a table where the column value is not present in a temporary
		    table with a single "val" column. The query size is independent of the number of
		    values that are kept, so it works for tables of any size

		Returns
		------
			int - 	number of rows deleted
			boolean - False if the query failed
		'''
		QUERY = "DELETE FROM {table} WHERE {column} not in (select val from {temp})".format(table=table,column=column,temp=temp_table)
		logger.debug(QUERY)
		try:
			res = self.cursor.execute(QUERY)
			num_deleted = res.rowcount
		except Exception as e:
			logger.error("Error removing rows from {table}: {err}".format(table=table,err=str(e)))
			return False
		logger.info("{num} rows removed from {table}".format(num=num_deleted,table=table))
		## Commit changes
		if not hold:
			logger.debug("Commit changes!")
			self.commit()
		return num_deleted

	def delete_genomes(self,nodes,genomes='',match_genome_only=False,table="genomes",hold=False):
		'''This function clean up genomes from given nodes

		Returns
		------
			boolean
		'''
		nodes = list(nodes)
		genomes = list(genomes) if genomes else []
		if not genomes and not nodes:
			return True
		conditions = []
		if genomes: # check if genomes were specified for targeted deletion (sometimes a node can hold multiple genomes)
			conditions.append(("genome",genomes))
			if not match_genome_only and nodes:
				conditions.append(("id",nodes))
		else:
			conditions.append(("id",nodes))

		if match_genome_only and genomes:
			logger.info("Deleting annotations of {ngenomes} genomes!".format(ngenomes=len(genomes)))
		else:
			logger.info("Deleting {nnodes} annotations!".format(nnodes=len(nodes)))
		if sum(len(values) for column,values in conditions) > MAX_BIND_PARAMS:
			## Too many values to pass as bind variables in one query, use temporary tables instead
			temp_tables = []
			where = []
			for column,values in conditions:
				temp_table = self._temp_value_table("ftx_temp_{table}_{column}".format(table=table,column=column),values)
				if not temp_table:
					self.drop_temp_tables(temp_tables)
					return False
				temp_tables.append(temp_table)
				where.append("{column} in (select val from {temp})".format(column=column,temp=temp_table))
			QUERY = "DELETE FROM {table} WHERE {where}".format(table=table,where=" AND ".join(where))
			logger.debug(QUERY)
			try:
				res = self.cursor.execute(QUERY)
				num_deleted = res.rowcount
			except Exception as e:
				logger.error("Error deleting genomes from {table}: {err}".format(table=table,err=str(e)))
				self.drop_temp_tables(temp_tables)
				return False
			self.drop_temp_tables(temp_tables)
		else:
			where = " AND ".join(["{column} in({ph})".format(column=column,ph=",".join(["?"]*len(values))) for column,values in conditions])
			QUERY = "DELETE FROM {table} WHERE {where}".format(table=table,where=where)
			values = tuple([value for column,vals in conditions for value in vals])
			logger.debug(QUERY)
			res = self.query(QUERY, insert_val=values,error=True)
			if isinstance(res,Exception):
				logger.error("Error deleting genomes from {table}: {err}".format(table=table,err=str(res)))
				return False
			num_deleted = res.rowcount
		logger.info("{num} genome annotations deleted!".format(num=num_deleted))

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

		col=1 gives rank_i -> rank name; col=2 gives both directions in one dictionary (rank_i -> name and
		name -> rank_i), which is what the tree modification code works with. Ordered by rank_i, so a
		database that still holds duplicate rank rows from before issue #138 answers the same way every run.

		Returns
		------
			dict - rank to rank_index and vice versa
		'''
		QUERY = "SELECT rank_i,rank FROM rank ORDER BY rank_i"
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
		if not parents:
			return children
		QUERY = '''SELECT child FROM tree WHERE parent in({nodes})'''.format(nodes=",".join(map(str,list(parents))))

		if selected:
			QUERY = '''SELECT child FROM tree WHERE parent in({nodes}) and rank_i = {rank}'''.format(nodes=",".join(map(str,list(parents))), rank=selected)
		logger.debug(QUERY)
		res = self.select(QUERY)
		if (len(res) + len(children)) != 0:
			children = set([child_i[0] for child_i in res])
			if level < maxdepth:
				try:
					children |= self.get_children(parents=children,level=level+1,maxdepth=maxdepth)
				except RecursionError:
					logger.debug('Maximum depth reached, halting tree traversal')
		return children

	def get_genome_nodes(self, table="genomes"):
		'''Get all node ids that carry at least one genome annotation

		Returns
		------
			list - node ids
		'''
		QUERY = "SELECT DISTINCT id FROM {table}".format(table=table)
		logger.debug(QUERY)
		return [row[0] for row in self.select(QUERY)]

	def get_ancestor_nodes(self, nodes, table="tree", max_depth=50):
		'''Collect all parents above the given nodes, walking the tree in the database
		    instead of one query per level

		Returns
		------
			list - the given nodes and all their parent node ids
			boolean - False if the query failed
		'''
		seed = self._temp_value_table("ftx_temp_seed_nodes", nodes)
		if not seed:
			return False
		QUERY = ("WITH RECURSIVE path(id,depth) AS ("
				" SELECT val,0 FROM {seed}"
				" UNION"
				" SELECT {table}.parent,path.depth+1 FROM {table} JOIN path ON {table}.child = path.id"
				" WHERE path.depth < {maxdepth})"
				" SELECT DISTINCT id FROM path").format(seed=seed, table=table, maxdepth=max_depth)
		logger.debug(QUERY)
		try:
			res = [row[0] for row in self.cursor.execute(QUERY).fetchall()]
		except Exception as e:
			logger.error("Could not collect parent nodes: {err}".format(err=str(e)))
			res = False
		self.drop_temp_tables([seed])
		return res

	def get_descendant_nodes(self, nodes, table="tree", max_depth=50):
		'''Collect all nodes below the given nodes, walking the tree in the database

		Returns
		------
			list - all node ids below the given nodes (the given nodes themselves are not included)
			boolean - False if the query failed
		'''
		seed = self._temp_value_table("ftx_temp_seed_nodes", nodes)
		if not seed:
			return False
		QUERY = ("WITH RECURSIVE path(id,depth) AS ("
				" SELECT val,0 FROM {seed}"
				" UNION"
				" SELECT {table}.child,path.depth+1 FROM {table} JOIN path ON {table}.parent = path.id"
				" WHERE path.depth < {maxdepth})"
				" SELECT DISTINCT id FROM path WHERE id not in (select val from {seed})").format(seed=seed, table=table, maxdepth=max_depth)
		logger.debug(QUERY)
		try:
			res = [row[0] for row in self.cursor.execute(QUERY).fetchall()]
		except Exception as e:
			logger.error("Could not collect child nodes: {err}".format(err=str(e)))
			res = False
		self.drop_temp_tables([seed])
		return res

		return children

	def get_parent(self,name,all=False):
		'''Get parent from node id parent

		Returns
		------
			list - parent link and rank
		'''
		#QUERY = '''SELECT parent,child,rank FROM tree LEFT JOIN rank on (tree.rank_i = rank.rank_i) WHERE child = "{node}"'''.format(node=name)
		QUERY = '''SELECT parent,child,rank_i FROM tree WHERE child = ?'''
		logger.debug(QUERY)
		if all:
			res = self.select(QUERY,params=(name,))
		else:
			rows = self.select(QUERY,params=(name,))
			res = rows[0] if rows else None
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
		if not name:
			return parents
		if len(name & set([1])) > 0 and not ncbi:
			#print(name)
			return parents
		QUERY = '''SELECT parent,child{simple} FROM tree WHERE child in ({node})'''.format(simple=simple, node=",".join(map(str,name)))
		#if find_all: logger.info(QUERY)
		res = self.select(QUERY)
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
		QUERY = '''SELECT id FROM nodes WHERE name = ? COLLATE NOCASE'''
		try:
			res = self.query(QUERY,insert_val=(name,),error=True).fetchone()[0]
		except TypeError:
			logger.debug(QUERY)
			raise NameError("Name not found in the database! {name}".format(name=name))
		return res

	def get_name(self,id):
		'''get node name from id

		Returns
		------
			char - node id from node name
		'''
		QUERY = '''SELECT name FROM nodes WHERE id = ?'''
		try:
			res = self.query(QUERY,insert_val=(id,),error=True).fetchone()[0]
		except TypeError:
			logger.debug(QUERY)
			raise NameError("Name not found in the database! {id}".format(id=id))
		return res

	def update_table(self,data,table):
		'''Add annotation to table

		Returns
		------
			see update responses
		'''
		return self.update(data, table=table)

	def update_genome(self,data):
		'''Add genome annotation to nodes

		Returns
		------
			see update responses
		'''
		return self.update_table(data, table="genomes")
