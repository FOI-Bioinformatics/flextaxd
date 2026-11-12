#!/usr/bin/env python3 -c
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

'''Create SQL database'''

class CreateDatabase(object):
    """docstring for CreateDatabase"""
    def __init__(self, verbose=False):
        super().__init__()
        self.verbose=verbose
        self.sql_create_nodes_table = """ CREATE TABLE IF NOT EXISTS nodes (
                                            id integer PRIMARY KEY,
                                            name text NOT NULL
                                        ); """

        self.sql_create_tree_table = """CREATE TABLE IF NOT EXISTS tree (
                                        parent integer NOT NULL,
                                        child integer NOT NULL,
                                        rank_i integer,
                                        FOREIGN KEY (parent) REFERENCES nodes (id),
                                        FOREIGN KEY (child) REFERENCES nodes (id),
                                        FOREIGN KEY (rank_i) REFERENCES rank (rank_i),
                                        unique (parent, child,rank_i)
                                    );"""

        ## The constraint above cannot cover a link that carries no rank: NULL never equals NULL in a
        ## sqlite unique constraint, so importing a taxonomy into a database that already held those links
        ## wrote the identical rank-less link again on every run (issue #147). An index over
        ## ifnull(rank_i,-1) covers those links the way the constraint covers the ranked ones, -1 being the
        ## same stand-in for "no rank" that deleting links already uses (DatabaseConnection.fast_delete_links).
        self.sql_create_tree_link_index = """CREATE UNIQUE INDEX IF NOT EXISTS tree_link
                                            ON tree (parent, child, ifnull(rank_i,-1));"""

        self.sql_create_rank_table = """CREATE TABLE IF NOT EXISTS rank (
                                        rank_i integer PRIMARY KEY,
                                        rank VARCHAR(15)
                                    );"""

        self.sql_create_genomes_table = """ CREATE TABLE IF NOT EXISTS genomes (
                                            id integer NOT NULL,
                                            genome VARCHAR(255) NOT NULL UNIQUE,
                                            reference VARCHAR(255),
                                            FOREIGN KEY (id) REFERENCES nodes (id)
                                        ); """

        self.sql_create_meta_table = """ CREATE TABLE IF NOT EXISTS meta (
                                            key VARCHAR(255) PRIMARY KEY,
                                            value TEXT
                                        ); """

        self.sql_create_report_table = """ CREATE TABLE IF NOT EXISTS report (
                                            kind VARCHAR(31) NOT NULL,
                                            source VARCHAR(255) NOT NULL,
                                            entry TEXT NOT NULL,
                                            UNIQUE (kind,source,entry)
                                        ); """

    def create_connection(self,db_file):
        """ create a database connection to the SQLite database
            specified by db_file
        :param db_file: database file
        :return: Connection object or None
        """
        try:
            conn = sqlite3.connect(db_file)
            return conn
        except Exception as e:
            logger.debug(e)
        return None

    def create_table(self,create_table_sql):
        """ create a table from the create_table_sql statement
        :param conn: Connection object
        :param create_table_sql: a CREATE TABLE statement
        :return:
        """
        try:
            c = self.conn.cursor()
            c.execute(create_table_sql)
        except Exception as e:
            logger.debug(e)

    def add_table(self,table):
        self.create_table(table)
        self.conn.commit()
        return

    def mark_tree_link_index(self):
        """Record in the meta table that the unique tree link index is in place. An older database gets
        the index when it is opened (DatabaseConnection.ensure_tree_link_index), which costs a look at
        every link, so the mark keeps that work to the databases that still need it."""
        try:
            c = self.conn.cursor()
            c.execute("INSERT OR REPLACE INTO meta (key,value) VALUES ('tree_link_index','1')")
        except Exception as e:
            logger.debug(e)
        return

    def create_database(self,database=False):
        # create a database connection
        self.conn = self.create_connection(database)
        if self.conn is not None:
            # create nodes table
            self.create_table(self.sql_create_nodes_table)
            # create tree table
            self.create_table(self.sql_create_tree_table)
            self.create_table(self.sql_create_tree_link_index)
            # create genomes table
            self.create_table(self.sql_create_genomes_table)
            # create rank tables
            self.create_table(self.sql_create_rank_table)
            # create the meta table holding the flextaxd version and the import report
            self.create_table(self.sql_create_meta_table)
            # create the report table holding build time findings per input
            self.create_table(self.sql_create_report_table)

            ## record that this database has the link index from the start, so that opening it does not
            ## look for the duplicate rank-less links a database of an older version can hold
            self.mark_tree_link_index()

            self.conn.commit()
        else:
            raise ConnectionError("Error! cannot connet to the database {db}".format(db=database))
        return

if __name__ == '__main__':
    CDB = CreateDatabase()
    CDB.create_database(sys.argv[1])
