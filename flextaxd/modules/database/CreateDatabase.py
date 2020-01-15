#!/usr/bin/env python3 -c
import sys
import os
import sqlite3

class ConnectionError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

'''Create SQL database'''

class CreateDatabase(object):
    """docstring for CreateDatabase"""
    def __init__(self, database, verbose):
        super().__init__()
        self.verbose=verbose
        self.database_name = database

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
                                        unique (parent, child)
                                    );"""

        self.sql_create_rank_table = """CREATE TABLE IF NOT EXISTS rank (
                                        rank_i integer PRIMARY KEY,
                                        rank VARCHAR(15)
                                    );"""

        self.sql_create_genomes_table = """ CREATE TABLE IF NOT EXISTS genomes (
                                            id integer PRIMARY KEY,
                                            genome text NOT NULL,
                                            FOREIGN KEY (id) REFERENCES nodes (id)
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
            print(e)
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
            print(e)

    def add_table(self,table):
        self.create_table(table)
        self.conn.commit()
        return

    def create_database(self,database=False):
        # create a database connection
        self.conn = self.create_connection(database)
        if self.conn is not None:
            # create nodes table
            self.create_table(self.sql_create_nodes_table)
            # create tree table
            self.create_table(self.sql_create_tree_table)
            # create genomes table
            self.create_table(self.sql_create_genomes_table)
            # create rank tables
            self.create_table(self.sql_create_rank_table)
            
            self.conn.commit()
        else:
            raise ConnectionError("Error! cannot connet to the database {db}".format(db=database))
        return

if __name__ == '__main__':
    CDB = CreateDatabase()
    CDB.create_database(sys.argv[1])
