#!/usr/bin/env python3 -c
import sys
import os
import sqlite3

'''Create SQL database'''

sql_create_nodes_table = """ CREATE TABLE IF NOT EXISTS nodes (
                                    id integer PRIMARY KEY,
                                    name text NOT NULL
                                ); """

sql_create_tree_table = """CREATE TABLE IF NOT EXISTS tree (
                                parent integer NOT NULL,
                                child integer NOT NULL,
                                rank_i integer,
                                FOREIGN KEY (parent) REFERENCES nodes (id),
                                FOREIGN KEY (child) REFERENCES nodes (id),
                                FOREIGN KEY (rank_i) REFERENCES rank (rank_i),
                                unique (parent, child)
                            );"""

sql_create_rank_table = """CREATE TABLE IF NOT EXISTS rank (
                                rank_i integer PRIMARY KEY,
                                rank VARCHAR(15)
                            );"""

sql_create_genomes_table = """ CREATE TABLE IF NOT EXISTS genomes (
                                    id integer NOT NULL,
                                    genome text NOT NULL,
                                    FOREIGN KEY (id) REFERENCES nodes (id)
                                ); """


def create_connection(db_file):
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

def create_table(conn, create_table_sql):
    """ create a table from the create_table_sql statement
    :param conn: Connection object
    :param create_table_sql: a CREATE TABLE statement
    :return:
    """
    try:
        c = conn.cursor()
        c.execute(create_table_sql)
    except Exception as e:
        print(e)

def add_table(database, table):
    conn = create_connection(database)
    create_table(sql_create_genomes_table)
    return

def main(database="database.db"):
    # create a database connection
    conn = create_connection(database)
    if conn is not None:
        # create nodes table
        create_table(conn, sql_create_nodes_table)
        # create tree table
        create_table(conn, sql_create_tree_table)
        # create genomes table
        create_table(conn, sql_create_genomes_table)
        # create rank tables
        create_table(conn, sql_create_rank_table)
    else:
        print("Error! cannot create the database connection.")

if __name__ == '__main__':
    main(sys.argv[1])
