from seqlib import read_fastq
from finder import pull_AD
import sqlite3
from sqlite3 import Error, connect
import pandas as pd

def initialize_record(fastq, db_file):
    conn = sqlite3.connect(db_file)
    for line in read_fastq(fastq):
        read = line[1]
        tile, barcode = pull_AD(read)


def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by db_file
    :param db_file: database file
    :return: Connection object or None
    """

    conn = None
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)

    return conn

def create_table(conn, create_table_sql):
    """ create a table from the create_table_sql statement
    :param conn: Connection object
    :param create_table_sql: a CREATE TABLE statement
    :return:
    """

    try:
        c = conn.cursor()
        c.execute(create_table_sql)
    except Error as e:
        print(e)

def create_link(conn, link):
    """
    Create a new project into the projects table
    :param conn:
    :param link:
    :return: link id
    """
    sql = ''' INSERT INTO ADbarcodes(AD,barcode)
              VALUES(?,?) '''
    cur = conn.cursor()
    cur.execute(sql, link)
    conn.commit()
    return cur.lastrowid

def insert_from_file(conn, file):
    """Add all tile, barcode pairs from fastq file into db."""

    links = []
    for line in read_fastq(file):
        AD, barcode = pull_AD(line[1])
        if AD != None:
            links.append((AD, barcode))
    
    cur = conn.cursor()
    sql = """INSERT INTO ADbarcodes
                        (AD, barcode) 
                        VALUES (?, ?);"""

    cur.executemany(sql, links)
    conn.commit()
    print("Total", cur.rowcount, "Records inserted successfully into ADbarcodes table")
    conn.commit()

def count_unique(conn):
    """SELECT *, COUNT(barcode) FROM ADbarcodes GROUP BY AD, barcode;"""
    pass

def export_to_pandas_df(conn, sql = "SELECT *, COUNT(barcode) FROM ADbarcodes GROUP BY AD, barcode;"):
    """Make the table into a pandas dataframe."""

    df = pd.read_sql_query(sql, conn)

    return df

def export_to_csv(conn, name):
    """Make the table into csv file."""
    df = export_to_pandas_df(conn)

    df.to_csv(f"{name}.csv", index = False)

# SAVED QUERIES
# to get the counts of each AD, barcode pair
# SELECT *, COUNT(barcode) FROM ADbarcodes GROUP BY AD, barcode;
# to get the number of total reads parsed
# SELECT SUM(count) FROM (SELECT *, COUNT(barcode) as count FROM ADbarcodes GROUP BY AD, barcode);


# variables
def initialize_db(name):
    database = fr"C:\sqlite\db\{name}.db"

    sql_create_ADbarcodes_table = """ CREATE TABLE IF NOT EXISTS ADbarcodes (
                                        AD text NOT NULL,
                                        barcode text NOT NULL
                                    ); """

    # sql_create_AD_table = """CREATE TABLE IF NOT EXISTS ADs (
    #                                 AD text PRIMARY KEY,
    #                                 count integer,
    #                                 FOREIGN KEY (AD) REFERENCES ADbarcodes (AD)
    #                             );"""

    # sql_create_barcodes_table = """CREATE TABLE IF NOT EXISTS barcodes (
    #                                 barcode text PRIMARY KEY,
    #                                 count integer,
    #                                 FOREIGN KEY (barcode) REFERENCES ADbarcodes (barcode)
    #                             );""" 

    conn = create_connection(database)

    # create tables
    if conn is not None:
        # create AD_barcodes table
        create_table(conn, sql_create_ADbarcodes_table)

        # # create ADs table
        # create_table(conn, sql_create_AD_table)

        # # create barcodes table
        # create_table(conn, sql_create_barcodes_table)
    else:
        print("Error! cannot create the database connection.")

    return conn

def main():

    conn = initialize_db("rep1")

    insert_from_file(conn, "LC_M1_checkstep1_rep1.fastq")
    
    export_to_csv(conn, "rep1")


if __name__ == '__main__':
    main()