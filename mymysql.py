import pymysql
import numpy as np
import pandas as pd

def _df_nullint(df):
    '''Convert int to nullible int types.'''
    dtypes = dict(df.dtypes)
    for key in dtypes.keys():
        if 'int' in str(dtypes[key]): 
            dtypes[key] = pd.Int64Dtype()
    return df.astype(dtypes) 

def query_df(command,database,cnf_path):
    '''Query database and return results as dataframe.'''
    connection = pymysql.connect(
        database=database,
        read_default_file=cnf_path,
    )
    df = pd.read_sql(command, con=connection)
    connection.close()
    df = _df_nullint(df)
    return df

def query(command,database,cnf_path):
    '''Query database and return results as numpy array.'''
    connection = pymysql.connect(
        database=database,
        read_default_file=cnf_path,
    )
    cursor = connection.cursor()
    cursor.execute(command)
    data = np.asarray(cursor.fetchall())
    cursor.close()
    connection.close()
    return data

def submit(command,database,cnf_path):
    '''Submit command to database.'''
    connection = pymysql.connect(
        database=database,
        read_default_file=cnf_path,
        autocommit=True,
    )
    cursor = connection.cursor()
    cursor.execute(command)
    cursor.close()
    connection.close()
    return 