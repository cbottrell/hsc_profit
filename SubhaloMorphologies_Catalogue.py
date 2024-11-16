import pymysql
import pandas as pd

def query_df(command,database,cnf_path):
    '''Query database and return results as dataframe.'''
    connection = pymysql.connect(
        database=database,
        read_default_file=cnf_path,
    )
    df = pd.read_sql(command, con=connection)
    connection.close()
    #df = _df_nullint(df)
    return df

dbcmd = [
    'SELECT s.*, m.*',
    'FROM Subhalos as s JOIN Morphologies_ProFit_HSC_Sersic as m',
    'ON s.snapnum=m.snapnum AND s.subfindid=m.subfindid',
    'WHERE s.SubhaloMassType_stars>=9 and s.snapnum=91',
    'ORDER BY s.snapnum,s.subfindid,m.camera,m.band',
]

df = query_df(' '.join(dbcmd), database='IllustrisTNG50_1',cnf_path='/home/bottrell/.mysql/ningaloo.cnf')