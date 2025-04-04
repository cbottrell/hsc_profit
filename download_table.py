import mymysql
import os

database = 'IllustrisTNG50_1'
table = 'Subhalos'
cnf_path='~/.mysql/ningaloo.cnf'

dbcmd = f'select * from {table} where snapnum>=72 and snapnum<=91 and SubhaloMassType_stars>=9 order by snapnum,subfindid'
df = mymysql.query_df(dbcmd, database=database, cnf_path=cnf_path)

filename = f'catalogues/{table}.csv'
if os.access(filename,0): os.remove(filename)
df.to_csv(filename,index=False)
