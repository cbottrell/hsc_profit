import pandas as pd
import numpy as np
import os

seed = 12345
np.random.seed(seed)

pdr_usr = os.getenv('SSP_PDR_USR')
pdr_pwd = os.getenv('SSP_PDR_PWD')

bands = ['HSC-G','HSC-R','HSC-I','HSC-Z','HSC-Y']

df = pd.read_csv('catalogues/Final_Tracts-Patches_pdr3_wide.csv')

nsamples = 1000
index = np.random.choice(np.arange(len(df)),nsamples,replace=False)

df = df.loc[index]
tracts = df['Tracts']
patches = df['Patches']

for tract,patch in zip(tracts,patches):

    patch = list(f'{patch:03}')
    patch[1]=','
    patch = ''.join(patch)
    
    for band in bands:

        filename = f'images/{band}-{tract}-{patch}.fits'
        print(filename)

        if not os.access(filename,0):
            url = f'https://hsc-release.mtk.nao.ac.jp/archive/filetree/pdr3_wide/deepCoadd-results/{band}/{tract}/{patch}/calexp-{band}-{tract}-{patch}.fits'
            print(url)
            os.system(f'wget -nv --user={pdr_usr} --password {pdr_pwd} {url} -O {filename}')