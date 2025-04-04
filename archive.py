import os,sys

universe = 'Simba'
simulation = 'Simba100-1'

morph ='sersic'

file_path = f'/scratch/pawsey0119/bottrell/Simulations/virgotng/data/{universe}/{simulation}/postprocessing/skirt_images_hsc/realistic'
out_path = f'{file_path}/fits/{morph}'
if not os.access(out_path,0): 
    os.system(f'mkdir -p {out_path}')
os.chdir(file_path)

snapnum = int(os.environ['SLURM_ARRAY_TASK_ID'])

# os.system(f'tar -xf {universe}_{simulation}_{snapnum:03}_hsc.tar -C .')
# keep only v0 camera
os.system(f'tar -cf {out_path}/profuse_{snapnum:03}.tar {snapnum:03}/profuse*v0*.fits')
# os.system(f'rm {snapnum:03}/profuse*.fits')
