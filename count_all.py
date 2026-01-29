import pandas as pd
import os
import argparse
import tarfile
from count_interactions import count_interactions

sources = ['gnina', 'gnina-apo', 'gnina-af3', 'rosetta', 'crystal', 'boltz2']
pdb_csv = '/mnt/research/woldring_lab/Members/Eaves/plip-plop/AF3/casf2016_smiles_seqs_cleaned.csv'
base_dir = '/mnt/scratch/jeaves/plip'

int_types = ['hydrophobic_interaction', 'hydrogen_bond', 'water_bridge', 'salt_bridge', 'pi_stack', 'pi_cation_interaction', 'halogen_bond']

pdbids = pd.read_csv(pdb_csv)['pdb_id'].unique()

all_counts = []
for source in sources:
    tar_path = os.path.join(base_dir, f'{source}_plips.tar.gz')
    extract_path = os.path.join(base_dir, source)
    
    with tarfile.open(tar_path, 'r') as tar:
        tar.extractall(path=extract_path, filter='data')
    
    for pdbid in pdbids:
        pdbid = pdbid.lower()
        
        if source in ['crystal', 'boltz2']:
            count_dict = None
            pdb_dir = os.path.join(extract_path, pdbid)
            
            contents = os.listdir(pdb_dir)
            
            for item in contents:
                if str(item).endswith('.xml'):
                    count_dict = count_interactions(os.path.join(pdb_dir, item))
            for int_type in int_types:
                all_counts.append({'pdbid': pdbid, 'source': source, 'pose': 1, 'int_type': int_type, 'count': count_dict[int_type]})
        else:
            for n in range(1, 101):
                count_dict = None
                poseid = str(n).zfill(4)
                
                pdb_dir = os.path.join(extract_path, f'{source}_best{poseid}', pdbid)
            
                contents = os.listdir(pdb_dir)
                
                for item in contents:
                    if str(item).endswith('.xml'):
                        count_dict = count_interactions(os.path.join(pdb_dir, item))
                if count_dict is None:
                    continue
                
                for int_type in int_types:
                    all_counts.append({'pdbid': pdbid, 'source': source, 'pose': 1, 'int_type': int_type, 'count': count_dict[int_type]})
    
df = pd.DataFrame(all_counts)

df.to_csv('/mnt/research/woldring_lab/Members/Eaves/FAIR_PLBAP/results/plip_counts.csv', index=False)