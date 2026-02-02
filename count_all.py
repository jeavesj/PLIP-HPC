import pandas as pd
import os
import argparse
import tarfile
import subprocess
from count_interactions import count_interactions

sources = ['gnina', 'gnina-apo', 'gnina-af3', 'rosetta', 'crystal', 'boltz2']
int_types = ['hydrophobic_interaction', 'hydrogen_bond', 'water_bridge', 'salt_bridge', 'pi_stack', 'pi_cation_interaction', 'halogen_bond']

base_dir = '/mnt/scratch/jeaves/plip'

pdb_csv = '/mnt/research/woldring_lab/Members/Eaves/plip-plop/AF3/casf2016_smiles_seqs_cleaned.csv'
out_csv = '/mnt/research/woldring_lab/Members/Eaves/FAIR_PLBAP/results/plip_counts.csv'

pdbids = pd.read_csv(pdb_csv)['pdb_id'].unique()
done_df = pd.read_csv(out_csv)

for source in sources:
    tar_path = os.path.join(base_dir, f'{source}_plips.tar.gz')
    extract_path = os.path.join(base_dir, source)
    
    if not os.path.exists(extract_path):
        with tarfile.open(tar_path, 'r') as tar:
            tar.extractall(path=extract_path, filter='data')
        print(f'Extracted files for {source} to {extract_path})')
    
    if source in done_df['source'].unique():
        done_pdbs = done_df[done_df['source']==source]['pdbid'].unique()

    for pdbid in pdbids:
        pdbid = pdbid.lower()
        
        if pdbid in done_pdbs:
            continue
        
        all_counts = []
        if source in ['crystal', 'boltz2']:
            count_dict = None
            pdb_dir = os.path.join(extract_path, pdbid)
            if not os.path.exists(pdb_dir):
                print(f'[WARN] Skipping {pdbid} {source} pose {poseid}. Directory does not exist: {pdb_dir}')
                continue
            
            contents = os.listdir(pdb_dir)
            for item in contents:
                if str(item).endswith('.xml'):
                    count_dict = count_interactions(os.path.join(pdb_dir, item))
            for int_type in int_types:
                all_counts.append({'pdbid': pdbid, 'source': source, 'pose': int(1), 'int_type': int_type, 'count': count_dict[int_type]})
        else:
            for n in range(1, 101):
                count_dict = None
                poseid = str(n).zfill(4)
                
                pdb_dir = os.path.join(extract_path, f'{source}_best{poseid}', pdbid)
                if not os.path.exists(pdb_dir):
                    print(f'[WARN] Skipping {pdbid} {source} pose {poseid}. Directory does not exist: {pdb_dir}')
                    continue
                contents = os.listdir(pdb_dir)
                for item in contents:
                    if str(item).endswith('.xml'):
                        count_dict = count_interactions(os.path.join(pdb_dir, item))
                if count_dict is None:
                    continue
                
                for int_type in int_types:
                    all_counts.append({'pdbid': pdbid, 'source': source, 'pose': int(n), 'int_type': int_type, 'count': count_dict[int_type]})
        
        df = pd.DataFrame(all_counts)
        df.to_csv(out_csv, mode='a', index=False)
    
    subprocess.run(['rm', '-r', extract_path])
    print(f'Removed extracted files for {source} (path: {extract_path})')
