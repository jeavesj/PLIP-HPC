import argparse
import subprocess
import numpy as np
import pandas as pd
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True, help='path to input pdb file')
    parser.add_argument('-o', type=str, required=False, help='path to output xml file')
    parser.add_argument('--plip_simg_path', type=str, default='/mnt/research/woldring_lab/Software/plip_3.0.0.simg', help='path to input pdb file')

    args = parser.parse_args()
    if not args.o:
        out_path = os.getcwd()
    else:
        out_path = args.o
    subprocess.run([args.plip_simg_path, '-f', args.i, '-o', out_path, '-x'])

if __name__ == "__main__":
    main()