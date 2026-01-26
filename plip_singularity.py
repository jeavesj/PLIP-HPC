import argparse
import subprocess
import numpy as np
import pandas as pd
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True, help='Path to input pdb file (if receptor only, also specify ligand file path with -l)')
    parser.add_argument('-l', type=str, required=False, help='Optional path to ligand file (SDF, mol2, or PDB) if input PDB file is receptor only')
    parser.add_argument('-o', type=str, required=False, help='Path to output xml file')
    parser.add_argument('--plip_simg_path', type=str, default='/mnt/research/woldring_lab/Software/plip_3.0.0.simg', help='path to plip singularity image')

    args = parser.parse_args()
    if not args.o:
        out_path = os.getcwd()
    else:
        out_path = args.o
    
    if args.l:
        protein_path = os.path.abspath(args.i)
        ligand_path = os.path.abspath(args.l)

        if not os.path.isfile(protein_path):
            raise FileNotFoundError(f"Protein file not found: {protein_path}")
        if not os.path.isfile(ligand_path):
            raise FileNotFoundError(f"Ligand file not found: {ligand_path}")

        os.makedirs(out_path, exist_ok=True)

        def _read_text(path: str) -> str:
            with open(path, "r") as f:
                return f.read()

        def _strip_end_records(pdb_text: str) -> list[str]:
            lines = []
            for ln in pdb_text.splitlines():
                if ln.startswith(("END", "ENDMDL")):
                    continue
                lines.append(ln)
            return lines

        def _force_chain_id(pdb_lines: list[str], chain_id: str = "Z") -> list[str]:
            # PDB chain ID is column 22 (0-based index 21)
            out = []
            for ln in pdb_lines:
                if ln.startswith(("ATOM", "HETATM")) and len(ln) >= 22:
                    ln = ln[:21] + chain_id + ln[22:]
                out.append(ln)
            return out

        def _ligand_to_pdb_lines(lig_path: str) -> list[str]:
            ext = os.path.splitext(lig_path)[1].lower()

            # If ligand already PDB, just read it
            if ext in {".pdb", ".ent"}:
                return _strip_end_records(_read_text(lig_path))

            # Try RDKit for SDF/MOL2
            try:
                from rdkit import Chem
                from rdkit.Chem import AllChem

                mol = None
                if ext in {".sdf", ".sd", ".mol"}:
                    suppl = Chem.SDMolSupplier(lig_path, removeHs=False)
                    for m in suppl:
                        if m is not None:
                            mol = m
                            break
                elif ext == ".mol2":
                    mol = Chem.MolFromMol2File(lig_path, removeHs=False)

                if mol is None:
                    raise ValueError(f"RDKit could not parse ligand: {lig_path}")

                # Ensure we have 3D coords (PLIP needs geometry)
                if mol.GetNumConformers() == 0:
                    mol = Chem.AddHs(mol)
                    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                    AllChem.UFFOptimizeMolecule(mol)

                pdb_block = Chem.MolToPDBBlock(mol)
                return _strip_end_records(pdb_block)

            except Exception:
                # Fall back to OpenBabel if available: obabel ligand.xxx -O ligand.pdb
                import tempfile

                with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, dir=out_path) as tmp:
                    tmp_pdb = tmp.name

                try:
                    subprocess.run(["obabel", lig_path, "-O", tmp_pdb], check=True)
                    lines = _strip_end_records(_read_text(tmp_pdb))
                    return lines
                finally:
                    try:
                        os.remove(tmp_pdb)
                    except OSError:
                        pass
        # Load protein and ligand PDB lines
        prot_lines = _strip_end_records(_read_text(protein_path))
        lig_lines = _ligand_to_pdb_lines(ligand_path)

        # Force ligand onto a distinct chain ID to avoid residue/chain collisions
        lig_lines = _force_chain_id(lig_lines, chain_id="Z")

        # Only keep sensible record types from ligand PDB
        keep_prefix = ("ATOM", "HETATM", "CONECT", "TER", "ANISOU", "REMARK")
        lig_lines = [ln for ln in lig_lines if ln.startswith(keep_prefix)]

        # Write combined PDB
        prot_stem = os.path.splitext(os.path.basename(protein_path))[0]
        lig_stem = os.path.splitext(os.path.basename(ligand_path))[0]
        input_file = os.path.join(out_path, f"{prot_stem}__{lig_stem}__plip_input.pdb")

        with open(input_file, "w") as f:
            for ln in prot_lines:
                f.write(ln + "\n")
            # Add TER separator if not already present at end of protein block
            if prot_lines and not prot_lines[-1].startswith("TER"):
                f.write("TER\n")
            for ln in lig_lines:
                f.write(ln + "\n")
            f.write("END\n")
    else:
        input_file = args.i

    subprocess.run([args.plip_simg_path, '-f', input_file, '-o', out_path, '-x'])

if __name__ == "__main__":
    main()