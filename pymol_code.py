from pymol import cmd
import matplotlib.pyplot as plt
import os
import pandas as pd

#Global variables
data_path = 'data/TTR'

def visualize_electrostatics(image_name):
    """
    Visualize Poisson-Boltzmann electrostatics on the surface of a protein.
    """

    #cmd.show('surface', 'protein')
    # Generate vacuum electrostatics
    cmd.vacuum('protein')  # mode=2 for protein contact potential
    # For example, cmd.color('red', 'protein and elem O') to color oxygen atoms red
    cmd.png(image_name)
    cmd.reinitialize()


def calculate_sasa():
    """
    Calculate the solvent-accessible surface area (SASA) of a PDB structure.
    """
    # Set parameters for SASA calculation
    cmd.set('dot_solvent', 1)  # Enable SASA calculation
    # Calculate SASA
    sasa = cmd.get_area('protein')
    print(f"SASA: {sasa} Å²")
    
    return sasa

def individual_analysis(file):
    image_name = f'{data_path}/images/{file[:-4]}'
    cmd.load(f'{data_path}/{file}', 'protein')
    cmd.png(image_name) 

    sasa = calculate_sasa()

    # figure out
    # visualize_electrostatics(f'electrostatics_{image_name}.png')
    cmd.reinitialize()

    return sasa

def calculate_rmsd(protein1, protein2):
    cmd.align(protein1, protein2)
    rmsd = cmd.rms_cur(protein1, protein2)
    return rmsd

def main():
    global data_path
    pdb_files = [f for f in os.listdir(data_path) if f.endswith('.pdb')]

    # Calculate SASA for each file
    sasa_data = {}
    for file in pdb_files:
        print(file)
        sasa = individual_analysis(file)
        sasa_data[file] = sasa

    # Save SASA data to a CSV file
    sasa_df = pd.DataFrame.from_dict(sasa_data, orient='index', columns=['SASA'])
    sasa_df.to_csv(f'{data_path}/dfs/sasa_data.csv')

    # Calculate RMSD for each pair of files
    rmsd_scores = {}
    n_files = len(pdb_files)
    for i in range(n_files):
        for j in range(i + 1, n_files):
            file1 = pdb_files[i]
            file2 = pdb_files[j]
            cmd.load(os.path.join(data_path, file1), 'protein1')
            cmd.load(os.path.join(data_path, file2), 'protein2')
            rmsd_score = calculate_rmsd('protein1', 'protein2')
            rmsd_scores[(file1, file2)] = rmsd_score
            cmd.delete('protein1')
            cmd.delete('protein2')

    # Save RMSD data to a CSV file
    rmsd_df = pd.DataFrame.from_dict(rmsd_scores, orient='index', columns=['RMSD'])
    rmsd_df.to_csv(f'{data_path}/dfs/rmsd_data.csv')
    
    return sasa_df, rmsd_df


cmd.extend('main', main)