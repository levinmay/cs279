from pymol import cmd
import matplotlib.pyplot as plt
import os
import pandas as pd
from pmg_tk.startup.apbs_gui.creating import pdb2pqr_cli
from pmg_tk.startup.apbs_gui.electrostatics import map_new_apbs
import scipy
from scipy.spatial.distance import pdist, squareform


# Set the data name to the name of the folder in the data folder that has your PDB files
data_path = 'static_folder/data/Hemo'

def calculate_tmScore(file1, file2):
    """
    Calculate the TM-score between two protein structures using TMalign.

    Args:
        file1 (str): Path to the first protein structure file.
        file2 (str): Path to the second protein structure file.

    Returns:
        float: The TM-score between the two protein structures.
    """
    return tmscore(file1, file2, args='', exe='/Users/maylevin/Desktop/cs279/app/TMscore', quiet=0)


def contact_map(file):
    """
    Generate a contact map for a protein structure.

    Args:
        file (str): Path to the protein structure file.
    """
    cmd.select("protein_chain", "chain A")
    coords = cmd.get_coords("protein_chain and name CA")
    distances = squareform(pdist(coords))
    # Plot as a heatmap
    plt.figure()
    plt.imshow(distances, cmap='viridis', origin='upper', interpolation='none')
    plt.colorbar(label='Distance (A)')
    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')
    plt.title(f'Contact Map for {file[:-4]}')
    plt.savefig(f'{data_path}/contact/{file[:-4]}_contact.png')


def calculate_sasa():
    """
    Calculate the solvent-accessible surface area (SASA) of a PDB structure.

    Returns:
        float: The SASA value in square angstroms.
    """
    # Set parameters for SASA calculation
    cmd.set('dot_solvent', 1)  # Enable SASA calculation
    # Calculate SASA
    sasa = cmd.get_area('protein')
    print(f"SASA: {sasa} Å²")
    
    return sasa

def individual_analysis(file):
    """
    Perform individual analysis on a protein structure.

    Args:
        file (str): Path to the protein structure file.

    Returns:
        float: The SASA value of the protein structure.
    """
    image_name = f'{data_path}/images/{file[:-4]}'
    cmd.load(f'{data_path}/{file}', 'protein')
    cmd.png(image_name) 
    contact_map(file)
    sasa = calculate_sasa()
    apbs(file)
    cmd.reinitialize()

    return sasa

def calculate_rmsd(protein1, protein2):
    """
    Calculate the RMSD between two protein structures.

    Args:
        protein1 (str): Name of the first protein structure.
        protein2 (str): Name of the second protein structure.

    Returns:
        float: The RMSD value between the two protein structures.
    """
    cmd.align(protein1, protein2)
    rmsd = cmd.rms_cur(protein1, protein2)
    return rmsd

def apbs(file):
    """
    Generate an APBS electrostatic potential map for a protein structure.

    Args:
        file (str): Path to the protein structure file.
    """
    pdb2pqr_cli("prepared01", "protein", options=["--ff", "amber"])
    map_new_apbs("apbs_map01", "prepared01")
    cmd.ramp_new("apbs_ramp01", "apbs_map01", [-5, 0, 5])
    cmd.set("surface_ramp_above_mode", 1, "prepared01")
    cmd.set("surface_color", "apbs_ramp01", "prepared01")
    cmd.show("surface", "prepared01")
    cmd.png(f'{data_path}/apbs/{file[:-4]}_apbs.png')
    return

def main():
    """
    Main function to perform analysis on a set of protein structures.

    Returns:
        tuple: A tuple containing the SASA DataFrame, RMSD DataFrame, and TM scores dictionary.
    """
    global data_path
    pdb_files = [f for f in os.listdir(data_path) if f.endswith('.pdb')]

    # Calculate SASA for each file
    sasa_data = {}
    for file in pdb_files:
        print(file)
        sasa = individual_analysis(file)
        sasa_data[file] = sasa

    sasa_df = pd.DataFrame.from_dict(sasa_data, orient='index', columns=['SASA'])
    sasa_df.to_csv(f'{data_path}/dfs/sasa_data.csv')

    # Calculate RMSD/TM scores for each pair of files
    rmsd_scores = {}
    TM_scores = {}

    n_files = len(pdb_files)
    for i in range(n_files):
        for j in range(i + 1, n_files):
            file1 = pdb_files[i]
            file2 = pdb_files[j]
            cmd.load(os.path.join(data_path, file1), 'protein1')
            cmd.load(os.path.join(data_path, file2), 'protein2')
            tm_score = calculate_tmScore('protein1', 'protein2')
            TM_scores[(file1, file2)] = tm_score
            rmsd_score = calculate_rmsd('protein1', 'protein2')
            rmsd_scores[(file1, file2)] = rmsd_score
            cmd.delete('protein1')
            cmd.delete('protein2')


    #Save TM scores to a CSV file
    TM_df = pd.DataFrame.from_dict(TM_scores, orient='index', columns=['TM Score'])
    TM_df.to_csv(f'{data_path}/dfs/TM_data.csv')
    
    # Save RMSD data to a CSV file
    rmsd_df = pd.DataFrame.from_dict(rmsd_scores, orient='index', columns=['RMSD'])
    rmsd_df.to_csv(f'{data_path}/dfs/rmsd_data.csv')
    
    return sasa_df, rmsd_df, TM_scores


cmd.extend('main', main)




def tmalign(mobile, target, args='', exe='TMalign', ter=0, transform=1, object=None, quiet=0):
    '''
DESCRIPTION

    TMalign wrapper

    Reference: Y. Zhang and J. Skolnick, Nucl. Acids Res. 2005 33, 2302-9
    http://zhanglab.ccmb.med.umich.edu/TM-align/

USAGE

    tmalign mobile, target [, args [, exe ]]

ARGUMENTS

    mobile, target = string: atom selections

    args = string: Extra arguments like -d0 5 -L 100

    exe = string: Path to TMalign executable {default: TMalign}

    ter = 0/1: If ter=0, then ignore chain breaks because TMalign will stop
    at first TER record {default: 0}

SEE ALSO

    tmscore, mmalign
    '''
    import subprocess
    import tempfile
    import os
    import re

    ter, quiet = int(ter), int(quiet)

    mobile_filename = tempfile.mktemp('.pdb', 'mobile')
    target_filename = tempfile.mktemp('.pdb', 'target')
    matrix_filename = tempfile.mktemp('.txt', 'matrix')
    mobile_ca_sele = '(%s) and (not hetatm) and name CA and alt +A' % (mobile)
    target_ca_sele = '(%s) and (not hetatm) and name CA and alt +A' % (target)

    if ter:
        save = cmd.save
    else:
        save = save_pdb_without_ter
    save(mobile_filename, mobile_ca_sele)
    save(target_filename, target_ca_sele)

    exe = cmd.exp_path(exe)
    args = [exe, mobile_filename, target_filename, '-m', matrix_filename] + args.split()

    try:
        process = subprocess.Popen(args, stdout=subprocess.PIPE,
                universal_newlines=True)
        lines = process.stdout.readlines()
    except OSError:
        print('Cannot execute "%s", please provide full path to TMscore or TMalign executable' % (exe))
        raise CmdException
    finally:
        os.remove(mobile_filename)
        os.remove(target_filename)

    # TMalign >= 2012/04/17
    if os.path.exists(matrix_filename):
        lines += open(matrix_filename).readlines()
        os.remove(matrix_filename)

    r = None
    re_score = re.compile(r'TM-score\s*=\s*(\d*\.\d*)')
    rowcount = 0
    matrix = []
    line_it = iter(lines)
    headercheck = False
    alignment = []
    for line in line_it:
        if 4 >= rowcount > 0:
            if rowcount >= 2:
                a = list(map(float, line.split()))
                matrix.extend(a[2:5])
                matrix.append(a[1])
            rowcount += 1
        elif not headercheck and line.startswith(' * '):
            a = line.split(None, 2)
            if len(a) == 3:
                headercheck = a[1]
        elif line.lower().startswith(' -------- rotation matrix'):
            rowcount = 1
        elif line.startswith('(":" denotes'):
            alignment = [next(line_it).rstrip() for i in range(3)]
        else:
            match = re_score.search(line)
            if match is not None:
                r = float(match.group(1))
        if not quiet:
            print(line.rstrip())

    if not quiet:
        for i in range(0, len(alignment[0]) - 1, 78):
            for line in alignment:
                print(line[i:i + 78])
            print('')

    assert len(matrix) == 3 * 4
    matrix.extend([0.0, 0.0, 0.0, 1.0])

    if int(transform):
        cmd.transform_selection('byobject (%s)' % (mobile), matrix, homogenous=1)

    # alignment object
    if object is not None:
        mobile_idx, target_idx = [], []
        space = {'mobile_idx': mobile_idx, 'target_idx': target_idx}
        cmd.iterate(mobile_ca_sele, 'mobile_idx.append("%s`%d" % (model, index))', space=space)
        cmd.iterate(target_ca_sele, 'target_idx.append("%s`%d" % (model, index))', space=space)
        for i, aa in enumerate(alignment[0]):
            if aa == '-':
                mobile_idx.insert(i, None)
        for i, aa in enumerate(alignment[2]):
            if aa == '-':
                target_idx.insert(i, None)
        if (len(mobile_idx) == len(target_idx) == len(alignment[2])):
            cmd.rms_cur(
                ' '.join(idx for (idx, m) in zip(mobile_idx, alignment[1]) if m in ':.'),
                ' '.join(idx for (idx, m) in zip(target_idx, alignment[1]) if m in ':.'),
                cycles=0, matchmaker=4, object=object)
        else:
            print('Could not load alignment object')

    if not quiet and r is not None:
        print('Found in output TM-score = %.4f' % (r))

    return r


def tmscore(mobile, target, args='', exe='TMscore', quiet=0, **kwargs):
    '''
DESCRIPTION

    TMscore wrapper

    Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710
    http://zhanglab.ccmb.med.umich.edu/TM-score/

ARGUMENTS

    mobile, target = string: atom selections

    args = string: Extra arguments like -d 5

    exe = string: Path to TMscore executable {default: TMscore}

    ter = 0/1: If ter=0, then ignore chain breaks because TMscore will stop
    at first TER record {default: 0}

SEE ALSO

    tmalign, mmalign
    '''
    kwargs.pop('_self', None)
    return tmalign(mobile, target, args, exe, quiet=quiet, **kwargs)

def save_pdb_without_ter(filename, selection, **kwargs):
    '''
DESCRIPTION

    Save PDB file without TER records. External applications like TMalign and
    DynDom stop reading PDB files at TER records, which might be undesired in
    case of missing loops.
    '''
    v = cmd.get_setting_boolean('pdb_use_ter_records')
    if v:
        cmd.set('pdb_use_ter_records', 0)
    cmd.save(filename, selection, **kwargs)
    if v:
        cmd.set('pdb_use_ter_records')
