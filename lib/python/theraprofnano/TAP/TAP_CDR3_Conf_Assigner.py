"""
TAP CDR3 Compactness Assigner - antibody CDR3 conformation analysis.
Extended from TNP's CDR3_Conf_Assigner to correctly handle antibody PDBs with
both heavy and light chains (filters to heavy chain only for CDR3 analysis).
"""
import numpy as np
import os
import sys
from sklearn.decomposition import PCA

from theraprofnano.Hydrophobicity_and_Charge_Profiler.Hydrophobicity_and_Charge_Assigner import (
    create_temp_folder
)

# CDR definitions for both numbering schemes (heavy chain only for compactness)
definitions = {
    "chothia": {
        "H1": ["H24", "H25", "H26", "H27", "H28", "H29", "H30", "H31", "H32", "H33", "H34"],
        "H2": ["H50", "H51", "H52", "H53", "H54", "H55", "H56", "H57", "H58"],
        "H3": ["H93", "H94", "H95", "H96", "H97", "H98", "H99", "H100", "H101", "H102", "H103", "H104"],
    },
    "imgt": {
        "H1": ["H25", "H26", "H27", "H28", "H29", "H30", "H31", "H32", "H33", "H34",
               "H35", "H36", "H37", "H38", "H39", "H40"],
        "H2": ["H54", "H55", "H56", "H57", "H58", "H59", "H60", "H61", "H62", "H63",
               "H64", "H65", "H66", "H67"],
        "H3": ["H103", "H104", "H105", "H106", "H107", "H108", "H109", "H110", "H111",
               "H112", "H113", "H114", "H115", "H116", "H117", "H118", "H119"],
    },
}


def parse_antibody(file, numbering_scheme):
    """Parse antibody PDB for CDR3 compactness. Filters to heavy chain only."""

    backbone_atoms = ['CA', 'C', 'N']

    if numbering_scheme == 'imgt':
        loop_boundaries = {1: list(range(27, 38 + 1)),
                           2: list(range(56, 65 + 1)),
                           3: list(range(105, 117 + 1))}
        anchors = {1: [24, 25, 26, 39, 40, 41],
                   2: [53, 54, 55, 66, 67, 68],
                   3: [102, 103, 118, 119]}

    if numbering_scheme == 'chothia':
        loop_boundaries = {1: list(range(26, 32 + 1)),
                           2: list(range(52, 56 + 1)),
                           3: list(range(95, 102 + 1))}
        anchors = {1: [23, 24, 25, 33, 34, 35],
                   2: [49, 50, 51, 57, 58, 59],
                   3: [92, 93, 103, 104]}

    data = {}
    data['chain'] = ''
    data['anchors'] = {1: [], 2: [], 3: []}
    data['cdrs'] = {1: [], 2: [], 3: []}

    with open(file, 'r') as s:
        s = s.readlines()
        for atom in s:
            atom = atom.strip('\n').split()
            if atom[0] == 'ATOM':
                # Only process heavy chain for CDR3 compactness
                if atom[4] != 'H':
                    continue
                data['chain'] = atom[4]
                if atom[2] in backbone_atoms:
                    for i in [1, 2, 3]:
                        resnum = int(''.join(filter(str.isdigit, atom[5])))
                        if resnum in loop_boundaries[i]:
                            data['cdrs'][i].append(np.array(
                                [float(atom[6]), float(atom[7]), float(atom[8])]))
                        elif resnum in anchors[i]:
                            data['anchors'][i].append(np.array(
                                [float(atom[6]), float(atom[7]), float(atom[8])]))

                if atom[2] == 'CA' and int(''.join(filter(str.isdigit, atom[5]))) == 118:
                    data['h3_anchor_ref'] = np.array(
                        [float(atom[6]), float(atom[7]), float(atom[8])])

                if atom[2] == 'CA' and int(''.join(filter(str.isdigit, atom[5]))) == 3:
                    data['orientation_ref'] = np.array(
                        [float(atom[6]), float(atom[7]), float(atom[8])])

    return data


def get_H3_anchor_line(data):
    H3_anchor = np.array(data['anchors'][3])
    pca = PCA(n_components=1)
    pca.fit(H3_anchor)

    centroid = np.mean(H3_anchor, axis=0)
    p1 = pca.components_[0]
    p2 = -1 * pca.components_[0]
    d1 = np.linalg.norm((p1 + centroid) - data['h3_anchor_ref'])
    d2 = np.linalg.norm((p2 + centroid) - data['h3_anchor_ref'])

    if d1 < d2:
        return centroid, p2
    else:
        return centroid, p1


def get_spherical_coordinates(h, c):
    rho = np.linalg.norm(h - c)
    return rho


def main_compactness(input_file, numbering_scheme, verbose=True):
    """
    Calculate CDR3 compactness for an antibody structure.
    Returns rho (loop reach). Compactness = CDR3_length / rho.
    Filters to heavy chain only.
    """

    temp_dir = create_temp_folder()
    os.mkdir(os.path.join(temp_dir, 'compactness'))

    if verbose:
        print("Temporary results stored in", temp_dir)

    temp_file = os.path.join(temp_dir, 'input.pdb')
    os.system('cp ' + input_file + ' ' + temp_file)

    try:
        data = parse_antibody(temp_file, numbering_scheme)
        c, __ = get_H3_anchor_line(data)
        h = np.mean(np.array(data['cdrs'][3]), axis=0)
        rho = get_spherical_coordinates(h, c)
        return rho
    except Exception as e:
        print(e)
        return None

    os.system('rm -rf ' + temp_dir)
