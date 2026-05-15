"""
TAP Hydrophobicity and Charge Assigner - antibody surface property calculation.
Extended from TNP's Hydrophobicity_and_Charge_Assigner to support both VH and VL chains.
"""
import numpy as np
import os
import sys

from theraprofnano.Hydrophobicity_and_Charge_Profiler.Hydrophobicity_and_Charge_Assigner import (
    stripInsertionCode, parsePSA, parsePSAArea, runPSA, create_temp_folder,
    normalize, residue_classes, AAtable, Naive_HH, coloring
)
from theraprofnano.Hydrophobicity_and_Charge_Profiler.Common.PDBUtils import PDBchain
from os.path import join

local_path = os.path.dirname(os.path.realpath(__file__))

# Extended CDR definitions with light chain support (+/- 2 on either side)
definitions = {
    "chothia": {
        "H1": ["H24", "H25", "H26", "H27", "H28", "H29", "H30", "H31", "H32", "H33", "H34"],
        "H2": ["H50", "H51", "H52", "H53", "H54", "H55", "H56", "H57", "H58"],
        "H3": ["H93", "H94", "H95", "H96", "H97", "H98", "H99", "H100", "H101", "H102", "H103", "H104"],
        "L1": ["L24", "L25", "L26", "L27", "L28", "L29", "L30", "L31", "L32", "L33", "L34"],
        "L2": ["L50", "L51", "L52", "L53", "L54", "L55", "L56"],
        "L3": ["L89", "L90", "L91", "L92", "L93", "L94", "L95", "L96", "L97"],
    },
    "imgt": {
        "H1": ["H25", "H26", "H27", "H28", "H29", "H30", "H31", "H32", "H33", "H34",
               "H35", "H36", "H37", "H38", "H39", "H40"],
        "H2": ["H54", "H55", "H56", "H57", "H58", "H59", "H60", "H61", "H62", "H63",
               "H64", "H65", "H66", "H67"],
        "H3": ["H103", "H104", "H105", "H106", "H107", "H108", "H109", "H110", "H111",
               "H112", "H113", "H114", "H115", "H116", "H117", "H118", "H119"],
        "L1": ["L24", "L25", "L26", "L27", "L28", "L29", "L30", "L31", "L32", "L33",
               "L34", "L35", "L36", "L37", "L38", "L39", "L40", "L41"],
        "L2": ["L53", "L54", "L55", "L56", "L57", "L58", "L59", "L60", "L61", "L62",
               "L63", "L64", "L65", "L66", "L67", "L68"],
        "L3": ["L103", "L104", "L105", "L106", "L107", "L108", "L109", "L110", "L111",
               "L112", "L113", "L114", "L115", "L116", "L117", "L118", "L119"],
    },
}


def is_CDR(res, deff):
    _res, ins = stripInsertionCode(res)
    for CDR in definitions[deff]:
        if _res in definitions[deff][CDR]:
            return CDR
    return False


def is_partic_CDR(res, deff, cdr_id):
    _res, ins = stripInsertionCode(res)
    if _res in definitions[deff][str(cdr_id)]:
        return cdr_id
    return False


def CreateAnnotation(annotation_index, which_ph, input_file, chains, numbering_scheme,
                     verbose=True):
    """
    Calculate PSH, PPC, PNC for an antibody structure (H+L chains).
    Extended from TNP version to process both heavy and light chain CDRs.
    """

    # Load hydrophobicity data from TNP's data file
    tnp_dat_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
        'Hydrophobicity_and_Charge_Profiler', 'dat', 'Hydrophobics.txt'
    )
    hydrophobics = {}
    for line in open(tnp_dat_path):
        line = line.strip().split('\t')
        hydrophobics[line[0].upper()] = (float(line[1]), float(line[2]), float(line[3]),
                                          float(line[4]), float(line[5]))

    hydrophobicity_annotation = normalize(hydrophobics, annotation_index)

    temp_dir = create_temp_folder()
    os.mkdir(join(temp_dir, 'surface'))

    if verbose:
        print("Temporary results stored in", temp_dir)

    chothia_file = join(temp_dir, 'input.pdb')
    os.system('cp ' + input_file + ' ' + chothia_file)

    # Run PSA
    surface_file = join(temp_dir, 'surface', 'result.psa')
    runPSA(chothia_file, surface_file)
    asa = parsePSA(surface_file, verbose=verbose)
    asa_area = parsePSAArea(surface_file)

    # Load both H and L chains
    nb_structure = {'H': PDBchain(chothia_file, 'H'),
                    'L': PDBchain(chothia_file, 'L')}
    process_chains = ['H', 'L']

    # Build residue dictionary from both chains
    residues = {}
    for chain in process_chains:
        for elem in sorted(nb_structure[chain].residues):
            typ = nb_structure[chain].residues[elem].res_type
            r_elem = {'res': nb_structure[chain].residues[elem], 'typ': typ, 'charge': 0}
            res_asa = 0
            res_asa_area = 0

            try:
                res_asa = asa[chain][str(elem[0]) + elem[1]]
                res_asa_area = asa_area[chain][str(elem[0]) + elem[1]]
            except KeyError:
                sys.stderr.write("ASA ERROR: %s\n" % input_file)

            if res_asa < 7.5:  # buried residue threshold
                continue

            r_elem['asa'] = res_asa
            r_elem['asa_sum'] = res_asa_area
            if np.isnan(r_elem['asa_sum']):
                if verbose:
                    print("NAN Residue Detected: " + r_elem['typ'], chain, elem[0] + elem[1])
                break
            r_elem['hydrophobic'] = float(hydrophobicity_annotation[typ])

            residues[(chain, elem)] = r_elem

    # Build CDR vicinity (across both chains)
    cdr_vicinity = []
    num_cdr_residues = 0

    for r1 in residues:
        numbering = r1[0] + str(r1[1][0])
        cdr_match = is_CDR(numbering, numbering_scheme)

        if cdr_match is False:
            # Non-CDR residue - assign charge
            if r1[0] in process_chains:
                aa = AAtable[str(residues[r1]['typ'])]
                if aa == "D" or aa == "E":
                    residues[r1]['negative'] = -1
                    residues[r1]['charge'] = -1
                elif aa == "R" or aa == "K":
                    residues[r1]['positive'] = 1
                    residues[r1]['charge'] = 1
                elif aa == "H":
                    residues[r1]['positive'] = 0.1
                    residues[r1]['charge'] = 0.1
            continue
        else:
            # CDR residue
            num_cdr_residues += 1
            cdr_vicinity.append(r1)

            if r1[0] in process_chains:
                aa = AAtable[str(residues[r1]['typ'])]
                if aa == "D" or aa == "E":
                    residues[r1]['negative'] = -1
                    residues[r1]['charge'] = -1
                elif aa == "R" or aa == "K":
                    residues[r1]['positive'] = 1
                    residues[r1]['charge'] = 1
                elif aa == "H":
                    residues[r1]['positive'] = 0.1
                    residues[r1]['charge'] = 0.1

        # Find non-CDR neighbors within 4A
        for r2 in residues:
            if r1 == r2 or r2 in cdr_vicinity:
                continue
            res1 = residues[r1]['res']
            res2 = residues[r2]['res']
            d = res2.distance(res1, res2)
            if d < 4.0:
                if r2 not in cdr_vicinity:
                    cdr_vicinity.append(r2)

    # Salt bridge neutralization
    sb_list = []
    num_saltbridges = 0
    for r1 in residues:
        for r2 in residues:
            if r1 == r2:
                continue
            res1 = residues[r1]['res']
            res2 = residues[r2]['res']
            donors = ["LYS", "ARG"]
            acceptors = ["ASP", "GLU"]

            if res1.res_type in donors and res2.res_type in acceptors:
                d = res1.sb_distance(res1, res2)
                if d < 3.2:
                    residues[r1]['Positive'] = 0
                    residues[r1]['Charge'] = 0
                    residues[r1]['hydrophobic'] = hydrophobicity_annotation["GLY"]
                    residues[r2]['Negative'] = 0
                    residues[r2]['Charge'] = 0
                    residues[r2]['hydrophobic'] = hydrophobicity_annotation["GLY"]
                    sb_list.append((
                        str(r1[0]) + str(r1[1][0]) + str(r1[1][1]) + ": " + str(res1.res_type),
                        str(r2[0]) + str(r2[1][0]) + str(r2[1][1]) + ": " + str(res2.res_type),
                        str("{0:.4f}".format(d))
                    ))
                    num_saltbridges += 1

    # Adjacency matrix for CDR vicinity (both chains)
    adj_matrix_cdr = {}
    for i in residue_classes:
        adj_matrix_cdr[i] = {}
        for j in residue_classes:
            adj_matrix_cdr[i][j] = 0

    for r1 in cdr_vicinity:
        for r2 in cdr_vicinity:
            if r1 != r2:
                res1 = residues[r1]['res']
                res2 = residues[r2]['res']
                d = res2.distance(res1, res2)

                if d > 7.5:
                    continue

                for c1 in residue_classes:
                    for c2 in residue_classes:
                        if c1 == c2:
                            if c1 in residues[r1] and c2 in residues[r2]:
                                coulombic_score = (residues[r1][c1] * residues[r2][c2]) / (d * d)
                                adj_matrix_cdr[c1][c2] += coulombic_score

    stats = {}
    stats[annotation_index] = {}
    stats[annotation_index]['Patch_Hydrophob_CDR'] = float(
        "{0:.4f}".format(adj_matrix_cdr['hydrophobic']['hydrophobic']))
    stats[annotation_index]['Patch_Pos_Charge_CDR'] = float(
        "{0:.4f}".format(adj_matrix_cdr['positive']['positive']))
    stats[annotation_index]['Patch_Neg_Charge_CDR'] = float(
        "{0:.4f}".format(adj_matrix_cdr['negative']['negative']))

    os.system('rm -rf ' + temp_dir)

    return stats
