"""
TAP CDR Assigner - antibody CDR region assignment supporting both heavy and light chains.
Based on TNP's CDR_Assigner but extended for full antibody (VH+VL) support.
"""
from scripts.region_definitions import Accept
from anarci import run_anarci

# Instantiate Accept
a = Accept()
a.numbering_scheme = "imgt"


def find_all_regions(output_name, numbering, Chain):
    """Assign CDR/framework regions for a numbered sequence. Supports H and L chains."""

    definitions = ['imgt', 'north', 'chothia', 'kabat']

    regions = {}
    cdr_length_dict = {}

    for definition in definitions:
        a.definition = definition
        regions[definition] = {}
        cdr_length_dict[definition] = {}

        if Chain == "H":
            list_of_regions = ["cdrh3", "fwh1", "cdrh1", "fwh2", "cdrh2", "fwh3", "fwh4"]
            list_of_cdrs = ["cdrh1", "cdrh2", "cdrh3"]
        elif Chain == "L":
            list_of_regions = ["cdrl3", "fwl1", "cdrl1", "fwl2", "cdrl2", "fwl3", "fwl4"]
            list_of_cdrs = ["cdrl1", "cdrl2", "cdrl3"]
        else:
            print('ERROR: Chain was not H or L')
            return {}, {}

        for r in list_of_regions:
            region = ''
            a.set_regions([r])

            for n in range(len(numbering)):
                if a.accept(numbering[n][0], Chain) == 1:
                    if numbering[n][1] == '-':
                        continue
                    else:
                        region += numbering[n][1]

            regions[definition][r] = region
            if r in list_of_cdrs:
                cdr_length_dict[definition][r] = len(regions[definition][r])

    return regions, cdr_length_dict


def anarci_fun(name, sequence, ncores=1):
    """Run ANARCI numbering on a single sequence."""
    anarci_tuple = (name, sequence)
    output = run_anarci([anarci_tuple], scheme='imgt', assign_germline=True,
                        allowed_species=None, ncpu=ncores)
    return output


def main(output_name, sequence, chain, output_dest=None, ncpu=1, verbose=True):
    """
    Main entry point for TAP CDR assignment.
    Returns (primary_seqs, length_dict) for the given chain.
    """
    if chain == "H":
        chain_id = "Heavy"
    elif chain == "L":
        chain_id = "Light"
    else:
        print('ERROR: chain type not H or L')
        return None, None

    output = anarci_fun(output_name, sequence, ncpu)

    assert output[1] != [None]

    numbering = output[1][0][0][0]

    primary_seqs, length_dict = find_all_regions(output_name, numbering, chain)

    return primary_seqs, length_dict
