"""Simple script to compute information for given structure
"""

import sys
import argparse
from pathlib import Path
import time


sys.path.append(str(Path(__file__).parent.parent/'linearbpdesign'))
sys.path.append(str(Path(__file__).parent.parent))
sys.path.append(str(Path(__file__).parent.parent/'script'))
import designBP as design
from sampler import Sampler as LinearSampler
from samplerbiseparable import Sampler as BILinearSampler

def min_helix_length_of_ss(ss):
    """Get minimum helix length
    """
    return min_helix_of_tree(design.dbn_to_tree(design.ssparse(ss)))

def min_helix_length_of_tree(target_tree):
    """Get minimum helix length
    """
    lst = design.count_helix(target_tree)
    if len(lst) >= 1:
        return min(lst)
    return 0

def deg_of_tree(v):
    """Compute max multiloop degree for given tree
    """
    BP_children = [vv[0] for vv in v[1] if not design.is_leaf(vv)]
    resu = len(BP_children) + int(v[0] != "root")
    for vv in v[1]:
        resu = max(resu, deg_of_tree(vv))
    return resu

def deg_of_ss(ss):
    """Compute max multiloop degree for structure in dbn
    """
    return deg_of_tree(design.dbn_to_tree(design.ssparse(ss)))

def run_ss_bi(ss, withNeg, skipDeg, includeTime, minSol, modulo):
    """Get structure information for biseparability
    """
    if includeTime:
        start = time.time()
    deg = deg_of_ss(ss)
    precompute = withNeg and deg<=skipDeg
    upto = None if modulo <= 2 else modulo
    inst = BILinearSampler(ss, precompute=precompute, uptomoduloA=upto, uptomoduloC=upto)
    tree = inst.target_tree
    helix = min_helix_length_of_tree(tree)
    helix = min(3, helix)
    # deg = deg_of_tree(tree)
    isProper = design.filter(tree)

    res = [ss, helix, deg, isProper]
    nb_sol = None
    if precompute:
        min_moduloA = inst.min_moduloA if upto is None else max(inst.min_moduloA, upto)
        min_moduloC = inst.min_moduloC if upto is None else max(inst.min_moduloC, upto)
        if min_moduloA is None or min_moduloC is None:
            nb_sol = None
        # elif upto is not None:
        #     min_moduloA = max(upto, min_moduloA)
        #     min_moduloC = max(upto, min_moduloC)
        #     nb_sol = sum(inst.weights)
        # else:
            # nb_sol = inst.solution_number(min_moduloA, min_moduloC)
        if sum(inst.weights) < minSol:
            min_moduloA, min_moduloC, nb_sol = inst.auto_fill(fromModulo=upto, minSol=minSol)

        res += [min_moduloA, min_moduloC]
        nb_sol = sum(inst.weights)
        res.append(nb_sol)

    if includeTime:
        res.append(time.time()-start)

    print(*res, sep='\t')

def run_ss(ss, withNeg, skipDeg, includeTime, minSol, modulo):
    if includeTime:
        start = time.time()
    deg = deg_of_ss(ss)
    precompute = withNeg and deg<=skipDeg
    inst = LinearSampler(ss, precompute=precompute)
    # inst = BILinearSampler(ss, precompute=precompute)
    tree = inst.target_tree
    helix = min_helix_length_of_tree(tree)
    isProper = design.filter(tree)
    upto = None if modulo <= 2 else modulo

    res = [ss, helix, deg, isProper]
    min_modulo = None
    nb_sol = None
    if precompute:
        min_modulo = inst.min_modulo
        if min_modulo is None:
            nb_sol = None
        else:
            # nb_sol = inst.solution_number(min_modulo)
            min_modulo, nb_sol = inst.auto_fill(fromModulo=upto, minSol=minSol)

        res.append(min_modulo)
        res.append(nb_sol)

    if includeTime:
        res.append(time.time()-start)

    print(*res, sep='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description =
        """Script to give simple structure information
        """
    )
    parser.add_argument('--withNeg', action='store_true', help='Add linearBPdesign related')
    parser.add_argument('--minSol', type=int, default=0, help='Return if possilbe min modulo such that the solution counts is at least minSol. However it stops if modulo goes too high')
    parser.add_argument('-m', '--modulo', type=int, default=0, help='Count solution number up to modulo')
    parser.add_argument('--includeTime', action='store_true', help='Add running time')
    parser.add_argument('--onlyA', action='store_true', help='Use A-separable')
    parser.add_argument('--skipDeg', type=int, default=20, help='Max multi-loop degree to compute')

    args = parser.parse_args()
    HEADER = ["Structure", "Min Helix", "Max Multi Deg", "Proper"]
    if args.withNeg:
        if args.onlyA:
            HEADER += ["Min Modulo", "Nb Sol"]
        else:
            HEADER += ["Min ModuloA", "Min ModuloC", "Nb Sol"]
    if args.includeTime:
        HEADER += ["Time"]
    print(*HEADER, sep='\t')
    for line in sys.stdin.readlines():
        ss = line.strip()
        if args.onlyA:
            run_ss(ss, args.withNeg, args.skipDeg, args.includeTime, args.minSol, args.modulo)
        else:
            run_ss_bi(ss, args.withNeg, args.skipDeg, args.includeTime, args.minSol, args.modulo)
