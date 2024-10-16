"""Script to find immediate solutions of structures
"""

import sys
import argparse
from pathlib import Path
from multiprocessing import Pool

import RNA

sys.path.append(str(Path(__file__).parent.parent/'linearbpdesign'))
sys.path.append(str(Path(__file__).parent.parent/'script'))
sys.path.append(str(Path(__file__).parent.parent))
import designBP as design
from rundesign import create_sampler


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


def run_single(ss, nb, seed, onlyA=False, upto=None):
    """Find immediate solutions of one given straucture
    """
    inst = create_sampler(ss, seed, onlyA, modulo=upto, minCandidate=nb)
    for ind in range(1, nb+1):
        res = inst.sample()
        sol = res
        fc = RNA.fold_compound(sol)
        mfe_ss, mfe = fc.mfe()
        isSol = (mfe_ss == ss)
        # Check if solution is unique
        if isSol:
            sub = fc.subopt(1)
            isUnique = (len(sub)==1 or (sub[0].energy!=sub[1].energy))
        else:
            isUnique = False
        lst = [ss, sol, isSol, isUnique]
        flush = ind // 1000 == 0
        print(*lst, sep='\t', flush=flush)





if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description =
        """Script to find immediate solutions of structures
        If linearbp seed is chosen, the code will automatically increase modulo if possible such that solution count is greater than sample number
        """
    )
    parser.add_argument('input', metavar='input', type=str, help='Path to structures file')
    parser.add_argument('--seed', choices=['uniform', 'bpenergy', 'gcheavy', 'linearbp'], default='linearbp', help='Strategy to initiate seed sequences')
    parser.add_argument('-n', '--nbSol', type=int, default=10, help='Number of solutions for each structure')
    parser.add_argument('-s', '--start', type=int, default=1, help='Strat index of structure')
    parser.add_argument('-e', '--end', type=int, default=10, help='End index of structure')
    parser.add_argument('--upto', type=int, default=0, help='Strat index of structure')
    parser.add_argument('--withC', action='store_true', help='Mix with C')

    args = parser.parse_args()
    seed = args.seed
    nbSol = args.nbSol
    input_file = args.input
    start = args.start
    end = args.end
    upto = args.upto

    # Extract structures from start to end
    with open(input_file) as f:
        structures = [line.strip() for ind, line in enumerate(f.readlines()) if (start-1)<=ind<end]


    for ss in structures:
        run_single(ss, nbSol, seed, not args.withC, upto=upto)
