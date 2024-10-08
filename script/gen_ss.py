"""Script to generate secondary structure of required length
The structure is uniformly generated or be MFE of uniformly sampled sequences
Corresponds to sec 6 and sec 7.2 in manuscript
"""

import sys
import argparse
from pathlib import Path, PurePath

import RNA


sys.path.append(str(Path(__file__).parent.parent/'linearbpdesign'))
sys.path.append(str(Path(__file__).parent.parent))
import designBP as design


def min_helix_length(ss):
    """Get minimum helix length
    """
    target_tree = design.dbn_to_tree(design.ssparse(ss))
    lst = design.count_helix(target_tree)
    if len(lst) >= 1:
        return min(lst)
    return 0


def sample_uniform(size, nb, helix_len=1):
    """Uniformly sample nb secondary structure of given size with desired helix len
    """
    count, count_stacked = {}, {}
    sampled = 0
    while sampled < nb:
        ss = design.ssrandom(size, count, count_stacked, 3, helix_len)
        min_len = min_helix_length(ss)
        if min_len == helix_len == 1:
            print(ss)
            sampled += 1
        elif min_len == helix_len == 2:
            print(ss)
            sampled += 1
        elif min_len >= helix_len == 3:
            print(ss)
            sampled += 1



def sample_mfe(size, nb, verbose=False):
    """Return MFE structure of randomly generated sequences
    """
    for i in range(nb):
        w = RNA.random_string(size, "ACGU")
        fc = RNA.fold_compound(w)
        mfe, _ = fc.mfe()
        fc.pf()
        ens = fc.ensemble_defect(mfe)
        if verbose:
            print(mfe, ens, (w.count('G') + w.count('C'))/size, w, sep=',')
        else:
            print(mfe, ens, sep=',')



if __name__ == "__main__":
    # General options
    general_group = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Script to generate secondary structure of required length
        The structure is uniformly generated or be MFE of uniformly sampled sequences
        Corresponds to sec 6 and sec 7.2 in manuscript
        """
    )

    # Common argument
    general_group.add_argument('size', metavar='Size', type=int, help="size of structure to sample")
    general_group.add_argument('-n', type=int, default=5000, help="Number of sampled structures (default: %(default)d)")
    general_group.add_argument('--helix_length', type=int, default=1, help="Helix length used in uniform sampling (default: %(default)d)")
    general_group.add_argument('--mfe', action="store_true", help="Use MFE structures from sampled sequence")
    general_group.add_argument('--verbose', action="store_true", help="Enable verbose mode")

    args = general_group.parse_args()

    if args.mfe:
        sample_mfe(args.size, args.n, verbose=args.verbose)
    else:
        sample_uniform(args.size, args.n, helix_len=args.helix_length)
