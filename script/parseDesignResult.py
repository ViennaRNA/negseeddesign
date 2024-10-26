"""The script read thourgh a folder containing output of rundesign.py and put into one DataFrame
"""

import sys
import argparse
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd


import RNA

def compute_entropy(labels):
    """Positional Shannon entropy for given information
    """
    n_labels = len(labels)
    values, counts = np.unique(labels, return_counts=True)
    probs = counts / n_labels
    n_classes = np.count_nonzero(probs)
    if n_classes <= 1:
        return 0
    ent = 0.
    for i in probs:
        ent -= i * np.log2(i)
    return ent

def thermo_result(ss, seq):
    """Return energy, probability, and normalised ensemble defect for given pair of structure and sequence
    """
    fc = RNA.fold_compound(seq)
    _, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    fc.pf()
    e = fc.eval_structure(ss)
    p = fc.pr_structure(ss)
    ens = fc.ensemble_defect(ss)
    return [e, p, ens]


def pairwise_hamming(sequences):
    """Return average pair waise hamming distance
    """
    l = len(sequences)
    m = np.zeros((l, l))
    for i in range(l):
        for j in range(i+1, l):
            m[i][j] = sum([x!=y for x, y in zip(sequences[i], sequences[j])])
    return np.mean(m + m.T, axis=1)



def gc_count(ss, seq):
    paired = 0
    unpaired = 0
    gc_paired = 0
    gc_unpaired = 0
    for x, y in zip(ss, seq):
        if x == '.':
            unpaired += 1
            gc_unpaired += y in 'GC'
        else:
            paired += 1
            gc_paired += y in 'GC'
    return [gc_paired/paired, gc_unpaired/unpaired]


def summary_from_file(path):
    """Convert a design result output to a summary dataframe
    """
    names = ['ss', 'seq', 'mutation', 'hamming', 'rounds', 'CG paired', 'GC unpaired', 'energy', 'prob', 'ensemble', 'avg Hamming']
    res = []
    seedLst = []
    finalLst = []
    with Path(path).open() as f:
        for line in f.readlines():
            if line.startswith("slurmstepd: error:"):
                return None
            lst = line.split('\t')
            # We skip the line if the solution is not unique
            if lst[-2] == 'False':
                continue
            ss = lst[0]
            seed = lst[1]
            final = lst[2]
            # Fill thermo result for seed and final sequences
            res.append([ss, seed, 'before', None, None] +  gc_count(ss, seed) + thermo_result(ss, seed))
            res.append([ss, final, 'after', int(lst[3]), int(lst[5])] + gc_count(ss, final) + thermo_result(ss, final))
            seedLst.append(seed)
            finalLst.append(final)

    if len(res) == 0:
        return None

    for i, v in enumerate(pairwise_hamming(seedLst)):
        res[i*2].append(v)
    for i, v in enumerate(pairwise_hamming(finalLst)):
        res[i*2+1].append(v)

    df = pd.DataFrame(res, columns=names)
    entropy_seed = entropy_structures(ss, seedLst)
    entropy_final = entropy_structures(ss, finalLst)
    for x in ['Paired Entropy', 'Paired Entropy Di', 'Unpaired Entropy', 'Unpaired Entropy Di']:
        df.loc[df['mutation']=='before', x] = entropy_seed[x]
        df.loc[df['mutation']=='after', x] = entropy_final[x]

    return df

def summary_from_file_immediate(path):
    """Convert a immediate experiment result output to a summary dataframe
    """
    names = ['ss', 'seq', 'unique', 'CG paired', 'GC unpaired', 'energy', 'prob', 'ensemble']
    res = []
    cur_res = []
    cur_seedLst = []
    cur_ss = None
    cur_ind = 1
    path = Path(path)
    base = int(path.stem.split('_')[-1])
    print(f'Start {base}:')
    with path.open() as f:
        for line in f.readlines():
            if line.startswith("slurmstepd: error:"):
                return None
            lst = line.strip().split('\t')
            ss = lst[0]
            if cur_ss is None:
                cur_ss = ss
            # Dump result
            if ss != cur_ss:
                df = pd.DataFrame(cur_res, columns=names)
                entropy_seed = entropy_structures(cur_ss, cur_seedLst)
                for x in ['Paired Entropy', 'Paired Entropy Di', 'Unpaired Entropy', 'Unpaired Entropy Di']:
                    df.loc[:, x] = entropy_seed[x]
                df.loc[:, 'ind'] = (base-1)*10 + cur_ind
                res.append(df)

                cur_ss = ss
                cur_seedLst = []
                cur_res = []
                cur_ind += 1

            seed = lst[1]
            unique = True if lst[-1] == 'True' else False
            # Fill thermo result for seed and final sequences
            cur_res.append([ss, seed, unique] + gc_count(ss, seed) + thermo_result(ss, seed))
            cur_seedLst.append(seed)

        df = pd.DataFrame(cur_res, columns=names)
        entropy_seed = entropy_structures(cur_ss, cur_seedLst)
        for x in ['Paired Entropy', 'Paired Entropy Di', 'Unpaired Entropy', 'Unpaired Entropy Di']:
            df.loc[:, x] = entropy_seed[x]
        res.append(df)

    if len(res) == 0:
        return None

    return pd.concat(res, ignore_index=True)

def entropy_structures(ss, sequences):
    """Compute paired/unpaired entropy with mono or di nucleotides
    """
    paired_entropy = []
    unpaired_entropy = []
    paired_entropy_di = []
    unpaired_entropy_di = []
    for ind, c in enumerate(ss):
        labels = [x[ind] for x in sequences]
        ent = compute_entropy(labels)
        if c == '.':
            unpaired_entropy.append(ent)
        else:
            paired_entropy.append(ent)
        if ind == len(ss) - 1:
            continue
        di_labels = [x[ind]+x[ind+1] for x in sequences]
        di_ent = compute_entropy(di_labels)
        if c == '.':
            unpaired_entropy_di.append(di_ent)
        else:
            paired_entropy_di.append(di_ent)
    return {'Paired Entropy': np.mean(paired_entropy),
            'Paired Entropy Di': np.mean(paired_entropy_di),
            'Unpaired Entropy': np.mean(unpaired_entropy),
            'Unpaired Entropy Di': np.mean(unpaired_entropy_di)}




if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description =
        """The script read thourgh a folder containing output of rundesign.py and put into one DataFrame
        """
    )
    parser.add_argument('target', metavar='target', type=str, help='Target folder to parse')
    parser.add_argument('output', metavar='output', type=str, help='Where to store dataframe in pkl')
    parser.add_argument('-p', '--pool', type=int, default=1, help='Number of multithreading')
    # parser.add_argument('--seed', choices=['uniform', 'bpenergy', 'gcheavy', 'linearbp'], default='uniform', help='Strategy to initiate seed sequences')
    # # parser.add_argument('--time', type=int, default=3600, help='Maximal running time in second')
    # parser.add_argument('--print-any', action='store_true', help='Print any MFE and near result')
    parser.add_argument('--immediate', action='store_true', help='Parse results of immediate exp')
    # parser.add_argument('--insertC', choices=[0, 1, 50], type=int, default=0, help='Insert C when only A flag is set. Add C in none (0), one (1), or half of (50) loops for uniform and bpenergy seed.')
    # parser.add_argument('-m', '--modulo', type=int, default=0, help='Modulo')
    # parser.add_argument('--candidates', type=int, default=0, help='Increase LinearBP sampler modulo if possible such that candidates number is enough')

    args = parser.parse_args()

    folder = Path(args.target)

    parser = summary_from_file_immediate if args.immediate else summary_from_file
    res = []
    if args.pool == 1:
        for x in folder.iterdir():
            if not x.suffix == '.csv':
                continue
            ind = int(x.stem.split('_')[-1])

            if args.immediate:
                df = summary_from_file_immediate(x)
            else:
                print(f'{ind:05d}', end='\r')
                df = summary_from_file(x)
                if df is not None:
                    df.loc[:, 'puzzle'] = ind
            if df is not None:
                res.append(df)
    else:
        with Pool(args.pool) as pool:
            for df in pool.imap_unordered(parser, (x for x in folder.iterdir() if x.suffix == '.csv')):
                if df is not None:
                    res.append(df)

    print(len(res))
    # print(df)
    # print(df.iloc[0])
    # assert False
    summaryDf = pd.concat(res, ignore_index=True)
    summaryDf.to_pickle(args.output)
