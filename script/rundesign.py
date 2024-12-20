import time
START = time.time()

import sys
import argparse
from subprocess import Popen, PIPE
from pathlib import Path
import random
from math import exp

# import infrared as ir
# import infrared.rna as rna

import RNA

sys.path.append(str(Path(__file__).parent.parent/'linearbpdesign'))
sys.path.append(str(Path(__file__).parent.parent/'script'))
sys.path.append(str(Path(__file__).parent.parent))
from sampler import Sampler as LinearSampler
from samplerbiseparable import Sampler as BILinearSampler

# in kcal/mol
ENERGYWEIGHT = - 1.98717 * (273.15 + 37) / 1000
# BP threshold to print solution
NEAR = 10

## Default parameters for the base pair model (magic params from the Redprint paper)
BPLST = ['AU', 'CG', 'GC', 'GU', 'UA', 'UG']
params_bp_in = [exp(x/ENERGYWEIGHT) for x in [-0.52309, -2.10208, -2.10208, -0.88474, -0.52309, -0.88474]]
params_bp_term = [exp(x/ENERGYWEIGHT) for x in [1.26630, -0.09070, -0.09070, 0.78566, 1.26630, 0.78566]]



def dbn_to_bps(ss):
    """Return list of bps (0-index) of given nested basepairs
    """
    tmp = []
    res = []
    for ind, c in enumerate(ss):
        if c == '(':
            tmp.append(ind)
        elif c == ')':
            i = tmp.pop()
            res.append((i, ind))
    return res

def loops_from_dbn(ss):
    """Return list of loops for given dbn each loop is a list of unpaired positions
    """
    tmp = [[]]
    res = []
    for ind, c in enumerate(ss):
        if c == '(':
            tmp.append([])
        elif c == ')':
            res.append(tmp.pop())
        else:
            tmp[-1].append(ind)
    res.append(tmp.pop())
    return [x for x in res[-1::-1] if len(x)>0]

def is_unique(seq):
    fc = RNA.fold_compound(seq)
    sub = fc.subopt(1)
    return len(sub)==1 or (sub[0].energy!=sub[1].energy)


def create_sampler(target, strategy, onlyA, modulo=None, minCandidate=0, insertC=False):
    """Create different sampler
    """
    match strategy:
        case 'uniform':
            sampler = UniformSampler(target, onlyA)
        case 'bpenergy':
            sampler = BPEnergySampler(target, onlyA, insertC=insertC)
        case 'gcheavy':
            sampler = GCSampler(target)
        case 'linearbp':
            m = None if modulo == 0 else modulo
            # Fill cache upto modulo
            if onlyA:
                sampler = LinearSampler(target, uptomodulo=m)
            else:
                sampler = BILinearSampler(target, uptomoduloA=m, uptomoduloC=m)
            # Keep filling cache if not enough samplers in pool
            if sum(sampler.weights) < minCandidate:
                sampler.auto_fill(fromModulo=m, minSol=minCandidate)
    return sampler


def design(target, sampler, nSol=1, print_any=False):
    """Run RNAinverse and report result
    """
    nFound = 0
    nRound = 0
    current_best = (None, None, None, len(target))
    current_time = time.time()
    while nFound < nSol:
        nRound += 1
        seed = sampler.sample()
        # print(nRound, seed, time.time() - START, end='\r')
        final, bp = RNA.inverse_fold(seed, target)
        hamming = sum(x!=y for x, y in zip(seed, final))
        new_time = time.time()
        print_time = new_time - current_time
        if (bp == 0) and is_unique(final):
            nFound += 1
            print(target, seed, final, hamming, int(bp), nRound, True, f'{print_time:.3f}', sep='\t')
            # Reinit
            nRound = 0
            current_best = (seed, final, hamming, 0)
            current_time = new_time
        elif bp < current_best[3]:
            current_best = (seed, final, hamming, int(bp))
            if (0 < current_best[3] <= NEAR) and print_any:
                print(target, seed, final, hamming, int(bp), nRound, False, f'{print_time:.3f}', sep='\t')
        if print_any and (bp == 0) and (not is_unique(final)):
            print(target, seed, final, hamming, int(bp), nRound, False, f'{print_time:.3f}', sep='\t')
        # if time.time() - START >timeLimit:
        #      break

    # if current_best[0] is not None:
    #     print(target, *current_best, nRound, False, f'{time.time() - current_time:.3f}', sep='\t')


class GCSampler:
    """Simple sampler for GC heavy seed w/o using Infrared
    """
    def __init__(self, target):
        self.target = target
        self.length = len(target)
        self.bps = dbn_to_bps(target)

    def sample(self):
        seq = ['A'] * self.length
        for i, j in self.bps:
            x = random.choice(['CG', 'GC'])
            seq[i] = x[0]
            seq[j] = x[1]
        return ''.join(seq)

class UniformSampler:
    """Simple uniform seed sampler w/o using Infrared
    """
    def __init__(self, target, withA=False):
        self.target = target
        self.length = len(target)
        self.bps = dbn_to_bps(target)
        self.withA = withA

    def sample(self):
        if self.withA:
            seq = ['A'] * self.length
        else:
            seq = [x for x in RNA.random_string(self.length, 'ACGU')]
        for i, j in self.bps:
            x = random.choice(BPLST)
            seq[i] = x[0]
            seq[j] = x[1]
        return ''.join(seq)

class BPEnergySampler:
    """Simple seed sampler based on base pair energy model w/o using Infrared
    """
    def __init__(self, target, withA=False, insertC=False):
        self.target = target
        self.length = len(target)
        self.bps = dbn_to_bps(target)
        self.bps_in = [(i, j) for i, j in self.bps if (i-1, j+1) in self.bps]
        self.bps_term = [(i, j) for i, j in self.bps if not (i-1, j+1) in self.bps]
        self.loops = loops_from_dbn(target)
        self.withA = withA
        self.insertC = insertC

    def sample(self):
        if self.withA:
            seq = ['A'] * self.length
            if self.insertC:
                for lst in self.loops:
                    x = random.choice('AC')
                    for t in lst:
                        seq[t] = x
        else:
            seq = [x for x in RNA.random_string(self.length, 'ACGU')]
        for i, j in self.bps_in:
            x = random.choices(BPLST, weights=params_bp_in)[0]
            seq[i] = x[0]
            seq[j] = x[1]
        for i, j in self.bps_term:
            x = random.choices(BPLST, weights=params_bp_term)[0]
            seq[i] = x[0]
            seq[j] = x[1]
        return ''.join(seq)

# class ModelSampler:
#     """Simple sequence sampler from given infrared model
#     """
#     def __init__(self, model):
#         self.sampler = ir.Sampler(model)
#         # Force precomputation
#         self.sampler.sample()

#     def sample(self):
#         return rna.ass_to_seq(self.sampler.sample())

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description =
        """Run RNAinverse with defined initial sequence strategy for given structure.
        The script calls RNAlib inverse_fold function repeatedly and reports the number of tries needed to find the solution.
        This is equivalent to running RNAinverse with negative value for -R option (assuming seed sequence is changed at each try)
        """
    )
    parser.add_argument('target', metavar='target', type=str, help='Target PK-free secondary structure in dbn')
    parser.add_argument('-n', '--number', type=int, default=1, help='(Maximum) number of solutions. The script will try to find up to n solutions within time limit')
    parser.add_argument('--seed', choices=['uniform', 'bpenergy', 'gcheavy', 'linearbp'], default='uniform', help='Strategy to initiate seed sequences')
    # parser.add_argument('--time', type=int, default=3600, help='Maximal running time in second')
    parser.add_argument('--print-any', action='store_true', help='Print any MFE and near result')
    parser.add_argument('--onlyA', action='store_true', help='Only A in unpaired region for Uniform, bpenergy, LinearBPDesign')
    parser.add_argument('-m', '--modulo', type=int, default=0, help='Modulo')
    parser.add_argument('--insertC', action='store_true', help='For each loop, fill with C in 50/50')
    parser.add_argument('--candidates', type=int, default=1000, help='Increase LinearBP sampler modulo if possible such that candidates number is enough')

    args = parser.parse_args()

    target = args.target
    sampler = create_sampler(target, args.seed, args.onlyA, modulo=args.modulo, minCandidate=args.candidates, insertC=args.insertC)

    # for _ in range(10):
    #     print(sampler.sample())
    # assert False


    design(target, sampler, nSol=args.number, print_any=args.print_any)
