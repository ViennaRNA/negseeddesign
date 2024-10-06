import time
START = time.time()

import sys
import argparse
from subprocess import Popen, PIPE
from pathlib import Path

import infrared as ir
import infrared.rna as rna

import RNA

sys.path.append(str(Path(__file__).parent/'linearbpdesign'))
from linearbpdesign.sampler import Sampler as LinearSampler

# in kcal/mol
ENERGYWEIGHT = - 1.98717 * (273.15 + 37) / 1000


def is_unique(seq):
    fc = RNA.fold_compound(seq)
    sub = fc.subopt(10)
    return len(sub)==1 or (sub[0].energy!=sub[1].energy)


def create_model_uniform(target):
    """Create infrared model from given target that supports uniform sequence sampling. The model pre-restricts nucleotides in paired positions
    """
    model = ir.Model(len(target), 4)
    model.add_constraints(rna.BPComp(i, j) for (i, j) in rna.parse(target))
    return model


def create_model_incarnation(target):
    """Create infrared model in incarnation way with targeted gc value if given
    Sequence weight is defined by basepair energy
    """
    model = ir.Model(len(target), 4)
    bps = rna.parse(target)
    model.add_constraints(rna.BPComp(i, j) for (i, j) in bps)

    # Add function
    model.add_functions([rna.BPEnergy(i, j, (i-1, j+1) not in bps) for (i, j) in bps], 'energy')
    model.add_functions([rna.GCCont(i) for i in range(len(target))], 'gc')

    # Set weight
    model.set_feature_weight(ENERGYWEIGHT, 'energy')
    model.set_feature_weight(0, 'gc')
    return model


def create_model_GC(target):
    """Create infrared model wigh random GC at unpaired region and A at unpaired positions
    """
    model = ir.Model(len(target), 4)
    model.add_constraints(rna.BPComp(i, j) for (i, j) in rna.parse(target))

    # Restrict domain
    for ind, c in enumerate(target):
        if c == '.':
            model.restrict_domains(ind, (0, 0))
        else:
            model.restrict_domains(ind, (1, 2))
    return model



def design(target, sampler, nSol=1, timeLimit=30):
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
        if (bp == 0) and is_unique(final):
            nFound += 1
            new_time = time.time()
            print(target, seed, final, hamming, int(bp), nRound, True, f'{new_time-current_time:.3f}', sep='\t')
            # Reinit
            nRound = 0
            current_best = (None, None, None, len(target))
            current_time = new_time
        elif bp < current_best[3]:
            current_best = (seed, final, hamming, int(bp))
        if time.time() - START >timeLimit:
             break

    if current_best[0] is not None:
        print(target, *current_best, nRound, False, f'{time.time() - current_time:.3f}', sep='\t')

class ModelSampler:
    """Simple sequence sampler from given infrared model
    """
    def __init__(self, model):
        self.sampler = ir.Sampler(model)
        # Force precomputation
        self.sampler.sample()

    def sample(self):
        return rna.ass_to_seq(self.sampler.sample())

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
    parser.add_argument('--seed', choices=['uniform', 'incarnation', 'gcheavy', 'linearbp'], default='uniform', help='Strategy to initiate seed sequences')
    parser.add_argument('--time', type=int, default=3600, help='Maximal running time in second')

    args = parser.parse_args()

    target = args.target
    model = None
    sampler = None
    match args.seed:
        case 'uniform':
            model = create_model_uniform(target)
        case 'incarnation':
            model = create_model_incarnation(target)
        case 'gcheavy':
            model = create_model_GC(target)
        case 'linearbp':
            sampler = LinearSampler(target)

    if model is not None:
        sampler = ModelSampler(model)

    design(target, sampler, nSol=args.number, timeLimit=args.time)
