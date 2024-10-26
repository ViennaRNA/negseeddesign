This project contains the script to run the experiments and reproduce the figures in the manuscript
__Old dog, new tricks: Exact seeding strategy improves RNA design performances__

# Dependencies

The project is heavily based on `ViennaRNA (RNAlib)` and [LinearBPDesign](https://gitlab.inria.fr/amibio/linearbpdesign). The latter one is installed as a git submodule. Run `git submodule update --init --recursive` to initiate.

Other required dependencies for figures are `numpy`, `pandas`, `seaborn`, `matplotlib`, `varnaapi`, `logomaker` ...


# Script

The folder `script` contains different scripts to run experiments described in the manuscript

- `rundesign.py`: design given target structure using `RNAinverse` with different seeds. For example, the command below returns 10 solutions of structure `..(((((..((.((((.......)))).))..)))))...` with Biseparable seed with modulo up to 3
	```bash
	python script/rundesign.py -n 10 --seed linearbp -m 3 "..(((((..((.((((.......)))).))..)))))..."
	```
- `immediate_solutions.py`: sample seeds and check whether they are T-design (w/o `RNAinverse`). The command below sample each 20 Boltzmann sampled seeds for 1st to 10th structures in xxx.txt (one line for each structure in dot-bracket notation)
	```bash
	python script/immediate_solutions.py xxx.txt -n 20 -s 1 -e 10 --seed bpenergy
	```
- `gen_ss.py`: generate uniform or MFE structures. The command below generates uniformly 1,000 structures of size 100 nts with helix length of 3+
	```bash
	python script/gen_ss.py 100 -n 1000 --helix_length 3
	```
- `parseDesignResult.py`: evaluate (ensemble defect, diversity ...) resulting design produced by `rundesign.py` and stored as pandas DataFrame in pickle. The script assumes the result of each target indiced i is stored in `puzzle_i.csv`. The command below parses all `puzzle_*.csv` results under `xxx` and stores the DataFrame in `yyy.pkl`.
	```bash
	python script/parseDesignResult.py xxx yyy.pkl
	```

# Notebooks

The folder `notebooks` contains jupyter notebooks as indicated individually by the file name to reproduce figures in the manuscript
