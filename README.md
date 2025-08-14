
# asedias

`asedias` is designed to streamline the activation strain analysis workflow and to provide a flexible calculator interface integrated with ASE (Atomic Simulation Environment). It also includes built-in plotting functionality to automate complex post-processing tasks.

| 💬 Conversational AI Docs  | 🔗 Try it on Google Colab |
|---------------------------|---------------------------|
| [![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/kangmg/asedias) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/kangmg/asedias/blob/main/notebook/asedias_examples.ipynb) |


## Installation

You can install the `asedias` package using either of the following methods:

1. Install directly from GitHub using pip:

```bash
pip install git+https://github.com/kangmg/asedias.git
```

2. Clone the repository and install locally:

```bash
git clone https://github.com/kangmg/asedias.git
cd asedias
pip install .
cd ..
rm -r asedias
```

## Example Usage

```python
from aimnet.calculators import AIMNet2ASE
from asedias.core import Fragment, aseDIAS
from ase.io import read

images = read('rxn.traj', index=':')

def attach_aimnet2(atoms, spin, charge):
    """Attach a calculator to the atoms."""
    calc = AIMNet2ASE('aimnet2nse', charge=charge, mult=(spin*2+1))
    atoms.calc = calc

# Define fragments as Fragment objects
frag1 = Fragment(index=[2, 7], charge=0, spin=0, frag_name='HF')
frag2 = Fragment(index=[0, 1, 3, 4, 5, 6], charge=0, spin=0, frag_name='ethane')

analysis = aseDIAS(images=images, fragments=[frag1, frag2], calc_attach=attach_aimnet2)
analysis.run()
plot = analysis.plot()
```
