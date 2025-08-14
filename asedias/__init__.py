
from asedias.core import aseDIAS, Fragment
from ase.optimize import BFGS

__version__ = '0.0.1'

class Parameters:
    def __init__(self):
        # distortion interaction settings
        self.optimizer = BFGS
        self.preopt_fmax = 0.05
        self.preopt_steps = 100
        self.opt_fmax = 0.05
        self.opt_steps = 30

        # plot settings
        self.marker = True
        self.linestyle = True
        self.hline = True
        self.include_fragments = True
        self.yaxis_unit = "eV"
        self.relative_idx = 0

    def __str__(self):
        max_key_len = max(len(k) for k in self.__dict__)
        lines = []
        for k, v in self.__dict__.items():
            lines.append(f"{k:<{max_key_len}} : {v}")
        return "\n".join(lines)


ParameterManager = Parameters()