import ase
from ase.optimize import BFGS
from ase.constraints import FixAtoms
from typing import Callable, Union, Optional
import time
import uuid
from typing import Callable, Union, Optional
import ase
from ase.optimize import BFGS
from ase.constraints import FixAtoms
import copy 
from ase.constraints import FixAtoms
from asedias.plot import plot_dias
from asedias.utils import json_dump
from tqdm import tqdm
import sys

# check if running in IPython
def in_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False

if in_ipython():
    from IPython.display import clear_output


class Progress:
    RED = '\033[31m'
    GREEN = '\033[32m'
    RESET = '\033[0m'

    def __init__(self):
        self.ipython = "ipykernel" in sys.modules
        self.steps = []
        self.status = {}

    def set_steps(self, step_names):
        self.steps = step_names
        for step in step_names:
            self.status[step] = False
        self.display_summary()

    def display_summary(self):
        summary = []
        for step in self.steps:
            color = self.GREEN if self.status[step] else self.RED
            summary.append(f"{color}{step}{self.RESET}")
        print("Progress Summary:")
        print(" -> ".join(summary))

    def start_step(self, step_name, total=1):
        self.current_step = step_name
        self.total = total
        self.current = 0
        print(f"\nStarting step: {step_name}")
        if self.ipython:
            self.pbar = tqdm(total=total, leave=True)
        else:
            self.pbar = None

    def update(self, increment=1):
        self.current += increment
        if self.ipython and self.pbar:
            self.pbar.update(increment)
        else:
            print(f"{self.current_step}: {self.current}/{self.total}", end="\r")

    def finish_step(self):
        if self.ipython and hasattr(self, 'pbar') and self.pbar:
            self.pbar.close()
        self.status[self.current_step] = True
        self.display_summary()


class Fragment:
    """
    Represents a molecular fragment with index, charge, spin, name, constraints,
    and stores optimized energy and optimized atoms.
    """
    def __init__(self, index: list[int], charge: Optional[int] = None, spin: Optional[int] = None,
                 frag_name: Optional[str] = None, constraints: Optional[list[int]] = None):
        self.index = index
        self.charge = charge
        self.spin = spin
        self.name = frag_name
        self.constraints = constraints
        self.optimized_energy = None
        self.optimized_atoms = None


class aseDIAS:
    """
    Main class to run distortion-interaction analysis with fragment optimization caching.
    """
    def __init__(self, fragments: list[Fragment], images: Union[list[ase.Atoms], ase.io.Trajectory],
                 calc_attach: Callable[[ase.Atoms], None],
                 precalc_attach: Optional[Callable[[ase.Atoms], None]] = None,
                 job_name: Optional[str] = None, total_spin: Optional[int] = 0):
        self.fragments = fragments
        self.images = images
        self.calc_attach = calc_attach
        self.precalc_attach = precalc_attach
        self.job_name = job_name if job_name else str(uuid.uuid4())[:6]
        self.resultDict = dict()

        self.total_charge = sum([frag.charge for frag in self.fragments])
        self.total_spin = total_spin

    def run(self):
        from asedias import ParameterManager

        # -------------------------
        # Initialize progress
        # -------------------------
        progress = Progress()
        progress.set_steps(["Total Energies", "Optimize Fragments", "Fragment Energies"])

        # -------------------------
        # Step 1: Total energies
        # -------------------------
        progress.start_step("Total Energies", total=len(self.images))
        self.resultDict["total_energies"] = []
        calc_copy = copy.deepcopy(self.images[0])

        self.calc_attach(calc_copy, charge=self.total_charge, spin=self.total_spin)

        for atoms in self.images:
            total_energy = calc_copy.calc.get_potential_energy(atoms)
            self.resultDict["total_energies"].append(total_energy)
            progress.update()
        progress.finish_step()

        # -------------------------
        # Step 2: Optimize fragments
        # -------------------------
        progress.start_step("Optimize Fragments", total=len(self.fragments))
        for frag_order, frag in enumerate(self.fragments):
            frag_copy = self.images[0][frag.index]
            frag_name = frag.name if frag.name else f"frag_{frag_order}"
            frag_charge = getattr(frag, "charge", None)
            frag_spin = getattr(frag, "spin", None)

            self.resultDict[frag_name] = {}
            self.resultDict[frag_name]['frag_energies'] = []

            # Apply constraints
            if getattr(frag, "constraints", None):
                frag_copy.set_constraint(FixAtoms(indices=frag.constraints))

            # Pre-optimization if any
            if self.precalc_attach:
                self.precalc_attach(frag_copy, spin=frag_spin, charge=frag_charge)
                preopt = ParameterManager.optimizer(frag_copy)
                preopt.run(steps=ParameterManager.opt_steps,
                           fmax=ParameterManager.opt_fmax)

            # Optimization
            self.calc_attach(frag_copy, spin=frag_spin, charge=frag_charge)
            opt = ParameterManager.optimizer(frag_copy)
            assert opt.run(steps=ParameterManager.opt_steps, fmax=ParameterManager.opt_fmax), \
                f"optimization failed in {frag_name}"

            frag.optimized_atoms = frag_copy
            frag.optimized_energy = frag_copy.get_potential_energy()
            self.resultDict[frag_name]['optimized_energy'] = frag.optimized_energy
            progress.update()
        progress.finish_step()

        # -------------------------
        # Step 3: Fragment energies
        # -------------------------
        progress.start_step("Fragment Energies", total=len(self.fragments) * len(self.images))
        num_structs = len(self.resultDict['total_energies'])
        total_distortion_energies = [0.0] * num_structs
        interaction_energies = [0.0] * num_structs

        for frag_order, frag in enumerate(self.fragments):
            frag_name = frag.name if frag.name else f"frag_{frag_order}"
            frag_atoms_list = [atoms[frag.index] for atoms in self.images]
            frag_energy_series = []

            for atoms in frag_atoms_list:
                atoms.info.update({
                    "spin": getattr(frag, "spin", None),
                    "charge": getattr(frag, "charge", None)
                })
                frag_energy = frag.optimized_atoms.calc.get_potential_energy(atoms)
                frag_energy_series.append(frag_energy)
                progress.update()

            self.resultDict[frag_name]['frag_energies'] = frag_energy_series

            # ΔE_distortion_frag[i] = E_frag[i] - E_frag^opt
            distortion_per_frag = [e - frag.optimized_energy for e in frag_energy_series]
            self.resultDict[frag_name]['distortion_energies'] = distortion_per_frag

            # accumulate
            for i, d in enumerate(distortion_per_frag):
                total_distortion_energies[i] += d
                interaction_energies[i] += frag_energy_series[i]

        # ΔE_interaction[i] = E_total[i] - Σ_j E_frag[j][i]
        for i, e_total in enumerate(self.resultDict['total_energies']):
            interaction_energies[i] = e_total - interaction_energies[i]

        # Store total distortion & interaction energies
        self.resultDict['distortion_energies'] = total_distortion_energies
        self.resultDict['interaction_energies'] = interaction_energies
        progress.finish_step()


    def plot(self, include_fragments=None, yaxis_unit=None, geometric_indices=None, **kwargs):
        """
        Plot the results of the distortion interaction analysis.

        Parameters:
        - include_fragments: Whether to include fragments in the plot.
        - yaxis_unit: Unit for the y-axis.
        - geometric_indices: Indices of the geometric fragments to plot.
        - **kwargs: Additional keyword arguments for customization.
        """
        from asedias import ParameterManager
        marker = kwargs.get('marker', ParameterManager.marker)
        linestyle = kwargs.get('linestyle', ParameterManager.linestyle)
        hline = kwargs.get('hline', ParameterManager.hline)\
        
        include_fragments = include_fragments if include_fragments is not None else ParameterManager.include_fragments
        yaxis_unit = yaxis_unit if yaxis_unit is not None else ParameterManager.yaxis_unit

        return plot_dias(
            self.images, self.resultDict, include_fragments=include_fragments,
            yaxis_unit=yaxis_unit, geometric_indices=geometric_indices,
            marker=marker, linestyle=linestyle, hline=hline
            )
        
    def save(self, metadata: Optional[Union[dict, str]] = None):
        """
        Save the results of the distortion interaction analysis to a JSON file.

        Parameters:
        - filename: The name of the file to save the results.
        - metadata: Optional metadata to include in the saved file.
        """
        metadata = "asedias result" if not metadata else metadata

        json_dump(self.resultDict, self.job_name, metadata=metadata)


            