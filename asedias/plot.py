import ase
import numpy as np
import matplotlib.pyplot as plt
from typing import Union
from ase import units


def relative_values(energy_series: Union[list, tuple], relative_index: Union[str, int] = "min") -> tuple:
    """
    Compute relative energy values based on a reference point.
    """
    _energy_series = np.array(energy_series)
    if relative_index == "min":
        _energy_series -= np.min(_energy_series)
    elif isinstance(relative_index, int):
        _energy_series -= _energy_series[relative_index]
    else:
        raise ValueError(f"Invalid relative_index value: {relative_index}")
    return tuple(_energy_series.tolist())


def geometric_analysis(images: list[ase.Atoms], geometric_indices: list[int]) -> tuple:
    """
    Compute geometric coordinates for x-axis (distance, angle, dihedral)
    """
    def _is_monotonic(lst: list):
        return all(lst[i] <= lst[i+1] for i in range(len(lst)-1)) or \
               all(lst[i] >= lst[i+1] for i in range(len(lst)-1))

    _index_diff = set(geometric_indices) - set(range(len(images[0])))
    assert _index_diff == set(), f"Invalid geometric_indices: {_index_diff}"

    if len(geometric_indices) == 2:
        xaxis = tuple(atoms.get_distance(*geometric_indices) for atoms in images)
        xaxis_type = "Distance"
    elif len(geometric_indices) == 3:
        xaxis = tuple(atoms.get_angle(*geometric_indices) for atoms in images)
        xaxis_type = "Angle"
    elif len(geometric_indices) == 4:
        xaxis = tuple(atoms.get_dihedral(*geometric_indices) for atoms in images)
        xaxis_type = "Dihedral Angle"
    else:
        raise ValueError(f"Expected 2/3/4 indices, got {len(geometric_indices)}")

    if not _is_monotonic(xaxis):
        import warnings
        warnings.warn("x-axis values are not monotonic. Plot may behave unexpectedly.", UserWarning)

    return xaxis_type, xaxis


def xaxis_formatter(images: list[ase.Atoms], geometric_indices: list = None) -> tuple[str, tuple]:
    """
    Returns x-axis label and values
    """
    _first_image = images[0]
    if not geometric_indices:
        xaxis_string = 'images'
        xaxis = tuple(range(len(images)))
    else:
        xaxis_type, xaxis = geometric_analysis(images, geometric_indices)
        xaxis_unit = {"Angle": "Degree", "Distance": "Å", "Dihedral Angle": "Degree"}
        xaxis_head = {"Angle": "∠", "Distance": "r", "Dihedral Angle": "∠"}
        symbols = _first_image.get_chemical_symbols()
        selected_symbols = [symbols[idx] for idx in geometric_indices]
        symbols_text = "-".join(selected_symbols)
        xaxis_string = f"{xaxis_head[xaxis_type]} {symbols_text} / {xaxis_unit[xaxis_type]}"
    return xaxis_string, xaxis


def husl_palette(n_colors: int) -> list[tuple]:
    """
    Generate HUSL-like palette without seaborn
    """
    import colorsys
    if n_colors < 2 or n_colors > 9:
        raise ValueError("n_colors must be in [2,9]")
    # evenly spaced hues in HUSL range, fixed saturation & lightness
    hues = np.linspace(0, 330, n_colors)
    palette = [colorsys.hls_to_rgb(h/360, 0.65, 0.65) for h in hues]
    return palette


def get_markers(n: int) -> list[str]:
    """
    Simple marker set for plotting
    """
    base = ['o','s','^','v','<','>','d','p','h']
    if n > len(base):
        base *= (n // len(base) + 1)
    return base[:n]


def plot_dias(images: list[ase.Atoms], resultDict: dict, include_fragments: bool = False,
              yaxis_unit: str = 'eV', relative_idx: Union[str,int] = 0,
              geometric_indices: list[int] = None, marker: bool = False,
              linestyle: bool = True, hline: bool = True) -> plt.Figure:
    """
    Optimized DIAS plotting (distortion, interaction, total energies already computed)
    Returns matplotlib Figure object.
    """

    # ASE unit conversion (eV -> desired unit)
    unit_map = {
        'EV': 1.0,
        'KCAL/MOL': 1 / (units.kcal / units.mol),
        'KJ/MOL': 1 / (units.kJ / units.mol),
        'HARTREE': 1 / units.Hartree
    }
    if yaxis_unit.upper() not in unit_map:
        raise ValueError(f"Unsupported yaxis_unit: {yaxis_unit}")
    unit_conv = unit_map[yaxis_unit.upper()]

    # x-axis
    xaxis_string, xaxis = xaxis_formatter(images, geometric_indices)
    yaxis_string = f"{'Rel. ' if relative_idx is not None else ''}Energy / {yaxis_unit}"

    fig, ax = plt.subplots()

    if hline:
        ax.plot(xaxis, [0]*len(xaxis), linewidth=0.5, color='black')

    # helper: unit conversion + relative
    def _prepare_energy(series):
        arr = np.array(series) * unit_conv
        if relative_idx is not None:
            arr = np.array(relative_values(arr, relative_idx))
        return arr

    # main energies
    main_energies = {
        'tot': resultDict['total_energies'],
        'dist': resultDict['distortion_energies'],
        'int': resultDict['interaction_energies']
    }
    colors = ['black', 'tab:blue', 'tab:red']
    markers_list = ['+', 'x', '4']
    linestyles_list = ['-', '-.', '--']

    for idx, (ename, series) in enumerate(main_energies.items()):
        ax.plot(
            xaxis, _prepare_energy(series),
            label=f"E_{ename[:3]}",
            color=colors[idx],
            marker=markers_list[idx] if marker else None,
            linestyle=linestyles_list[idx] if linestyle else None
        )

    # fragment energies
    if include_fragments:
        frag_keys = [k for k in resultDict.keys() if k not in ['total_energies','distortion_energies','interaction_energies']]
        frag_keys = sorted(frag_keys)
        n_frags = len(frag_keys)
        clr_list = husl_palette(n_frags)
        mkr_list = get_markers(n_frags)

        for idx, frag in enumerate(frag_keys):
            frag_series = _prepare_energy(resultDict[frag]['distortion_energies'])
            ax.plot(
                xaxis, frag_series,
                label=f"E_dis({frag})",
                color=clr_list[idx],
                marker=mkr_list[idx] if marker else None,
                linestyle=":" if linestyle else None
            )

    ax.set_xlabel(xaxis_string)
    ax.set_ylabel(yaxis_string)
    ax.legend()
    fig.tight_layout()

    return fig
