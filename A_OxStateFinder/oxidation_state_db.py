"""
oxidation_state_db.py

This module defines a dictionary of known oxidation states for the elements
in the periodic table.

The oxidation states are compiled based on common literature and textbook sources.
This file can be imported wherever oxidation state data is required.
"""

oxidation_state_db = {
    'H': [+1, -1],
    'He': [0],
    'Li': [+1],
    'Be': [+2],
    'B': [+3],
    'C': [-4, -2, 0, +2, +4],
    'N': [-3, -2, -1, +1, +2, +3, +4, +5],
    'O': [-2, -1, 0, +1, +2],
    'F': [-1],
    'Ne': [0],
    'Na': [+1],
    'Mg': [+2],
    'Al': [+3],
    'Si': [-4, +2, +4],
    'P': [-3, +3, +5],
    'S': [-2, +2, +4, +6],
    'Cl': [-1, +1, +3, +5, +7],
    'Ar': [0],
    'K': [+1],
    'Ca': [+2],
    'Sc': [+3],
    'Ti': [+2, +3, +4],
    'V': [+2, +3, +4, +5],
    'Cr': [+2, +3, +6],
    'Mn': [+2, +4, +7],
    'Fe': [+2, +3],
    'Co': [+2, +3],
    'Ni': [+2, +3],
    'Cu': [+1, +2],
    'Zn': [+2],
    'Ga': [+1, +3],
    'Ge': [+2, +4],
    'As': [-3, +3, +5],
    'Se': [-2, +4, +6],
    'Br': [-1, +1, +3, +5],
    'Kr': [0, +2],
    'Rb': [+1],
    'Sr': [+2],
    'Y': [+3],
    'Zr': [+4],
    'Nb': [+3, +5],
    'Mo': [+4, +6],
    'Tc': [+7],
    'Ru': [+2, +3, +4, +6, +8],
    'Rh': [+1, +3],
    'Pd': [+2, +4],
    'Ag': [+1, +2],
    'Cd': [+2],
    'In': [+1, +3],
    'Sn': [+2, +4],
    'Sb': [-3, +3, +5],
    'Te': [-2, +4, +6],
    'I': [-1, +1, +3, +5, +7],
    'Xe': [0, +2, +4, +6, +8],
    'Cs': [+1],
    'Ba': [+2],
    'La': [+3],
    'Ce': [+3, +4],
    'Pr': [+3, +4],
    'Nd': [+3],
    'Pm': [+3],
    'Sm': [+2, +3],
    'Eu': [+2, +3],
    'Gd': [+3],
    'Tb': [+3, +4],
    'Dy': [+3],
    'Ho': [+3],
    'Er': [+3],
    'Tm': [+3],
    'Yb': [+2, +3],
    'Lu': [+3],
    'Hf': [+4],
    'Ta': [+5],
    'W': [+4, +6],
    'Re': [+4, +6, +7],
    'Os': [+3, +4, +6, +8],
    'Ir': [+3, +4],
    'Pt': [+2, +4],
    'Au': [+1, +3],
    'Hg': [+1, +2],
    'Tl': [+1, +3],
    'Pb': [+2, +4],
    'Bi': [+3, +5],
    'Po': [-2, +2, +4, +6],
    'At': [-1, +1, +3, +5, +7],
    'Rn': [0],
    'Fr': [+1],
    'Ra': [+2],
    'Ac': [+3],
    'Th': [+4],
    'Pa': [+5],
    'U': [+3, +4, +5, +6],
    'Np': [+3, +4, +5, +6],
    'Pu': [+3, +4, +5, +6],
    'Am': [+3, +4, +5, +6],
    'Cm': [+3, +4],
    'Bk': [+3, +4],
    'Cf': [+3],
    'Es': [+3],
    'Fm': [+3],
    'Md': [+2, +3],
    'No': [+2, +3],
    'Lr': [+3],
    'Rf': [+4],
    'Db': [+5],
    'Sg': [+6],
    'Bh': [+7],
    'Hs': [+8]
}

# Optional utility (if needed)
def get_possible_oxidation_states(element_symbol):
    """
    Returns the list of known oxidation states for a given element.

    Parameters
    ----------
    element_symbol : str
        The atomic symbol of the element (e.g., 'Fe', 'O', 'Cl').

    Returns
    -------
    list of int or None
        A list of possible oxidation states, or None if the element is unknown.
    """
    return oxidation_state_db.get(element_symbol)

