from __future__ import print_function
import re
from itertools import product
from oxidation_states import oxidation_state_db as db

def parse_molecule_list_with_oxidation_states(reactant_str, product_str):
    """
    Parses lists of reactant and product molecules, assigning possible oxidation states to each element in the molecules.
    Args:
        reactant_str (str): A string of reactant molecules separated by ' + ' (e.g., "H2O + O2").
        product_str (str): A string of product molecules separated by ' + ' (e.g., "H2O2").
    Returns:
        tuple: Two lists (reactants, products), each containing tuples of the form (oxidation_state_dict, molecule_str),
               where oxidation_state_dict maps element symbols to their assigned oxidation states for that molecule.
               If no valid oxidation state assignment is found, oxidation_state_dict is None.
    Notes:
        - The function expects a global variable `db` mapping element symbols to lists of possible oxidation states.
        - Molecule strings may include charges in the format '^2-' or '^1+'.
        - The function attempts to assign oxidation states such that the sum of oxidation states times atom counts equals the molecule's charge.
        - Only the combination with the highest sum of weighted oxidation states is returned if multiple valid assignments exist.
    """
    
    def add_oxidation_state(molecule_list):
        """
        Calculates and assigns possible oxidation states to each element in a list of molecules.

        For each molecule in the input list, this function:
        - Parses the molecular formula to extract elements and their counts.
        - Determines the overall charge of the molecule if specified (e.g., '^2-').
        - Looks up possible oxidation states for each element from a database `db`.
        - Finds valid combinations of oxidation states that satisfy the molecule's charge.
        - Selects the combination with the highest sum of weighted oxidation states.
        - Returns a list of tuples, each containing a dictionary mapping elements to their assigned oxidation states (or None if no valid combination is found) and the original molecule string.

        Args:
            molecule_list (list of str): List of molecular formula strings.

        Returns:
            list of tuple: Each tuple contains:
                - dict or None: Mapping of elements to their oxidation states, or None if no valid combination exists.
                - str: The original molecule string.
        """
        result = []
        for molecule in molecule_list:
            charge_match = re.search(r'\^(\d+)([-+])', molecule)
            charge = int(charge_match.group(1)) * (-1 if charge_match and charge_match.group(2) == '-' else 1) if charge_match else 0
            matches = re.findall(r'([A-Z][a-z]?)(\d*)', molecule)
            elements = [e for e, _ in matches]
            counts = [int(n) if n else 1 for _, n in matches]
            states = [db.get(e, [0]) for e in elements]
            valid_combos = [
                (sum((len(combo)-i)*ox for i, ox in enumerate(combo) if ox > 0),
                 {elements[i]: combo[i] for i in range(len(combo))})
                for combo in product(*states)
                if sum(c*n for c, n in zip(combo, counts)) == charge
            ]
            result.append((max(valid_combos, key=lambda x: x[0])[1] if valid_combos else None, molecule))
        return result
    reactants = reactant_str.split(' + ')
    products = product_str.split(' + ')
    return add_oxidation_state(reactants), add_oxidation_state(products)

def list_of_half_reactions(reactants_list, products_list):
    """
    Generates a list of half-reactions based on changes in element counts between reactants and products.
    Args:
        reactants_list (list): A list of tuples, each containing a dictionary of element counts and a string representation of a reactant.
        products_list (list): A list of tuples, each containing a dictionary of element counts and a string representation of a product.
    Returns:
        list: A list of dictionaries. Each dictionary maps a tuple of (reactant_str, product_str) to a tuple containing:
            - element (str): The element whose count changes.
            - delta (int): The change in the element's count (product - reactant).
            - reactant_count (int): The count of the element in the reactant.
            - product_count (int): The count of the element in the product.
    Note:
        This function assumes that the input lists contain tuples of (dict, str), where the dict maps elements to their counts.
    """
    
    half_list = []
    for r, p in product(reactants_list, products_list):
        r_dict, r_str = r
        p_dict, p_str = p
        if isinstance(r_dict, dict) and isinstance(p_dict, dict):
            for elem in r_dict:
                if elem in p_dict and r_dict[elem] != p_dict[elem]:
                    half_list.append({(r_str, p_str): (elem, p_dict[elem] - r_dict[elem], r_dict[elem], p_dict[elem])})
    return half_list

def print_half_reactions(half_list):
    """
    Prints formatted half-reactions (oxidation and reduction) from a list of reaction entries.
    Each entry in `half_list` is expected to be a dictionary where:
        - The key is a tuple: (reactant_str, product_str)
        - The value is a tuple: (element, delta, reactant_oxidation_state, product_oxidation_state)
    The function prints each unique half-reaction, indicating whether it is an oxidation or reduction,
    and includes the number of electrons transferred.
    Args:
        half_list (list): A list of dictionaries representing half-reactions.
    Example:
        half_list = [
            {('Fe²⁺', 'Fe³⁺'): ('Fe', 1, 2, 3)},
            {('Cu²⁺', 'Cu'): ('Cu', 2, 2, 0)}
        ]
        print_half_reactions(half_list)
    """
    
    seen = set()
    for entry in half_list:
        for (r_str, p_str), (elem, delta, r_ox, p_ox) in entry.items():
            key = (elem, r_str, p_str)
            if key in seen:
                continue
            seen.add(key)
            if r_ox < p_ox:
                # Oxidation
                print(f"{r_str} ({elem}: {r_ox}) → {p_str} ({elem}: {p_ox}) + {abs(p_ox - r_ox)}e⁻   [Oxidation]")
            else:
                # Reduction
                print(f"{r_str} ({elem}: {r_ox}) + {abs(r_ox - p_ox)}e⁻ → {p_str} ({elem}: {p_ox})   [Reduction]")

def main():
    reaction = "Cr2O7^2- + 3SO3^2- = 2Cr^3+ + 3SO4^2-"
    print("Reaction:\n", reaction)
    # Split the reaction into reactants and products
    reactants, products = reaction.split(' = ')
    reactants_parsed, products_parsed = parse_molecule_list_with_oxidation_states(reactants, products)
    half_list = list_of_half_reactions(reactants_parsed, products_parsed)
    print("\nHalf-Reactions:")
    print_half_reactions(half_list)

if __name__ == "__main__":
    main()
