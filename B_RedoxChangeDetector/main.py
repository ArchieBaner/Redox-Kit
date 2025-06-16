import re
from itertools import product
from oxidation_states import oxidation_state_db

def parse_molecules(part):
    """
    Parses a string representing one or more chemical molecules, extracting their elemental composition,
    coefficients, and oxidation states.
    Args:
        part (str): A string containing one or more molecules separated by ' + '. Each molecule may have
            an optional leading coefficient and an optional charge (e.g., "2Fe^3+", "H2O").
    Returns:
        list: A list of molecules, where each molecule is represented as a list of tuples.
            Each tuple contains:
                - el (str): The element symbol.
                - count (int): The total count of the element in the molecule, multiplied by the molecule's coefficient.
                - ox_state (int or None): The assigned oxidation state for the element, or None if not assigned.
    Notes:
        - The function expects helper functions `extract_charge` and `assign_oxidation_states` to be defined elsewhere.
        - The molecule coefficient multiplies the total count of each element for reporting, but does not affect oxidation state assignment.
        - Charges are removed before parsing the elemental composition.
    """
    
    molecules = []
    for mol in part.split(' + '):
        m = re.match(r'(\d+)', mol.strip())
        coeff = int(m.group(1)) if m else 1
        # Remove charge part before parsing elements
        mol_no_charge = re.sub(r'\^.*', '', mol)
        elements = re.findall(r'([A-Z][a-z]*)(\d*)', mol_no_charge)
        # Only multiply the entire molecule by the coefficient, not each element's subscript
        composition = [(el, int(n) if n else 1) for el, n in elements]
        charge = extract_charge(mol)
        ox_states = assign_oxidation_states(composition, charge)
        # Now apply the molecule coefficient to the count for reporting, but not for oxidation state assignment
        molecules.append([(el, count * coeff, ox_states.get(el)) for el, count in composition])
    return molecules

def extract_charge(molecule):
    """
    Extracts the ionic charge from a molecule string formatted with a caret (^) followed by a number and a sign.

    Args:
        molecule (str): The string representation of the molecule, e.g., "Fe^3+" or "Cl^1-".

    Returns:
        int: The charge of the molecule as an integer. Positive for '+' and negative for '-'. Returns 0 if no charge is found.

    Examples:
        >>> extract_charge("Fe^3+")
        3
        >>> extract_charge("Cl^1-")
        -1
        >>> extract_charge("H2O")
        0
    """
    if (m := re.search(r'\^(\d+)([-+])', molecule)):
        return -int(m.group(1)) if m.group(2) == '-' else int(m.group(1))
    return 0

def assign_oxidation_states(species, total_charge):
    """
    Assigns oxidation states to elements in a chemical species based on possible oxidation states and total charge.

    Args:
        species (list of tuple): A list of tuples, where each tuple contains an element symbol (str) and its count (int) in the species.
        total_charge (int): The overall charge of the species.

    Returns:
        dict: A dictionary mapping each element symbol to its assigned oxidation state. If no valid assignment is found, returns a dictionary mapping each element to None.

    Notes:
        - Requires a global dictionary `oxidation_state_db` mapping element symbols to lists of possible oxidation states.
        - Uses a brute-force approach to try all combinations of possible oxidation states for the elements.
        - If only one element is present, assigns the total charge as its oxidation state.
    """
    if len(species) == 1:
        return {species[0][0]: total_charge}
    elements, counts = zip(*species)
    options = [oxidation_state_db[el] for el in elements]
    valid = [(sum(o * (len(combo) - i) for i, o in enumerate(combo) if o > 0), dict(zip(elements, combo)))
             for combo in product(*options) if sum(c * o for c, o in zip(counts, combo)) == total_charge]
    return max(valid, key=lambda x: x[0])[1] if valid else {el: None for el in elements}

def compare_oxidation_states(reactants, products):
    """
    Compares the oxidation states of elements in reactants and products to determine if they undergo oxidation or reduction.

    Args:
        reactants (list of list of tuples): Each reactant is represented as a list of tuples (element, count, oxidation_state).
        products (list of list of tuples): Each product is represented as a list of tuples (element, count, oxidation_state).

    Returns:
        list of tuples: Each tuple contains (element, reactant_oxidation_state, product_oxidation_state, 'Reduction' or 'Oxidation'),
        for elements whose oxidation state changes between reactants and products.
    """
    def collect_states(info): 
        return {el: ox for mol in info for el, _, ox in mol}
    
    r_states, p_states = collect_states(reactants), collect_states(products)
    print("Reactant Oxidation States:", r_states)
    print("Product Oxidation States:", p_states)
    return [(el, r_states[el], p_states[el], 'Reduction' if r_states[el] < p_states[el] else 'Oxidation')
            for el in r_states if el in p_states and r_states[el] != p_states[el]]

def analyze_and_print_redox(reaction):
    """
    Analyzes a redox reaction and prints the oxidation state changes for each element.

    Args:
        reaction (str): A string representing the redox reaction, with reactants and products
            separated by ' = '. Example: "Fe2+ + Cu2+ = Fe3+ + Cu+".

    The function parses the reactants and products, compares the oxidation states of each element,
    and prints the changes in oxidation states in a readable format.
    """
    reactants, products = map(parse_molecules, reaction.split(' = '))
    changes = compare_oxidation_states(reactants, products)
    print("Redox Changes:")
    for el, r_ox, p_ox, change in changes:
        print(f"{el}: {r_ox} -> {p_ox} ({change})")

def main():
    """
    Executes the main logic for analyzing a redox reaction.

    This function defines a sample redox reaction as a string and passes it to the
    analyze_and_print_redox function for analysis and output.
    """
    reaction = "12Cr2O7^2- + 5SO3^2- = Cr^3+ + SO4^2-"
    analyze_and_print_redox(reaction)

if __name__ == "__main__":
    main()
