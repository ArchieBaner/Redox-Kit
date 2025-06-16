from __future__ import print_function
import re
from itertools import product
from oxidation_state_db import get_possible_oxidation_states

def parse(reactants):
    """
    Parses a string of chemical reactants, extracting coefficients, element symbols, and their counts,
    then computes and prints the oxidation states for each reactant.
    Args:
        reactants (str): A string representing chemical reactants separated by ' + ' (e.g., "2H2 + O2").
    Returns:
        dict or None: Oxidation states for each element, or None if not found.
    Side Effects:
        None
    Notes:
        - Assumes the existence of `oxidation_state` and `find_charge` functions.
        - Coefficients are optional; if absent, defaults to 1.
        - Element symbols are parsed using regular expressions.
    """
    # Extract leading coefficient if present
    coeff_match = re.match(r'(^\d+)', reactants.strip())
    coeff = int(coeff_match.group(1)) if coeff_match else 1

    # Identify elements and their counts
    elements = [
        [symbol, int(count) * coeff if count else coeff]
        for symbol, count in re.findall(r'([A-Z][a-z]?)(\d*)', reactants)
    ]

    # Compute oxidation states
    return oxidation_state(elements, find_charge(reactants))

def oxidation_state(species, charge):
    """
    Determines the most likely oxidation states for elements in a chemical species given its composition and overall charge.
    Args:
        species (list of tuple): A list of (element_symbol, count) pairs representing the chemical species.
        charge (int): The overall charge of the species.
    Returns:
        dict or None: A dictionary mapping element symbols to their assigned oxidation states if a valid combination is found,
        otherwise None.
    Notes:
        - The function uses `get_possible_oxidation_states(symbol)` to retrieve possible oxidation states for each element.
        - It considers all combinations of oxidation states and selects the one that matches the total charge and maximizes a scoring heuristic.
        - Requires `product` from itertools to generate combinations.
    """
    # Get possible oxidation states for each element; replace None with empty list
    state_options = [
        get_possible_oxidation_states(symbol) or []
        for symbol, _ in species
    ]
    element_counts = [count for _, count in species]

    valid_results = []
    # Try all combinations of oxidation states
    for state_combo in product(*state_options):
        total = sum(s * c for s, c in zip(state_combo, element_counts))
        if total == charge:
            # Heuristic: prefer higher oxidation states for earlier elements
            score = sum((len(species) - i) * ox for i, ox in enumerate(state_combo) if ox > 0)
            valid_results.append((score, {species[i][0]: state_combo[i] for i in range(len(species))}))

    if valid_results:
        # Return the combination with the highest score
        return max(valid_results, key=lambda x: x[0])[1]
    return None

def find_charge(reactant_species):
    """
    Extracts and returns the ionic charge from a chemical species string.
    The function searches for a charge notation in the format '^<number><sign>' 
    (e.g., '^2+', '^3-') within the input string. If found, it returns the 
    charge as an integer, with the appropriate sign. If no charge is found, 
    returns 0.
    Args:
        reactant_species (str): The chemical species string to parse.
    Returns:
        int: The ionic charge of the species. Positive for cations, negative for anions, 0 if no charge is specified.
    """
    num = re.search(r'\^(\d+)([-+])', reactant_species)
    if num:
        charge = int(num.group(1))
        return -charge if num.group(2) == '-' else charge
    return 0

def main():
    """
    Main function to parse a chemical formula and display the oxidation states of its elements.
    This function sets a sample chemical formula, parses it using the `parse` function,
    and prints the oxidation states of the elements if a valid combination is found.
    If no valid oxidation state combination is found, it notifies the user.
    """
   
    chemical = "SO4^2-"
    print(f"Parsing chemical: \n  {chemical}\n")
    result = parse(chemical)
    if result:
        print("Oxidation States:")
        for element, state in result.items():
            print(f"  {element}: {state}")
    else:
        print("No valid oxidation state combination found.")

if __name__ == "__main__":
    main()
