import re
from oxidation_states import oxidation_state_db as db  # Dict: element symbol -> list of possible oxidation states
from itertools import product  # Cartesian product for generating combinations
from math import gcd  # Greatest common divisor
from collections import Counter  # Multiset for counting species

def parse_molecule_list_with_oxidation_states(reactant_str, product_str):
    """
    Parses lists of reactant and product molecules, annotating each with possible oxidation states for their constituent elements.
    Args:
        reactant_str (str): A string of reactant molecules separated by ' + ', e.g., "Fe^2+ + O2".
        product_str (str): A string of product molecules separated by ' + ', e.g., "Fe2O3".
    Returns:
        tuple: Two lists (reactants, products), each containing tuples of the form:
            ({element: oxidation_state, ...} or None, molecule_str)
            - The dictionary maps each element in the molecule to its assigned oxidation state, or None if no valid assignment is found.
            - molecule_str is the original molecule string.
    Notes:
        - The function attempts to assign oxidation states to each element in a molecule based on a database `db` of possible states.
        - It considers the molecule's overall charge (if specified, e.g., "^2-") when determining valid oxidation state combinations.
        - If multiple valid combinations exist, the one with the highest score (most likely) is chosen.
        - Requires `re` and `product` from `itertools`, and a `db` dictionary mapping elements to possible oxidation states.
    """
    
    def add_oxidation_state(molecule_list):
        """
        Determines and assigns possible oxidation states to elements in a list of chemical molecules.

        For each molecule in the input list, this function:
            - Extracts the overall charge from the molecule's string representation (e.g., "^2-").
            - Parses the chemical formula to identify elements and their counts.
            - Looks up possible oxidation states for each element from a database (`db`).
            - Tries all combinations of oxidation states for the elements such that the sum of oxidation numbers (weighted by atom counts) equals the molecule's charge.
            - Selects the combination with the highest score (favoring positive oxidation states).
            - Returns a list of tuples, each containing a dictionary mapping elements to their assigned oxidation states (or None if no valid combination is found) and the original molecule string.

        Args:
            molecule_list (list of str): List of molecule strings (e.g., ["Fe^3+", "H2O"]).

        Returns:
            list of tuple: Each tuple contains:
                - dict or None: Mapping of elements to their oxidation states, or None if no valid assignment.
                - str: The original molecule string.
        """
        result = []
        for molecule in molecule_list:
            # Extract charge (e.g., "^2-")
            charge_match = re.search(r'\^(\d+)([-+])', molecule)
            charge = int(charge_match.group(1)) * (-1 if charge_match and charge_match.group(2) == '-' else 1) if charge_match else 0
            # Extract elements and their counts
            matches = re.findall(r'([A-Z][a-z]?)(\d*)', molecule)
            elements = [e for e, _ in matches]
            counts = [int(n) if n else 1 for _, n in matches]
            states = [db.get(e, [0]) for e in elements]
            # Try all combinations of oxidation states for the elements
            valid_combos = [
                (sum((len(combo)-i)*ox for i, ox in enumerate(combo) if ox > 0),
                 {elements[i]: combo[i] for i in range(len(combo))})
                for combo in product(*states)
                if sum(c*n for c, n in zip(combo, counts)) == charge
            ]
            # Pick the combination with the highest score (most likely)
            result.append((max(valid_combos, key=lambda x: x[0])[1] if valid_combos else None, molecule))
        return result
    # Split reactants and products, then annotate with oxidation states
    reactants = reactant_str.split(' + ')
    products = product_str.split(' + ')
    return add_oxidation_state(reactants), add_oxidation_state(products)

def list_of_half_reactions(reactants_list, products_list):
    """
    Generates a list of half-reactions based on differences in element counts between reactants and products.
    Args:
        reactants_list (list): A list of tuples, each containing a dictionary of element counts and a string representation of a reactant.
        products_list (list): A list of tuples, each containing a dictionary of element counts and a string representation of a product.
    Returns:
        list: A list of dictionaries where each key is a tuple of (reactant_str, product_str) and each value is a tuple of (element, change_in_count).
    """
    
    half_list = []
    for r, p in product(reactants_list, products_list):
        r_dict, r_str = r
        p_dict, p_str = p
        if isinstance(r_dict, dict) and isinstance(p_dict, dict):
            for elem in r_dict:
                if elem in p_dict and r_dict[elem] != p_dict[elem]:
                    half_list.append({(r_str, p_str): (elem, p_dict[elem] - r_dict[elem])})
    return half_list

def balance_main_atoms(half_reaction):
    """
    Balances the main atoms (excluding O and H) in a list of half-reactions.
    Each half-reaction is expected to be a dictionary with a tuple of (lhs, rhs) as the key,
    and a tuple of (element, delta_ox) as the value. The function computes the lowest common
    multiple (LCM) of the atom counts for the specified element on both sides, and multiplies
    the species accordingly to balance the element.
    Args:
        half_reaction (list): A list of dictionaries, each representing a half-reaction in the form:
            {
                (lhs: str, rhs: str): (element: str, delta_ox: int)
            }
    Returns:
        list: A list of dictionaries with balanced species, in the same format as the input, but with
              coefficients adjusted to balance the specified element.
    """
    
    def count_element(molecule, elem):
        """
        Counts the number of atoms of a specific element in a chemical formula.

        Args:
            molecule (str): The chemical formula as a string (e.g., 'H2O', 'C6H12O6').
            elem (str): The element symbol to count (e.g., 'H', 'O', 'C').

        Returns:
            int: The total number of atoms of the specified element in the molecule.

        Example:
            count_element('C6H12O6', 'H')  # Returns 12
        """
        return sum(int(n) if n else 1 for e, n in re.findall(r'([A-Z][a-z]?)(\d*)', molecule) if e == elem)
    
    def get_leading_coeff(molecule):
        """
        Extracts the leading integer coefficient from a chemical molecule string.

        Args:
            molecule (str): The chemical molecule string, which may start with an integer coefficient.

        Returns:
            int: The leading coefficient if present; otherwise, 1.
        """
        m = re.match(r'^(\d+)', molecule)
        return int(m.group(1)) if m else 1
    
    def multiply_species(molecule, mult):
        """
        Multiplies the coefficient of a chemical species in its string representation.

        If the molecule string already starts with a coefficient (e.g., '2H2O'), the function multiplies this coefficient by `mult`.
        If there is no leading coefficient, it prepends the molecule with `mult`.

        Args:
            molecule (str): The string representation of the chemical species, possibly with a leading coefficient.
            mult (int): The multiplier to apply to the species.

        Returns:
            str: The updated molecule string with the coefficient multiplied or prepended as needed.
        """
        if mult == 1: return molecule
        m = re.match(r'^(\d+)(.*)', molecule)
        return f"{mult*int(m.group(1))}{m.group(2)}" if m else f"{mult}{molecule}"
    
    balanced = []
    for reaction in half_reaction:
        (lhs, rhs), (element, delta_ox) = list(reaction.items())[0]
        lhs_count = count_element(lhs, element) * get_leading_coeff(lhs)
        rhs_count = count_element(rhs, element) * get_leading_coeff(rhs)
        lcm_atoms = abs(lhs_count * rhs_count) // gcd(lhs_count, rhs_count) if lhs_count and rhs_count else 1
        lhs_mult, rhs_mult = lcm_atoms // lhs_count if lhs_count else 1, lcm_atoms // rhs_count if rhs_count else 1
        balanced.append({(multiply_species(lhs, lhs_mult), multiply_species(rhs, rhs_mult)): (element, delta_ox)})
    return balanced

def balance_oxygen(half_reaction, medium='acidic'):
    """
    Balances the oxygen atoms in a list of half-reactions by adding H2O molecules to the appropriate side.
    Args:
        half_reaction (list): A list of dictionaries, each representing a half-reaction in the form 
            { (lhs, rhs): (element, delta_ox) }, where lhs and rhs are strings representing the left-hand 
            and right-hand sides of the reaction, element is the element being balanced, and delta_ox is 
            the change in oxidation state.
        medium (str, optional): The medium in which the reaction occurs ('acidic' or 'basic'). 
            Currently only 'acidic' is supported. Defaults to 'acidic'.
    Returns:
        list: A list of dictionaries representing the balanced half-reactions, with H2O added to balance 
            the oxygen atoms.
    """
    
    def count_oxygen(molecule):
        """
        Counts the total number of oxygen atoms in a given chemical equation string.

        The function parses a chemical equation (which may contain multiple molecules separated by '+')
        and sums the number of oxygen atoms, taking into account any leading stoichiometric coefficients.

        Args:
            molecule (str): A string representing a chemical equation or molecule, e.g., "2H2O + O2".

        Returns:
            int: The total number of oxygen atoms present in the input string.

        Example:
            count_oxygen("2H2O + O2")  # Returns 4
        """
        total = 0
        for part in re.split(r'\s*\+\s*', molecule):
            m = re.match(r'^(\d+)?([A-Za-z0-9\^\-\+]+)$', part.strip())
            coeff = int(m.group(1)) if m and m.group(1) else 1
            formula = m.group(2) if m else part.strip()
            o_count = 0
            for e, n in re.findall(r'([A-Z][a-z]?)(\d*)', formula):
                if e == 'O':
                    o_count += int(n) if n else 1
            total += coeff * o_count
        return total

    balanced = []
    for reaction in half_reaction:
        (lhs, rhs), (element, delta_ox) = list(reaction.items())[0]
        lhs_o = count_oxygen(lhs)
        rhs_o = count_oxygen(rhs)
        diff = lhs_o - rhs_o
        if diff > 0:
            new_rhs, new_lhs = f"{rhs} + {diff}H2O", lhs
        elif diff < 0:
            new_lhs, new_rhs = f"{lhs} + {-diff}H2O", rhs
        else:
            new_lhs, new_rhs = lhs, rhs
        balanced.append({(new_lhs, new_rhs): (element, delta_ox)})
    return balanced

def balance_hydrogen(half_reaction, medium='acidic'):
    """
    Balances the hydrogen atoms in a list of half-reactions for redox equations, 
    by adding H⁺ (in acidic medium) or H₂O/OH⁻ (in basic medium) to the appropriate side.
    Args:
        half_reaction (list): 
            A list of dictionaries, each representing a half-reaction in the form 
            { (lhs, rhs): (element, delta_ox) }, where lhs and rhs are strings 
            representing the left-hand side and right-hand side of the reaction, 
            element is the element being balanced, and delta_ox is the change in oxidation state.
        medium (str, optional): 
            The medium in which the reaction occurs. Can be 'acidic' or 'basic'. 
            Defaults to 'acidic'.
    Returns:
        list: 
            A list of dictionaries representing the hydrogen-balanced half-reactions, 
            in the same format as the input.
    Notes:
        - In acidic medium, H⁺ ions are added to balance hydrogen atoms.
        - In basic medium, H₂O and OH⁻ are used to balance hydrogen atoms.
        - The function assumes that the input half-reactions are in the correct format.
    """
    
    def count_hydrogen(molecule):
        """
        Counts the total number of hydrogen (H) atoms in a given chemical molecule string.

        The molecule string can contain multiple parts separated by '+', each optionally
        preceded by a coefficient (e.g., '2H2O + H2SO4'). The function parses each part,
        extracts the coefficient and chemical formula, and sums the number of hydrogen atoms,
        taking into account the coefficients.

        Args:
            molecule (str): A string representing one or more chemical formulas, possibly with coefficients,
                            separated by '+' (e.g., '2H2O + H2SO4').

        Returns:
            int: The total number of hydrogen atoms in the molecule string.
        """
        total = 0
        for part in re.split(r'\s*\+\s*', molecule):
            m = re.match(r'^(\d+)?([A-Za-z0-9\^\-\+]+)$', part.strip())
            coeff = int(m.group(1)) if m and m.group(1) else 1
            formula = m.group(2) if m else part.strip()
            # Count H atoms in the formula
            h_count = 0
            for e, n in re.findall(r'([A-Z][a-z]?)(\d*)', formula):
                if e == 'H':
                    h_count += int(n) if n else 1
            total += coeff * h_count
        return total

    balanced = []
    for reaction in half_reaction:
        (lhs, rhs), (element, delta_ox) = list(reaction.items())[0]
        lhs_h = count_hydrogen(lhs)
        rhs_h = count_hydrogen(rhs)
        diff = lhs_h - rhs_h
        if medium == 'acidic':
            if diff > 0:
                new_rhs, new_lhs = f"{rhs} + {diff}H^+", lhs
            elif diff < 0:
                new_lhs, new_rhs = f"{lhs} + {-diff}H^+", rhs
            else:
                new_lhs, new_rhs = lhs, rhs
        elif medium == 'basic':
            if diff > 0:
                new_rhs, new_lhs = f"{rhs} + {diff}H2O", f"{lhs} + {diff}OH^-"
            elif diff < 0:
                new_lhs, new_rhs = f"{lhs} + {-diff}H2O", f"{rhs} + {-diff}OH^-"
            else:
                new_lhs, new_rhs = lhs, rhs
        else:
            new_lhs, new_rhs = lhs, rhs
        balanced.append({(new_lhs, new_rhs): (element, delta_ox)})
    return balanced

def balance_charge_with_electrons(half_reaction):
    """
    Balances the charge in a list of half-reactions by adding electrons to the appropriate side.

    Each half-reaction is represented as a dictionary with a tuple of (lhs, rhs) as the key and a tuple of (element, delta_ox) as the value.
    The function determines the number of electrons needed to balance the charge based on the change in oxidation state (delta_ox) and the number of atoms of the specified element.
    Electrons are added to the right-hand side (RHS) for oxidation (delta_ox > 0) and to the left-hand side (LHS) for reduction (delta_ox < 0).

    Args:
        half_reaction (list): A list of dictionaries, each representing a half-reaction in the form:
            { (lhs: str, rhs: str): (element: str, delta_ox: int) }

    Returns:
        list: A list of dictionaries, each representing the balanced half-reaction with electrons added, in the form:
            { (new_lhs: str, new_rhs: str): (element: str, delta_ox: int, reaction_type: str) }
        where reaction_type is "oxidation", "reduction", or "none".
    """
    def count_element(molecule, elem):
        """
        Counts the total number of atoms of a specified element in a chemical equation string.

        Args:
            molecule (str): A string representing a chemical equation or formula, possibly containing multiple compounds separated by '+' and optional coefficients.
            elem (str): The chemical symbol of the element to count (e.g., 'H', 'O', 'Na').

        Returns:
            int: The total count of the specified element in the given molecule string.

        Example:
            count_element("2H2O + O2", "O")  # Returns 4
        """
        total = 0
        for part in re.split(r'\s*\+\s*', molecule):
            m = re.match(r'^(\d+)?([A-Za-z0-9\^\-\+]+)$', part.strip())
            coeff = int(m.group(1)) if m and m.group(1) else 1
            formula = m.group(2) if m else part.strip()
            total += coeff * sum(int(n) if n else 1 for e, n in re.findall(r'([A-Z][a-z]?)(\d*)', formula) if e == elem)
        return total
    balanced = []
    for reaction in half_reaction:
        (lhs, rhs), (element, delta_ox) = list(reaction.items())[0]
        n_atoms = max(count_element(lhs, element), count_element(rhs, element))
        n_electrons = abs(delta_ox) * n_atoms
        if delta_ox > 0:  # Oxidation: electrons on RHS
            new_lhs, new_rhs, reaction_type = lhs, f"{rhs} + {n_electrons}e^-", "oxidation"
        elif delta_ox < 0:  # Reduction: electrons on LHS
            new_lhs, new_rhs, reaction_type = f"{lhs} + {n_electrons}e^-", rhs, "reduction"
        else:
            new_lhs, new_rhs, reaction_type = lhs, rhs, "none"
        balanced.append({(new_lhs, new_rhs): (element, delta_ox, reaction_type)})
    return balanced

def equalize_and_combine(oxid_half, red_half):
    """
    Combine two redox half-reactions into a balanced overall reaction.

    This function takes two half-reactions (oxidation and reduction), equalizes the number of electrons transferred by scaling the half-reactions appropriately, removes electrons from both sides, and cancels out identical species that appear on both sides of the equation. The result is a single balanced redox equation.

    Args:
        oxid_half (dict): The oxidation half-reaction, represented as a dictionary with a tuple of (lhs, rhs) as the key and any value.
        red_half (dict): The reduction half-reaction, represented similarly as a dictionary.

    Returns:
        str: The balanced overall redox reaction as a string in the form "lhs = rhs".

    Notes:
        - Electrons (e^-) are automatically removed from both sides after balancing.
        - Identical species present on both sides are canceled out.
        - The function assumes the input half-reactions are formatted as expected.
        - After combining, the function checks for a highest common factor (HCF) among all coefficients on both sides and scales the equation down if possible.
    """

    def electron_count(half):
        """
        Calculates the net electron count in a half-reaction.

        Args:
            half (dict): A dictionary representing a half-reaction, where the key is a tuple (lhs, rhs)
                with lhs and rhs being strings representing the left-hand side and right-hand side of the reaction.

        Returns:
            int: The net number of electrons, calculated as (electrons on lhs) - (electrons on rhs).

        Example:
            half = {('2e^- + Fe^3+', 'Fe^2+'): None}
            electron_count(half)  # Returns 2
        """
        (lhs, rhs), _ = list(half.items())[0]

        def count_e(side):
            """
            Counts the total number of electrons ('e^-') present on a given side of a redox equation.

            Args:
                side (str): A string representing one side of a redox equation, with electron terms
                            possibly prefixed by a coefficient (e.g., '2e^- + e^-').

            Returns:
                int: The total count of electrons on the given side.

            Example:
                count_e('2e^- + e^-')  # Returns 3
            """
            return sum(int(m.group(1)) if m.group(1) else 1
                       for part in side.split('+')
                       if (m := re.match(r'(\d+)?e\^\-$', part.strip())))
        return count_e(lhs) - count_e(rhs)
    n_e_oxid, n_e_red = abs(electron_count(oxid_half)), abs(electron_count(red_half))
    lcm_e = abs(n_e_oxid * n_e_red) // gcd(n_e_oxid, n_e_red) if n_e_oxid and n_e_red else 1

    def scale_half(half, mult):
        """
        Scales the coefficients of a half-reaction by a given multiplier.

        Args:
            half (dict): A dictionary representing a half-reaction, where the key is a tuple of two strings (lhs, rhs)
                representing the left-hand side and right-hand side of the reaction, respectively.
            mult (int): The multiplier to scale the coefficients by.

        Returns:
            tuple: A tuple containing two strings, representing the scaled left-hand side and right-hand side of the reaction.

        Example:
            half = {('2 H2O', 'O2 + 4 H+ + 4 e-'): None}
            scale_half(half, 3)
            # Returns: ('6 H2O', '3 O2 + 12 H+ + 12 e-')
        """
        (lhs, rhs), _ = list(half.items())[0]

        def scale_side(side):
            """
            Scales each component in a chemical equation side by a given multiplier.

            Args:
                side (str): A string representing one side of a chemical equation, with components separated by ' + '.

            Returns:
                str: The scaled side of the equation, where each component's coefficient is multiplied by 'mult'.
                      If a component does not have an explicit coefficient, 'mult' is used as its coefficient (unless mult == 1).

            Note:
                This function expects a variable 'mult' to be defined in the enclosing scope.
            """
            result = []
            for part in filter(None, map(str.strip, side.split(' + '))):
                m = re.match(r'^(\d+)(.*)', part)
                result.append(f"{int(m.group(1))*mult}{m.group(2)}" if m else (f"{mult}{part}" if mult > 1 else part))
            return ' + '.join(result)
        return scale_side(lhs), scale_side(rhs)
    oxid_lhs, oxid_rhs = scale_half(oxid_half, lcm_e // n_e_oxid if n_e_oxid else 1)
    red_lhs, red_rhs   = scale_half(red_half, lcm_e // n_e_red if n_e_red else 1)

    def merge_and_cancel(lhs, rhs):
        """
        Merges and cancels identical chemical species from the left-hand side (lhs) and right-hand side (rhs) of a chemical equation.

        This function:
        - Parses the lhs and rhs strings into species counts.
        - Removes any electron terms (e^-) from both sides.
        - Cancels out identical species present on both sides by subtracting the minimum count from each.
        - Returns the resulting lhs and rhs as formatted strings.

        Args:
            lhs (str): The left-hand side of the equation, with species separated by ' + '.
            rhs (str): The right-hand side of the equation, with species separated by ' + '.

        Returns:
            tuple[str, str]: The formatted lhs and rhs strings after merging and cancellation.
        """
        lhs_counts, rhs_counts = Counter(map(str.strip, lhs.split(' + '))), Counter(map(str.strip, rhs.split(' + ')))
        # Remove electrons
        for k in list(lhs_counts):
            if re.fullmatch(r'(\d+)?e\^\-', k): lhs_counts.pop(k)
        for k in list(rhs_counts):
            if re.fullmatch(r'(\d+)?e\^\-', k): rhs_counts.pop(k)
        # Cancel identical species
        for sp in set(lhs_counts) & set(rhs_counts):
            n = min(lhs_counts[sp], rhs_counts[sp])
            lhs_counts[sp] -= n
            rhs_counts[sp] -= n

        def fmt(counter):
            """
            Format a Counter-like dictionary into a string representation.

            Each key-value pair is formatted as 'vkey' if the value v > 1, or just 'key' if v == 1.
            Only items with a value greater than 0 are included.
            The resulting terms are joined with ' + '.

            Args:
                counter (dict): A dictionary (typically from collections.Counter) mapping keys to integer counts.

            Returns:
                str: A formatted string representing the nonzero items in the counter.
            """
            return ' + '.join(f"{v}{k}" if v > 1 else k for k, v in counter.items() if v > 0)
        
        return fmt(lhs_counts), fmt(rhs_counts)
    
    final_lhs, final_rhs = merge_and_cancel(oxid_lhs + ' + ' + red_lhs, oxid_rhs + ' + ' + red_rhs)

    # --- HCF reduction block ---
    def extract_coeffs(side):
        """
        Extracts all integer coefficients from a chemical equation side string.
        Returns a list of coefficients (default 1 if not present).
        """
        coeffs = []
        for part in filter(None, map(str.strip, side.split(' + '))):
            m = re.match(r'^(\d+)', part)
            coeffs.append(int(m.group(1)) if m else 1)
        return coeffs

    lhs_coeffs = extract_coeffs(final_lhs) if final_lhs else []
    rhs_coeffs = extract_coeffs(final_rhs) if final_rhs else []
    all_coeffs = lhs_coeffs + rhs_coeffs
    hcf = all_coeffs[0] if all_coeffs else 1
    for c in all_coeffs[1:]:
        hcf = gcd(hcf, c)
    if hcf > 1:
        def reduce_side(side):
            result = []
            for part in filter(None, map(str.strip, side.split(' + '))):
                m = re.match(r'^(\d+)(.*)', part)
                if m:
                    new_coeff = int(m.group(1)) // hcf
                    result.append(f"{new_coeff}{m.group(2)}" if new_coeff > 1 else m.group(2))
                else:
                    result.append(part)
            return ' + '.join(result)
        final_lhs = reduce_side(final_lhs)
        final_rhs = reduce_side(final_rhs)
    # --- end HCF reduction block ---

    # --- H2O, H+, OH- simplification block ---
    def simplify_species(lhs, rhs):
        """
        Simplifies the chemical species on both sides of a redox reaction by removing equal amounts of H2O, H^+, and OH^- that appear on both the left-hand side (lhs) and right-hand side (rhs).
        Args:
            lhs (str): The left-hand side of the reaction, as a string with species separated by ' + ' (e.g., "2H2O + H^+").
            rhs (str): The right-hand side of the reaction, as a string with species separated by ' + '.
        Returns:
            tuple: A tuple (lhs, rhs) of the simplified left-hand and right-hand side strings, with common H2O, H^+, and OH^- species removed from both sides.
        """
        
        def parse_counts(side):
            """
            Parses a chemical equation side and counts occurrences of specific species.

            Args:
                side (str): A string representing one side of a chemical equation, with species separated by ' + '.

            Returns:
                collections.Counter: A counter mapping each recognized species ('H2O', 'H^+', 'OH^-') to its total count.

            Notes:
                - Only counts 'H2O', 'H^+', and 'OH^-' species, optionally preceded by a coefficient.
                - Ignores unrecognized species or malformed parts.
                - Coefficient defaults to 1 if not specified.
            """
            counts = Counter()
            for part in filter(None, map(str.strip, side.split(' + '))):
                m = re.match(r'^(\d+)?(H2O|H\^\+|OH\^-)$', part)
                if m:
                    coeff = int(m.group(1)) if m.group(1) else 1
                    species = m.group(2)
                    counts[species] += coeff
            return counts

        lhs_counts = parse_counts(lhs)
        rhs_counts = parse_counts(rhs)
        for sp in ['H2O', 'H^+', 'OH^-']:
            n = min(lhs_counts[sp], rhs_counts[sp])
            if n > 0:
                # Remove n from both sides
                def remove_n(side, sp, n):
                    parts = []
                    for part in filter(None, map(str.strip, side.split(' + '))):
                        m = re.match(r'^(\d+)?(' + re.escape(sp) + r')$', part)
                        if m:
                            coeff = int(m.group(1)) if m.group(1) else 1
                            new_coeff = coeff - n
                            if new_coeff > 1:
                                parts.append(f"{new_coeff}{sp}")
                            elif new_coeff == 1:
                                parts.append(sp)
                            # If new_coeff <= 0, skip
                        else:
                            parts.append(part)
                    return ' + '.join(parts)
                lhs_ = remove_n(lhs, sp, n)
                rhs_ = remove_n(rhs, sp, n)
                lhs, rhs = lhs_, rhs_
        return lhs, rhs

    final_lhs, final_rhs = simplify_species(final_lhs, final_rhs)
    # --- end H2O, H+, OH- simplification block ---

    return f"{final_lhs} = {final_rhs}"

def main():
    """
    Main function to balance a redox reaction using the half-reaction method.

    Steps performed:
    1. Parses a sample redox reaction string into reactants and products.
    2. Annotates each species with oxidation states.
    3. Generates half-reactions for oxidation and reduction processes.
    4. Balances the main atoms involved in redox changes.
    5. Balances oxygen atoms by adding H2O molecules.
    6. Balances hydrogen atoms by adding H+ (for acidic solutions).
    7. Balances the charge by adding electrons to each half-reaction.
    8. Prints the balanced half-reactions with electron counts.
    9. Combines the half-reactions into the final balanced redox equation and prints it.

    Note: The function uses helper functions for parsing, balancing, and combining reactions.
    """
    # Example reaction to balance
    reaction_str = '4NO3^1- + 6I^1- = NO^0 + I2^0'
    
    print(f"\nUnbalanced reaction: \n{reaction_str}\n")
    # Step 1: Parse the string to get half reactions
    reactant_str, product_str = reaction_str.split(' = ')
    # Parse reactants and products, annotate with oxidation states
    list_of_reactants_w_oxs, list_of_products_w_oxs = parse_molecule_list_with_oxidation_states(reactant_str, product_str)
    # Step 1.5: Generate half reactions
    half_reactions = list_of_half_reactions(list_of_reactants_w_oxs, list_of_products_w_oxs)
    # Balance main atoms (the element being oxidized/reduced)
    balanced_main = balance_main_atoms(half_reactions)
    # Balance oxygen by adding H2O
    balanced_oxygen = balance_oxygen(balanced_main)
    # Balance hydrogen by adding H+ (acidic) or H2O/OH- (basic)
    balanced_hydrogen = balance_hydrogen(balanced_oxygen)
    # Balance charge by adding electrons
    balanced_electrons = balance_charge_with_electrons(balanced_hydrogen)
    print("Balanced Half Reactions with Electrons:")
    for half in balanced_electrons:
        (lhs, rhs), (element, delta_ox, reaction_type) = list(half.items())[0]
        print(f"{lhs} -> {rhs} | {element} ox state change: {delta_ox} | Type: {reaction_type}")
    print("\nFinal Balanced Reaction:")
    # Combine the two half-reactions into the final balanced equation
    final_reaction = equalize_and_combine(balanced_electrons[0], balanced_electrons[1])
    print(final_reaction)
    # Step 2: Balance non O_H atoms (not implemented here)

if __name__ == "__main__":
    main()