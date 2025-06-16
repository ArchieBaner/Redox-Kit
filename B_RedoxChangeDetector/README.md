Algorithm
---------

1. **Define the Redox Reaction**
    - Set the variable `reaction` to the redox reaction string to be analyzed (e.g., `"12Cr2O7^2- + 5SO3^2- = Cr^3+ + SO4^2-"`).

2. **Split Reaction into Reactants and Products**
    - Use the string `split(' = ')` method to separate the reaction into reactant and product parts.

3. **Parse Reactant Molecules**
    - For each reactant molecule string:
        - Extract any leading coefficient (e.g., `12` in `12Cr2O7^2-`).
        - Remove the charge annotation (e.g., `^2-`).
        - Use regex to extract element symbols and their counts (e.g., `Cr2O7` → `Cr:2`, `O:7`).
        - Multiply each element's count by the molecule's coefficient.
        - Extract the molecule's charge using `extract_charge`.
        - Assign oxidation states to each element using `assign_oxidation_states`, passing the element counts and total charge.
        - Store each molecule as a list of tuples: `(element, total_count, oxidation_state)`.

4. **Parse Product Molecules**
    - Repeat the same parsing steps as for reactants on the product molecules.

5. **Compare Oxidation States**
    - For each element present in both reactants and products:
        - Collect the assigned oxidation states from reactants and products.
        - Compare the oxidation states:
            - If the oxidation state increases, label as "Oxidation".
            - If the oxidation state decreases, label as "Reduction".
        - Record the element, its oxidation state change, and the type of change.

6. **Print Redox Changes**
    - For each element whose oxidation state changed:
        - Print the element symbol, its oxidation state in reactants and products, and whether it was oxidized or reduced.

7. **Program Entry Point**
    - The `main()` function defines the reaction string and calls the analysis function.
    - If the script is run as the main module, execute `main()`.

Variables Used
--------------
- `reaction`: The redox reaction string to be analyzed.
- `molecules`: List of parsed molecules, each as a list of (element, count, oxidation state) tuples.
- `composition`: List of (element, count) tuples for a molecule.
- `coeff`: Molecule coefficient (multiplies element counts).
- `charge`: The overall charge of a molecule.
- `ox_states`: Dictionary mapping elements to their assigned oxidation states.
- `changes`: List of tuples describing oxidation state changes for each element.
- `oxidation_state_db`: External dictionary mapping element symbols to lists of possible oxidation states.
- `species`: List of (element, count) tuples for a chemical species.
- `total_charge`: The overall charge of a species.
- `reactants`, `products`: Lists of parsed reactant and product molecules.
- `r_states`, `p_states`: Dictionaries mapping elements to their oxidation states in reactants and products.

Functions Used
--------------
- `main()`: Entry point; defines the reaction and calls the analysis function.
- `analyze_and_print_redox(reaction)`: Parses reactants and products, compares oxidation states, and prints changes.
- `parse_molecules(part)`: Parses a string of molecules, extracting element counts, coefficients, and oxidation states.
- `extract_charge(molecule)`: Extracts the ionic charge from a molecule string (e.g., "^2-").
- `assign_oxidation_states(species, total_charge)`: Assigns oxidation states to elements based on possible values and total charge.
- `compare_oxidation_states(reactants, products)`: Compares oxidation states between reactants and products, identifying redox changes.

Summary
-------
The `main()` function defines a redox reaction and triggers its analysis. The program parses each molecule, assigns oxidation states using a database, compares reactant and product states, and prints out which elements are oxidized or reduced. Helper functions handle parsing, charge extraction, and oxidation state assignment. All variables and functions used in the process are listed above for clarity.