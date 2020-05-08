import re
import numpy as np
import pubchempy as pcp
from matchms.importing import load_adducts
from matchms.utils import is_valid_inchi, is_valid_inchikey, is_valid_smiles


def lookup_metadata_completion(spectrum_in, search_depth=10):
    """Look for missing metadata and try to complete by PubChem lookups."""

    if spectrum_in is None:
        return None

    spectrum = spectrum_in.clone()

    # Read key annotation fields
    smiles = spectrum.get("smiles")
    inchi = spectrum.get("inchi")
    inchikey = spectrum.get("inchikey")
    compound_name = spectrum.get("name")
    parent_mass = spectrum.get("parent_mass")
    if inchi:
        formula = inchi.strip().split("/")[1]
    else:
        formula = None

    # Case - -- annotation metadata looks complete
    if smiles and inchikey and inchi and not inchi.startswith("issue"):
        #print("Annotation metadata seems complete (no in-depth checks done).")
        return spectrum

    # Case 1 -- Inchikey match
    if is_valid_inchikey(inchikey):
        try:
            results = pcp.get_properties(["InChIKey", "InChI", "IsomericSMILES",
                                          "MolecularFormula", "MolecularWeight"],
                                         inchikey[:14], "inchikey")
        except:
            results = None
        if results:
            if len(results) == 1:
                # Unique Inchikey match
                print("Found match based on inchikey.")
                add_entries(spectrum, results[0])
                return spectrum

            result = find_matches(results, inchi, inchikey)
            if result:
                add_entries(spectrum, result)
            else:
                print("No matches found for inchikey:", inchikey)

    # Case 2 -- Unique smiles match
    if smiles:
        try:
            results = pcp.get_properties(["InChIKey", "InChI", "IsomericSMILES",
                                          "MolecularFormula", "MolecularWeight"],
                                         smiles, "smiles")
        except:
            results = None
        if results:
            if len(results) == 1:
                print("Found match based on smiles.")
                add_entries(spectrum, results[0])
                return spectrum

            print("Found several matches based on smiles!")

    # Case 3 -- Inchi match
    if is_valid_inchi(inchi):
        try:
            results = pcp.get_properties(["InChIKey", "InChI", "IsomericSMILES",
                                          "MolecularFormula", "MolecularWeight"],
                                         inchi, "inchi")
        except:
            results = None
        if results:
            # Unique inchi match
            if len(results) == 1:
                print("Found match based on inchi.")
                add_entries(spectrum, results[0])
                return spectrum
            # Several inchi matches
            result = find_matches(results, inchi, inchikey)
            if result:
                add_entries(spectrum, result)
            else:
                print("No matches found for inchi:", inchi)

    # Case 4 --  Name match
    if compound_name:
        result, log_entry = lookup_by_name(compound_name, inchikey, inchi, formula,
                                           parent_mass, search_depth=search_depth)
        if result:
            add_entries(spectrum, result, log_entry)
            return spectrum

    # Case 4 --  Formula match
    if formula:
        result, log_entry = lookup_by_formula(formula, inchikey, inchi, search_depth=search_depth)
        if result:
            add_entries(spectrum, result, log_entry)
            return spectrum

    return spectrum


def add_entries(spectrum, lookup_result, log_entry=None):
    """Add found entries to metadata (unless already present)."""
    entry_change = False
    smiles = spectrum.get("smiles")
    inchi = spectrum.get("inchi")
    inchikey = spectrum.get("inchikey")
    smiles_pubchem = lookup_result.get("IsomericSMILES")
    inchi_pubchem = lookup_result.get("InChI")
    inchikey_pubchem = lookup_result.get("InChIKey")
    if not is_valid_smiles(smiles):
        spectrum.set("smiles", smiles_pubchem)
        entry_change = True
    if not is_valid_inchi(inchi):
        spectrum.set("inchi", inchi_pubchem)
        entry_change = True
    if not is_valid_inchikey(inchikey):
        spectrum.set("inchikey", inchikey_pubchem)
        entry_change = True
    # Add entry to log the reason for the metadata change
    if entry_change and log_entry:
        spectrum.set("lookup_completion", log_entry)


def lookup_by_name(name, inchikey=None, inchi=None, formula=None,
                   parent_mass=None, mass_tolerance=1.0, search_depth=10):
    """Return PubChem match if found."""

    compound_name = extract_compound_name(name)
    if len(compound_name) <= 4:
        return None, None

    # Do PubChem lookup
    try:
        results_name = pcp.get_properties(["InChIKey", "InChI", "IsomericSMILES",
                                           "MolecularFormula", "MolecularWeight"],
                                          compound_name, 'name',
                                          listkey_count=search_depth)
    except:
        print("No matches found.")
        return None, None

    if len(results_name) == 0:
        return None, None

    unique_match = len({x.get("InChIKey") for x in results_name}) == 1

    if inchikey and unique_match:
        # GOOD ENOUGH: one exact name match AND inchikey match
        if likely_inchikey_match(inchikey, results_name[0].get("InChIKey")):
            log_entry = "Match by name and inchikey"
            print(log_entry, "(Name:", compound_name, ")")
            return results_name[0], log_entry

    # GOOD ENOUGH: one exact name match AND inchi match
    if likely_inchi_match(inchi, results_name[0].get("InChI")) and unique_match:
        log_entry = "Match by name and inchi"
        print(log_entry, "(Name:", compound_name, ")")
        return results_name[0], log_entry

    # GOOD ENOUGH: one exact name match AND same molecular formula
    if formula and results_name[0].get("MolecularFormula") == formula and unique_match:
        log_entry = "Match by name and molecular formula"
        print(log_entry, "(Name:", compound_name, ")")
        return results_name[0], log_entry

    # GOOD ENOUGH: one exact name match AND same molecular weight
    if parent_mass:
        weight_difference = results_name[0].get("MolecularWeight") - parent_mass
        if np.abs(weight_difference) < mass_tolerance and unique_match:
            log_entry = "Match by name and molecular mass"
            print(log_entry, "(Name:", compound_name, ")")
            return results_name[0], log_entry

    # More than one match found
    if not unique_match:
        result = find_matches(results_name, inchi, inchikey)
        return result, None  #TODO add log_entry

    return None, None


def lookup_by_formula(formula, inchikey=None, inchi=None, mass_tolerance=1.0,
                      search_depth=10):
    """Return PubChem match if found."""

    # Do PubChem lookup
    try:
        results_formula = pcp.get_properties(["InChIKey", "InChI", "IsomericSMILES",
                                              "MolecularFormula", "MolecularWeight"],
                                             formula, "formula",
                                             listkey_count=search_depth)
    except:
        print("No matches found.")
        return None, None

    if len(results_formula) == 0:
        return None, None

    # Matches found
    if len(results_formula) > 0:
        result = find_matches(results_formula, inchi, inchikey)
        return result, None  # TODO add log entry


def extract_compound_name(name):
    """Clean compound name. Includes removing potential adduct."""
    known_adducts = load_adducts()
    adduct_original = name.split(' ')[-1]
    adduct = adduct_original.replace('\n', '')\
                            .replace(' ', '')\
                            .replace('[', '')\
                            .replace(']', '')\
                            .replace('*', '')
    if adduct in known_adducts["adducts_negative"] or adduct in known_adducts["adducts_positive"]:
        name = name.replace(adduct_original, '')

    # Clean string further:
    # remove type NCGC00180417-03_C31H40O16_
    name = re.split(r"[A-Z]{3,}[0-9]{8,}-[0-9]{2,}_[A-Z,0-9]{4,}_", name)[-1]
    # remove type NCGC00160232-01! or MLS001142816-01!
    name = re.split(r"[A-Z]{3,}[0-9]{8,}-[0-9]{2,}\!", name)[-1]
    # remove type Massbank:EA008813 or MassbankEU:EA008813
    name = re.split(r"[A-Z]{2,}[0-9]{5,}", name)[-1]
    # remove type HMDB:HMDB00943-1336
    name = re.split(r"HMDB:HMDB[0-9]{4,}-[0-9]{1,}", name)[-1]
    # remove type MoNA:662599
    name = re.split(r"HMDB:HMDB[0-9]{4,}-[0-9]{1,}", name)[-1]
    # ReSpect:PS013405 option1|option2|option3... [M+H]
    #TODO non ideal --> here just pick the last given compound name
    name = re.split(r"ReSpect:[A-Z]{2,}[0-9]{5,}.*\|", name)[-1]
    # remove type 0072_2-Mercaptobenzothiaz
    name = re.split(r"[0-9]{4}_", name)[-1]

    # Remove further non compound-name parts
    parts_remove = ["Spectral Match to",
                    "from NIST14",
                    "Massbank:"]
    for part in parts_remove:
        name = name.replace(part, "")

    return name.rstrip()


def find_matches(pubchem_results, inchi=None, inchikey=None,
                 min_inchi_match=3, min_inchikey_match=1):
    """Find matching compound from pubchem search."""
    assert is_valid_inchikey(inchikey) or inchi is not None, "Not enough input."
    if not is_valid_inchikey(inchikey):
        inchikey = "not provided"

    for result in pubchem_results:
        inchi_pubchem = result["InChI"]
        inchikey_pubchem = result["InChIKey"]

        match_inchikey = likely_inchikey_match(inchikey, inchikey_pubchem,
                                               min_agreement=min_inchikey_match)

        if inchi:
            match_inchi = likely_inchi_match(inchi, inchi_pubchem,
                                             min_agreement=min_inchi_match)
        else:
            match_inchi = False

        # GOOD ENOUGH: first time match between inchi and/or inchikey
        if match_inchi or match_inchikey:
            print("Match based on: ", match_inchi * "InChI", ", ",
                  match_inchikey * "InchiKey", ".")
            return result

    return None


def likely_inchi_match(inchi_1, inchi_2, min_agreement=3):
    """Try to match defective inchi to non-defective ones.

    Compares inchi parts seperately. Match is found if at least the first
    'min_agreement' parts are a good enough match.
    The main 'defects' this method accounts for are missing '-' in the inchi.
    In addition, differences between '-', '+', and '?'will be ignored.

    Args:
    --------
    inchi_1: str
        inchi of molecule.
    inchi_2: str
        inchi of molecule.
    min_agreement: int
        Minimum number of first parts that MUST be a match between both input
        inchi to finally consider it a match. Default is min_agreement=3.
    """
    if min_agreement < 2:
        print("Warning! 'min_agreement' < 2 has no discriminative power. Should be => 2.")
    if min_agreement == 2:
        print("Warning! 'min_agreement' == 2 has little discriminative power",
              "(only looking at structure formula. Better use > 2.")
    agreement = 0

    # Remove spaces and '"' to account for different notations.
    # Remove everything with little discriminative power.
    ignore_lst = ['"', ' ', '-', '+', '?']
    for ignore in ignore_lst:
        inchi_1 = inchi_1.replace(ignore, '')
        inchi_2 = inchi_2.replace(ignore, '')

    # Split inchi in parts.
    inchi_1_parts = inchi_1.split('/')
    inchi_2_parts = inchi_2.split('/')

    # Check if both inchi have sufficient parts (seperated by '/')
    if len(inchi_1_parts) >= min_agreement and len(
            inchi_2_parts) >= min_agreement:
        # Count how many parts agree well
        for i in range(min_agreement):
            agreement += (inchi_1_parts[i] == inchi_2_parts[i])

    return agreement == min_agreement


def likely_inchikey_match(inchikey_1, inchikey_2, min_agreement=2):
    """Try to match inchikeys.

    Compares inchikey parts seperately. Match is found if at least the first
    'min_agreement' parts are a good enough match.

    Args:
    --------
    inchikey_1: str
        inchikey of molecule.
    inchikey_2: str
        inchikey of molecule.
    min_agreement: int
        Minimum number of first parts that MUST be a match between both input
        inchikey to finally consider it a match. Default is min_agreement=2.
    """
    if min_agreement not in [1, 2, 3]:
        print("Warning! 'min_agreement' should be 1, 2, or 3.")
    agreement = 0

    # Harmonize strings
    inchikey_1 = inchikey_1.upper().replace('"', '').replace(' ', '')
    inchikey_2 = inchikey_2.upper().replace('"', '').replace(' ', '')

    # Split inchikey in parts.
    inchikey_1_parts = inchikey_1.split('-')
    inchikey_2_parts = inchikey_2.split('-')

    # Check if both inchikey have sufficient parts (seperated by '/')
    if len(inchikey_1_parts) >= min_agreement and len(
            inchikey_2_parts) >= min_agreement:
        # Count how many parts mostly agree
        for i in range(min_agreement):
            agreement += (inchikey_1_parts[i] == inchikey_2_parts[i])

    return agreement == min_agreement
