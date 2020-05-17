import numpy as np
import pubchempy as pcp
from func_timeout import func_set_timeout
from matchms.utils import is_valid_inchi, is_valid_inchikey, is_valid_smiles
import time

def lookup_metadata_completion(spectrum_in, search_depth=10):
    """Look for missing metadata and try to complete by PubChem lookups."""

    if spectrum_in is None:
        return None

    spectrum = spectrum_in.clone()

    # Read key annotation fields
    smiles = spectrum.get("smiles")
    inchi = spectrum.get("inchi")
    inchikey = spectrum.get("inchikey")
    compound_name = spectrum.get("compound_name")
    names_with_no_matches = []  # Avoid doing failed searches twice

    # Case 0 -- annotation metadata looks complete
    if smiles and inchikey and inchi:
        #print("Annotation metadata seems complete (no in-depth checks done).")
        return spectrum

    # Case 1 -- Inchikey match
    if is_valid_inchikey(inchikey):
        result, log_entry = lookup_by_inchikey(spectrum, search_depth)
        if result:
            add_entries(spectrum, result, log_entry)
            return spectrum
        print("No matches found for inchikey:", inchikey)

    # Case 2 -- Smiles match
    if smiles:
        result, log_entry = lookup_by_smiles(spectrum, search_depth)
        if result:
            add_entries(spectrum, result, log_entry)
            return spectrum
        print("Found no matches based on smiles.")

    # Case 3 -- Inchi match TODO: should be covered by derive functions!
    if is_valid_inchi(inchi):
        result, log_entry = lookup_by_inchi(spectrum, search_depth)
        if result:
            add_entries(spectrum, result, log_entry)
            return spectrum

        print("No matches found based on inchi:", inchi)

    # Case 4 -- Name match
    if compound_name and compound_name not in names_with_no_matches:
        result, log_entry = lookup_by_name(spectrum, search_depth)
        if result:
            add_entries(spectrum, result, log_entry)
            return spectrum
        if log_entry:
            if "No matches for compound name" in log_entry:
                # Avoid doing failed searches twice
                names_with_no_matches.append(compound_name)

    # Case 5 -- Formula match
    if inchi:
        result, log_entry = lookup_by_formula(spectrum, search_depth)
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
        spectrum.set("metadata_added_based_on", "PubChem match:" + log_entry)
        print("Added metadata based on: \n", log_entry)


@func_set_timeout(10)
def pubchem_get_properties(property_lst, search_item, search_item_name, search_depth):
    """Search PubChem."""
    results = pcp.get_properties(property_lst, search_item,
                                 search_item_name, search_depth=search_depth)
    return results


def lookup_by_inchikey(spectrum, search_depth):
    """Search for match by PubChem lookup based on inchikey."""
    inchikey = spectrum.get("inchikey")
    search_item_name = "inchikey"
    tstart = time.time()
    try:
        results = pubchem_get_properties(["InChIKey", "InChI", "IsomericSMILES",
                                          "MolecularFormula", "MolecularWeight"],
                                         inchikey[:14], search_item_name, search_depth=search_depth)
    except:
        if (time.time() - tstart) > 10:
            print("Timeout.")#" or no match found for search based on", search_item_name)
        results = []
    # Look for match with inchi AND inchikey
    result, log_entry = find_matches(results, spectrum,
                                     search_criteria=["inchi", "inchikey"],
                                     min_matches=2)
    if result:
        return result, log_entry

    # Look for match with inchi AND inchikey
    result, log_entry = find_matches(results, spectrum,
                                     search_criteria=["weight", "inchikey"],
                                     min_matches=2)
    if result:
        return result, log_entry

    # Look for match with full inchikey
    result, log_entry = find_matches(results, spectrum,
                                     search_criteria=["inchikey_full"],
                                     min_matches=1)
    if result:
        return result, log_entry
    return None, None


def lookup_by_smiles(spectrum, search_depth):
    """Search for match by PubChem lookup based on smiles."""
    smiles = spectrum.get("smiles")
    search_item_name = "smiles"
    try:
        results = pubchem_get_properties(["InChIKey", "InChI", "IsomericSMILES",
                                          "MolecularFormula", "MolecularWeight"],
                                         smiles, search_item_name, search_depth=search_depth)
    except:
        #print("Timeout or no match found for search based on", search_item_name)
        results = []

    # Look for match with inchikey
    result, log_entry = find_matches(results, spectrum, search_criteria=["inchikey"])
    if result:
        return result, log_entry + "Smiles match."

    # Accept any of the following: inchi, weight, formula
    result, log_entry = find_matches(results, spectrum,
                                     search_criteria=["inchi", "weight", "formula"],
                                     min_matches=1)
    if result:
        return result, log_entry + "Smiles match."

    # Accept unique match with Smiles
    if len({x.get("InChIKey", "")[:14] for x in results}) == 1:
        return results[0], "Unique Smiles match ({}).".format(smiles)
    return None, None


def lookup_by_inchi(spectrum, search_depth):
    """Search for match by PubChem lookup based on inchi."""
    inchi = spectrum.get("inchi")
    try:
        results = pcp.get_properties(["InChIKey", "InChI", "IsomericSMILES",
                                      "MolecularFormula", "MolecularWeight"],
                                     inchi, "inchi", search_depth=search_depth)
    except:
        # print("Timeout or no match found for search based on", search_item_name)
        results = []
    # Accept any of the following: inchi, weight, formula
    result, log_entry = find_matches(results, spectrum,
                                     search_criteria=["inchikey", "weight", "formula"],
                                     min_matches=1)
    if result:
        return result, log_entry

    # Accept unique match with InChI
    if len({x.get("InChIKey", "")[:14] for x in results}) == 1:
        return results[0], "Unique InChI match ({}).".format(inchi)
    return None, None


def lookup_by_name(spectrum, search_depth):
    """Search for match by PubChem lookup based on name."""
    compound_name = spectrum.get("compound_name")
    search_item_name = "name"
    if len(compound_name) <= 4:  # no meaningful name
        return None, None

    # Do PubChem lookup
    tstart = time.time()
    try:
        results = pubchem_get_properties(["InChIKey", "InChI", "IsomericSMILES",
                                          "MolecularFormula", "MolecularWeight"],
                                         compound_name, "name", search_depth=search_depth)
    except:
        if (time.time() - tstart) > 10:
            print("Timeout.")#
        #print("Timeout or no match found for search based on", search_item_name)
        results = []

    # Accept unique name match with two of the following
    if len({x.get("InChIKey", "")[:14] for x in results}) == 1:
        result, log_entry = find_matches(results, spectrum,
                                         search_criteria=["inchikey", "inchi", "weight", "formula"],
                                         min_matches=2)
        if result:
            return result, log_entry + "Matching compound name ({}).".format(compound_name)

    # Accept unique name match with any of the following
    if len({x.get("InChIKey", "")[:14] for x in results}) == 1:
        result, log_entry = find_matches(results, spectrum,
                                         search_criteria=["inchikey", "inchi", "weight", "formula"],
                                         min_matches=1)
        if result:
            return result, log_entry + "Matching compound name ({}).".format(compound_name)

    # Accept match with two of the following: inchi, inchikey, weight, formula
    result, log_entry = find_matches(results, spectrum,
                                     search_criteria=["inchikey", "inchi", "weight", "formula"],
                                     min_matches=2)
    if result:
        return result, log_entry + "Matching compound name ({}).".format(compound_name)

    # Accept any of the following: inchi, inchikey
    result, log_entry = find_matches(results, spectrum,
                                     search_criteria=["inchikey", "inchi"],
                                     min_matches=1)
    if result:
        return result, log_entry + "Matching compound name ({}).".format(compound_name)

    # # Accept one of <=5 name matches with one of the following (if unambiguous)
    # if len({x.get("InChIKey", "")[:14] for x in results}) <= 5:
        
    #     result, log_entry = find_matches(results, spectrum,
    #                                      search_criteria=["weight", "formula"],
    #                                      min_matches=1)
    #     if result:
    #         num_name_matches = len({x.get("InChIKey", "")[:14] for x in results})
    #         return result, log_entry + "Among {} entries matching compound name ({}).".format(num_name_matches,
    #                                                                                           compound_name)

    return None, "No matches for compound name."


def lookup_by_formula(spectrum, search_depth):
    """Search for match by PubChem lookup based on molecular formula."""
    formula = spectrum.get("formula")
    inchi = spectrum.get("inchi")
    search_item_name = "formula"
    if inchi and not formula:
        formula = inchi.strip().split("/")[1]

    if formula is None:
        return None, None

    try:
        results = pubchem_get_properties(["InChIKey", "InChI", "IsomericSMILES",
                                          "MolecularFormula", "MolecularWeight"],
                                         formula, search_item_name,
                                         listkey_count=search_depth)
    except:
        #print("Timeout or no match found for search based on", search_item_name)
        results = []

    # Accept any of the following: inchi, inchikey
    result, log_entry = find_matches(results, spectrum,
                                     search_criteria=["inchikey", "inchi"],
                                     min_matches=1)

    if result:
        return result, log_entry

    return None, None


def find_matches(pubchem_results, spectrum, search_criteria=["inchi", "inchikey"],
                 min_matches=2):
    """Find matching compound from pubchem search."""
    min_matches = min(min_matches, len(search_criteria))
    matching_functions = {"inchi": inchi_match,
                          "inchikey": inchikey_match,
                          "inchikey_full": inchikey_full_match,
                          "formula": formula_match,
                          "weight": weight_match}
    for result in pubchem_results:
        log_entry = ""
        matches = []
        for criterium in search_criteria:
            match_bool, log = matching_functions[criterium](result, spectrum)
            matches.append(match_bool)
            log_entry += log
        if sum(matches) >= min_matches:
            return result, log_entry
    return None, ""


def inchi_match(pubchem_result, spectrum, min_inchi_match=3):
    """True when inchi match is found"""
    inchi = spectrum.get("inchi")
    if inchi is None:
        return False, ""  #"No inchi."
    inchi_pubchem = pubchem_result.get("InChI")
    if likely_inchi_match(inchi, inchi_pubchem, min_agreement=min_inchi_match):
        log_entry = "Matching inchi (>= {} parts).".format(min_inchi_match)
        return True, log_entry
    return False, ""  #"No found inchi match."


def inchikey_match(pubchem_result, spectrum, min_inchikey_match=1):
    """True when inchikey match is found"""
    inchikey = spectrum.get("inchikey")
    inchikey_pubchem = pubchem_result.get("InChIKey")
    if inchikey is None or inchikey_pubchem is None:
        return False, ""

    if likely_inchikey_match(inchikey, inchikey_pubchem, min_agreement=min_inchikey_match):
        log_entry = "Matching inchikey (>= {} parts).".format(min_inchikey_match)
        return True, log_entry
    return False, ""  #"No found inchikey match."


def inchikey_full_match(pubchem_result, spectrum, min_inchikey_match=3):
    """True when full inchikey match is found"""
    inchikey = spectrum.get("inchikey")
    if inchikey is None:
        return False, ""  #"No inchikey."
    inchikey_pubchem = pubchem_result.get("InChIKey")
    if likely_inchikey_match(inchikey, inchikey_pubchem, min_agreement=min_inchikey_match):
        log_entry = "Matching inchikey (>= {} parts).".format(min_inchikey_match)
        return True, log_entry
    return False, ""  #"No found inchikey match."


def formula_match(pubchem_result, spectrum):
    """True when formula match is found"""
    inchi = spectrum.get("inchi")
    if inchi:
        formula = inchi.strip().split("/")[1]
    else:
        formula = None

    formula_pubchem = pubchem_result.get("MolecularFormula")

    if formula and formula_pubchem:
        if formula.upper().strip() == formula_pubchem.upper().strip():
            log_entry = "Matching formula."
            return True, log_entry
        return False, ""  #"No found formula match."
    return False, ""


def weight_match(pubchem_result, spectrum, mass_tolerance=2.0):
    """True when parent mass matches molecular weight."""
    weight = pubchem_result.get("MolecularWeight")
    parent_mass = spectrum.get("parent_mass")
    if parent_mass is None or parent_mass == 0:
        return False, ""  #"No parent mass found."
    if parent_mass and weight:
        weight_difference = weight - parent_mass
        if np.abs(weight_difference) <= mass_tolerance:
            log_entry = "Matching molecular weight ({:.1f} vs parent mass {:.1f}).".format(weight, parent_mass)
            return True, log_entry
    return False, ""  #"No found weight match."


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
    if inchi_1 is None or inchi_2 is None:
        return False

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

    if inchikey_1 is None or inchikey_2 is None:
        return False

    # Harmonize strings
    inchikey_1 = inchikey_1.upper().strip('"').replace(' ', '')
    inchikey_2 = inchikey_2.upper().strip('"').replace(' ', '')

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
