"""
PubMed Phrase Comparison Script (Expanded with Custom Mode)
---------------------------------------------------------

This script helps researchers evaluate the effectiveness and uniqueness of PubMed search phrases.

It operates in two modes:
  A. PHRASE COMBINATION MODE:
     - Accepts a list of phrases.
     - Generates all combinations (AND/OR).
     - Compares overlap across them.
  
  B. CUSTOM SEARCH MODE:
     - Accepts a single, specific complex search string.
     - Retrieves PMIDs for that specific string only.

Output Files:
  1. *_counts.csv:     Result counts for the search.
  2. *_overlap.csv:    Pairwise overlap (only useful if running multiple phrase combos).
  3. *_pmid_lists.csv: The complete list of PMIDs for every search run.

Usage:
  - Configure the settings below (USER CONFIGURATION).
  - Run: python pubmed_term_overlap.py
"""

import requests
import time
import csv
from itertools import combinations

# ================= USER CONFIGURATION =================

# --- MODE SELECTION ---
# Set to True to ignore the PHRASES list and use the CUSTOM_SEARCH_STRING instead.
USE_CUSTOM_SEARCH = False

# If USE_CUSTOM_SEARCH is True, put your exact query here:
CUSTOM_SEARCH_STRING = '("long covid"[Title/Abstract] OR "post-acute sequelae"[Title]) AND "vaccine"[Title/Abstract]'

# If USE_CUSTOM_SEARCH is False, the script will combine these phrases:
PHRASES = [
  "home",
  "ventilation",
  "oxygen therapy"
]

# --- SEARCH SETTINGS ---
# Set to True to automatically append a filter for children/adolescents (infant-18yrs)
# NOTE: Set this to False if your Custom String already handles demographics.
APPLY_AGE_FILTER = True

MAX_COMBO_SIZE = 10              # Max number of phrases to combine (only for Phrase Mode)
OPERATOR = "AND"                # "AND" for strict match, "OR" for broad (only for Phrase Mode)

# --- API & OUTPUT SETTINGS ---
TOOL_NAME = "term_overlap_analyzer"
EMAIL = "your_email@example.com"       # Required by NCBI
USE_API_KEY = False                    # Set True if you have a key
API_KEY = "your_api_key_here"

SLEEP_BETWEEN_CALLS = 0.34             # Throttle to avoid NCBI 429 errors
OUTPUT_PREFIX = "pubmed_term_analysis"

# ======================================================


def generate_phrase_combinations(phrases, max_combination_size=2, operator="AND"):
    """
    Build Boolean combinations of phrases using AND or OR.
    Returns a list of tuples: (combined_search_string, [component_phrases])
    """
    operator = operator.upper()
    assert operator in {"AND", "OR"}, "Operator must be 'AND' or 'OR'."
    all_combos = []
    for r in range(1, max_combination_size + 1):
        for combo in combinations(phrases, r):
            combined = f' {operator} '.join(f'"{p}"' for p in combo)
            all_combos.append((combined, combo))
    return all_combos


def fetch_pmids(search_term, apply_age_filter=True):
    """
    Query PubMed via the esearch endpoint and return a set of PMIDs.
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    if apply_age_filter:
        age_filter = '("infant"[MeSH Terms] OR "child, preschool"[MeSH Terms] OR "child"[MeSH Terms] OR "adolescent"[MeSH Terms])'
        search_term = f"({search_term}) AND {age_filter}"
    
    params = {
        "db": "pubmed",
        "term": search_term,
        "retmax": 100000,  # Retreive up to 100k PMIDs
        "retmode": "json",
        "tool": TOOL_NAME,
        "email": EMAIL
    }
    if USE_API_KEY:
        params["api_key"] = API_KEY

    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()
        
        # Parse result
        if "esearchresult" in data and "idlist" in data["esearchresult"]:
            return set(data["esearchresult"]["idlist"])
        else:
            return set()
            
    except Exception as e:
        print(f"Error fetching for term: {search_term}\n{e}")
        return set()


def save_term_counts(results_dict, filename):
    with open(filename, "w", newline='', encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Search Term", "Component Phrases", "Result Count"])
        for term, data in results_dict.items():
            components_str = "; ".join(data["components"])
            writer.writerow([term, components_str, len(data["pmids"])])


def save_overlap_matrix(results_dict, filename):
    terms = list(results_dict.keys())
    # If only 1 term exists (Custom Mode), we can't do overlaps
    if len(terms) < 2:
        return 

    with open(filename, "w", newline='', encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Term 1", "Term 2", "Shared Count", "Union Count", "Jaccard %"])
        for t1, t2 in combinations(terms, 2):
            s1, s2 = results_dict[t1]["pmids"], results_dict[t2]["pmids"]
            shared = s1 & s2
            union = s1 | s2
            jaccard = 100 * len(shared) / len(union) if union else 0
            writer.writerow([t1, t2, len(shared), len(union), f"{jaccard:.1f}"])


def save_pmid_lists(results_dict, filename):
    with open(filename, "w", newline='', encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Search Term", "Count", "PMIDs (Sequential ->)"])
        for term, data in results_dict.items():
            pmid_list = sorted(list(data["pmids"]))
            row = [term, len(pmid_list)] + pmid_list
            writer.writerow(row)


# -------------- MAIN SCRIPT ------------------

if __name__ == "__main__":
    
    # 1. Determine Search Strategy
    if USE_CUSTOM_SEARCH:
        print(f"🔧 Custom Search Mode ENABLED.")
        print(f"   Query: {CUSTOM_SEARCH_STRING}")
        # Create a single entry list for the loop
        all_combos = [(CUSTOM_SEARCH_STRING, ["Custom Query"])]
    else:
        print("🔍 Generating phrase combinations...")
        all_combos = generate_phrase_combinations(PHRASES, max_combination_size=MAX_COMBO_SIZE, operator=OPERATOR)

    # 2. Run Searches
    results = {}
    print(f"🔎 Running {len(all_combos)} search(es) (Age Filter: {APPLY_AGE_FILTER})...")
    
    for search_term, components in all_combos:
        pmids = fetch_pmids(search_term, apply_age_filter=APPLY_AGE_FILTER)
        results[search_term] = {"components": components, "pmids": pmids}
        print(f"✓ Found {len(pmids)} PMIDs")
        time.sleep(SLEEP_BETWEEN_CALLS)

    # 3. Save Outputs
    print("💾 Saving term counts...")
    save_term_counts(results, f"{OUTPUT_PREFIX}_counts.csv")

    if len(results) > 1:
        print("💾 Saving overlap matrix...")
        save_overlap_matrix(results, f"{OUTPUT_PREFIX}_overlap.csv")
    else:
        print("ℹ️ Skipping overlap matrix (only 1 search term used).")

    print("💾 Saving full PMID records...")
    save_pmid_lists(results, f"{OUTPUT_PREFIX}_pmid_lists.csv")

    print("\n✅ Done!")