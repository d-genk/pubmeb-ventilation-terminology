"""
PubMed Phrase Comparison Script (Verbose UID-Slicing Edition)
-------------------------------------------------------------

This script helps researchers evaluate the effectiveness and uniqueness of PubMed search phrases.

**Crucial Update:** NCBI limits standard API pagination to 10,000 results for PubMed. 
To bypass this limitation robustly, this script "slices" the search across 
ranges of PMIDs (UIDs), fetching results in chronological chunks to retrieve 
datasets of any size (e.g., 28k+).

Output Files:
  1. *_counts.csv:     Result counts and Term Mapping.
  2. *_pmid_lists.csv: The complete list of PMIDs.
  3. *_overlap.csv:    (Optional) Overlap matrix.

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
USE_CUSTOM_SEARCH = True

# Your exact query:
CUSTOM_SEARCH_STRING = 'Chronic invasive mechanical ventilation OR Chronic invasive ventilation OR Chronic mechanical ventilation OR Chronic ventilation via tracheostomy OR Chronic ventilator dependence OR Home invasive ventilation OR Home mechanical ventilation OR Home ventilation OR Long-term mechanical ventilation'

# Phrase list (used only if USE_CUSTOM_SEARCH = False)
PHRASES = [
  "home",
  "ventilation",
  "oxygen therapy"
]

# --- SEARCH SETTINGS ---
# Set to True to automatically append a filter for children/adolescents.
# Set to False to match raw PubMed GUI searches exactly.
APPLY_AGE_FILTER = False

MAX_COMBO_SIZE = 10
OPERATOR = "AND"

# --- API & OUTPUT SETTINGS ---
TOOL_NAME = "term_overlap_analyzer"
EMAIL = "your_email@example.com"       # Required by NCBI
USE_API_KEY = False                    # Set True if you have a key
API_KEY = "your_api_key_here"

# We slice the UID space. 2 million PMIDs usually yields < 10,000 results per chunk.
# PubMed IDs currently go up to roughly 40 million.
UID_SLICE_SIZE = 2000000
MAX_PMID = 45000000                    
SLEEP_BETWEEN_CALLS = 0.35             # Throttle to avoid NCBI 429 errors
OUTPUT_PREFIX = "pubmed_term_analysis"

# ======================================================

def generate_phrase_combinations(phrases, max_combination_size=2, operator="AND"):
    operator = operator.upper()
    all_combos = []
    for r in range(1, max_combination_size + 1):
        for combo in combinations(phrases, r):
            combined = f' {operator} '.join(f'"{p}"' for p in combo)
            all_combos.append((combined, combo))
    return all_combos


def fetch_total_count_and_translation(base_url, base_params, search_term):
    """
    Does a single lightweight query to get the overall count and the GUI mapping.
    """
    params = base_params.copy()
    params["retmax"] = 0
    params["term"] = search_term
    
    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()
        
        if "esearchresult" not in data:
            return 0, ""
            
        total_count = int(data["esearchresult"].get("count", 0))
        translation = data["esearchresult"].get("querytranslation", "")
        return total_count, translation
        
    except Exception as e:
        print(f"   [!] Error getting baseline stats: {e}")
        return 0, ""


def fetch_pmids_by_slicing(search_term, apply_age_filter=True):
    """
    Query PubMed by slicing the query by PMID range to avoid the 10k limit.
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    if apply_age_filter:
        age_filter = '("infant"[MeSH Terms] OR "child, preschool"[MeSH Terms] OR "child"[MeSH Terms] OR "adolescent"[MeSH Terms])'
        search_term = f"({search_term}) AND {age_filter}"
    
    base_params = {
        "db": "pubmed",
        "retmode": "json",
        "tool": TOOL_NAME,
        "email": EMAIL
    }
    if USE_API_KEY:
        base_params["api_key"] = API_KEY

    print(f"   ...Contacting NCBI to initialize: {search_term[:50]}...")

    # 1. Get the official total count and mapping translation
    total_count, translation = fetch_total_count_and_translation(base_url, base_params, search_term)
    print(f"   ...Target: {total_count} total records.")
    
    if total_count == 0:
        return set(), translation, 0

    all_pmids = set()
    
    # 2. Iterate through UID space
    start_uid = 1
    slice_num = 1
    
    while start_uid <= MAX_PMID:
        end_uid = start_uid + UID_SLICE_SIZE - 1
        
        # Create a modified query constrained by a specific UID range
        sliced_query = f"({search_term}) AND {start_uid}:{end_uid}[UID]"
        
        slice_params = base_params.copy()
        slice_params["term"] = sliced_query
        slice_params["retmax"] = 9999  # Safe limit under 10k
        slice_params["retstart"] = 0   # Always 0, since we rely on UID filtering
        
        # Verbose progress indicator
        print(f"   [Slice {slice_num}] Scanning PMIDs {start_uid}-{end_uid}...", end=" ", flush=True)
        
        try:
            time.sleep(SLEEP_BETWEEN_CALLS)
            response = requests.get(base_url, params=slice_params)
            response.raise_for_status()
            data = response.json()
            
            if "esearchresult" in data:
                slice_count = int(data["esearchresult"].get("count", 0))
                
                # Check if we got IDs
                new_ids = []
                if "idlist" in data["esearchresult"] and data["esearchresult"]["idlist"]:
                    new_ids = data["esearchresult"]["idlist"]
                    all_pmids.update(new_ids)
                
                print(f"Found {len(new_ids)} hits. (Total so far: {len(all_pmids)})")

                # Safety check: If a single slice has >9999 results, our slice size is too big.
                if slice_count > 9999:
                    print(f"      [!] WARNING: Slice limit exceeded ({slice_count} > 9999). Some results in this range were truncated.")
                    
        except Exception as e:
            print(f"\n      [!] Error fetching slice {start_uid}-{end_uid}: {e}")
            
        start_uid += UID_SLICE_SIZE
        slice_num += 1
        
        # Optimization: Stop early if we have found all expected results
        if len(all_pmids) >= total_count:
            print(f"   ✓ Reached total count ({len(all_pmids)}/{total_count}). Stopping early.")
            break

    return all_pmids, translation, total_count


def save_term_counts(results_dict, filename):
    with open(filename, "w", newline='', encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Search Term", "Component Phrases", "Result Count", "Fetched Count", "Search Details"])
        for term, data in results_dict.items():
            components_str = "; ".join(data["components"])
            writer.writerow([
                term, 
                components_str, 
                data["total_count"], 
                len(data["pmids"]), 
                data.get("translation", "")
            ])


def save_overlap_matrix(results_dict, filename):
    terms = list(results_dict.keys())
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
        all_combos = [(CUSTOM_SEARCH_STRING, ["Custom Query"])]
    else:
        print("🔍 Generating phrase combinations...")
        all_combos = generate_phrase_combinations(PHRASES, max_combination_size=MAX_COMBO_SIZE, operator=OPERATOR)

    # 2. Run Searches
    results = {}
    print(f"🔎 Running {len(all_combos)} search(es) (Age Filter: {APPLY_AGE_FILTER})...")
    
    for search_term, components in all_combos:
        # Fetch using the UID slicing method
        pmids, translation, total_count = fetch_pmids_by_slicing(search_term, apply_age_filter=APPLY_AGE_FILTER)
        
        results[search_term] = {
            "components": components, 
            "pmids": pmids, 
            "translation": translation,
            "total_count": total_count
        }
        
        print(f"   ✓ Finished. Retrieved {len(pmids)} PMIDs (Total expected: {total_count})")
        time.sleep(SLEEP_BETWEEN_CALLS)

    # 3. Save Outputs
    print("💾 Saving term counts (with Search Details)...")
    save_term_counts(results, f"{OUTPUT_PREFIX}_counts.csv")

    if len(results) > 1:
        print("💾 Saving overlap matrix...")
        save_overlap_matrix(results, f"{OUTPUT_PREFIX}_overlap.csv")
    else:
        print("ℹ️ Skipping overlap matrix (only 1 search term used).")

    print("💾 Saving full PMID records...")
    save_pmid_lists(results, f"{OUTPUT_PREFIX}_pmid_lists.csv")

    print("\n✅ Done!")