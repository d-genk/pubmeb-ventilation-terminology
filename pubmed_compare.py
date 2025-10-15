"""
PubMed Phrase Comparison Script
-------------------------------

This script helps researchers evaluate the effectiveness and uniqueness of PubMed search phrases.
It:
  1. Accepts a list of predefined medical search phrases.
  2. Constructs combinations of those phrases using Boolean logic ("AND" or "OR").
  3. Searches PubMed using the E-utilities API to retrieve PMIDs for each phrase.
  4. Compares overlap across search results for all combinations.
  5. Saves:
     - a CSV file listing how many articles were found for each phrase or combo
     - a CSV file showing pairwise overlaps (shared count and Jaccard index)

Usage Instructions:
- You must have Python 3 installed.
- Install the required Python library by running:
    pip install requests
- Customize the PHRASES list or OPERATOR (AND/OR) as needed.
- Add your email address (required by NCBI) and optionally an API key to increase rate limits.
- Run the script:
    python pubmed_term_overlap.py
- Output files will be saved in the same directory.

Next Steps (for researchers):
- Review `*_counts.csv` to identify phrases with the highest and lowest recall.
- Use `*_overlap.csv` to assess whether similar phrases retrieve the same literature.
- Identify redundant phrases or missed coverage areas.
- Refine your terminology strategy and test new candidate phrases iteratively.

"""

import requests
import time
import csv
from itertools import combinations

# -------------- USER CONFIGURATION ------------------

# Replace with your own list of phrases
PHRASES = [
    "pediatric tracheostomy",
    "mechanical ventilation",
    "chronic respiratory failure",
    "home ventilator support",
    "children with tracheostomy",
    "prolonged mechanical ventilation",
    "technology-dependent children",
    "pediatric home ventilation",
    "respiratory technology dependence",
    "long-term ventilation in children",
    "ventilator-dependent children",
    "chronic lung disease in pediatrics",
    "pediatric respiratory support",
    "invasive mechanical ventilation",
    "non-invasive ventilation",
    "tracheostomy-dependent children",
    "pediatric intensive care",
    "children with chronic illness",
    "pediatric pulmonary disease",
    "children with tracheostomy and ventilator",
    "long-term tracheostomy care"
]

TOOL_NAME = "term_overlap_analyzer"         # Name of your tool (for NCBI tracking)
EMAIL = "your_email@example.com"            # Required by NCBI — enter a valid email
USE_API_KEY = False                         # Set to True if you have a PubMed API key
API_KEY = "your_api_key_here"               # Optional: increases request rate limit

MAX_COMBO_SIZE = 2                          # Max number of phrases to combine (1 to N)
OPERATOR = "AND"                            # Use "AND" for strict match, "OR" for broad
SLEEP_BETWEEN_CALLS = 0.34                  # Seconds to wait between API calls
OUTPUT_PREFIX = "pubmed_term_analysis"      # Base name for output files

# ----------------------------------------------------


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


def fetch_pmids(search_term):
    """
    Query PubMed via the esearch endpoint and return a set of PMIDs.
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "pubmed",
        "term": search_term,
        "retmax": 1000,  # Max results per query; raise with caution
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
        return set(data["esearchresult"]["idlist"])
    except Exception as e:
        print(f"Error fetching for term: {search_term}\n{e}")
        return set()


def save_term_counts(results_dict, filename):
    """
    Save a CSV with search terms and their result counts.
    """
    with open(filename, "w", newline='', encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Search Term", "Component Phrases", "Result Count"])
        for term, data in results_dict.items():
            writer.writerow([term, "; ".join(data["components"]), len(data["pmids"])])


def save_overlap_matrix(results_dict, filename):
    """
    Save a CSV matrix of pairwise overlaps between all term results.
    """
    terms = list(results_dict.keys())
    with open(filename, "w", newline='', encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Term 1", "Term 2", "Shared Count", "Union Count", "Jaccard %"])
        for t1, t2 in combinations(terms, 2):
            s1, s2 = results_dict[t1]["pmids"], results_dict[t2]["pmids"]
            shared = s1 & s2
            union = s1 | s2
            jaccard = 100 * len(shared) / len(union) if union else 0
            writer.writerow([t1, t2, len(shared), len(union), f"{jaccard:.1f}"])


# -------------- MAIN SCRIPT ------------------

if __name__ == "__main__":
    print("🔍 Generating search phrase combinations...")
    all_combos = generate_phrase_combinations(PHRASES, max_combination_size=MAX_COMBO_SIZE, operator=OPERATOR)
    
    results = {}
    print(f"🔎 Running {len(all_combos)} PubMed searches...")
    for search_term, components in all_combos:
        pmids = fetch_pmids(search_term)
        results[search_term] = {"components": components, "pmids": pmids}
        print(f"✓ {search_term} → {len(pmids)} results")
        time.sleep(SLEEP_BETWEEN_CALLS)

    print("💾 Saving term counts...")
    save_term_counts(results, f"{OUTPUT_PREFIX}_counts.csv")

    print("💾 Saving overlap matrix...")
    save_overlap_matrix(results, f"{OUTPUT_PREFIX}_overlap.csv")

    print("\n✅ Done! You can now explore:")
    print(f"  • {OUTPUT_PREFIX}_counts.csv — Number of results per phrase")
    print(f"  • {OUTPUT_PREFIX}_overlap.csv — Overlap between search results (pairwise)")
