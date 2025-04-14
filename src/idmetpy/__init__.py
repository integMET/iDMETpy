import numpy as np
import pandas as pd
import json
from scipy.stats import chi2_contingency
import networkx as nx

def load_json_data(fin_json_path, gclass_json_path, fname_json_path):
    """Load data from JSON files."""
    with open(fin_json_path, 'r', encoding='utf-8') as file:
        FIN = json.load(file)
    with open(gclass_json_path, 'r', encoding='utf-8') as file:
        gclass = json.load(file)
    with open(fname_json_path, 'r', encoding='utf-8') as file:
        FNAME = json.load(file)
    return FIN, gclass, FNAME

def process_metabolites(FIN, threshold=1.5):
    """Process metabolites based on thresholds."""
    D = []
    for i in FIN:
        df = pd.DataFrame.from_dict(FIN[i])
        A = df[df['ratio'] >= threshold]['name']
        B = df[df['ratio'] < (1 / threshold)]['name']
        D.append((A, B))
    return D

def calculate_odds_ratios(D, FNAME):
    """Calculate odds ratios and perform chi-square tests."""
    results = []
    
    for i in range(len(D)):
        for j in range(i + 1, len(D)):
            x1, x2 = D[i]
            y1, y2 = D[j]

            a = len(set(x1) & set(y1))
            b = len(set(x1) & set(y2))
            c = len(set(x2) & set(y1))
            d = len(set(x2) & set(y2))

            an = a
            bn = b
            cn = c
            dn = d

            A_metabolites = "|".join(list(set(x1) & set(y1)))
            B_metabolites = "|".join(list(set(x1) & set(y2)))
            C_metabolites = "|".join(list(set(x2) & set(y1)))
            D_metabolites = "|".join(list(set(x2) & set(y2)))

            if a * b * c * d == 0:
                a += 0.5
                b += 0.5
                c += 0.5
                d += 0.5

            contingency_table = np.array([[a, c], [b, d]])
            chi2, p_value, _, _ = chi2_contingency(contingency_table, correction=False)
            odds_ratio = np.log2((a * d) / (b * c))

            results.append({
                "i": i,
                "j": j,
                "odds_ratio": odds_ratio,
                "p_value": p_value,
                "a_count": an,
                "b_count": bn,
                "c_count": cn,
                "d_count": dn,
                "a": A_metabolites,
                "b": B_metabolites,
                "c": C_metabolites,
                "d": D_metabolites,
                "ip": FNAME[i],
                "jp": FNAME[j]
            })

    return pd.DataFrame(results)

def create_adjacency_matrix(results, D_length):
    """Create adjacency matrix for graph visualization."""
    OR = np.zeros((D_length, D_length))
    P = np.zeros((D_length, D_length))

    for _, row in results.iterrows():
        i, j = int(row["i"]), int(row["j"])
        OR[i, j] = row["odds_ratio"]
        P[i, j] = row["p_value"]

    return OR, P

def hello() -> str:
    return "Hello from idmetpy!"

def main():
    # Example usage
    FIN, gclass, FNAME = load_data("FIN.pkl", "gclass.pkl", "FNAME.pkl")
    D = process_metabolites(FIN)
    results = calculate_odds_ratios(D, FNAME)
    OR, P = create_adjacency_matrix(results, len(D))

    print("Odds Ratios:", OR)
    print("P-values:", P)

if __name__ == "__main__":
    main()
