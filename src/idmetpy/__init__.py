import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
import networkx as nx

def load_data(fin_path, gclass_path, fname_path):
    """Load data from files."""
    FIN = pd.read_pickle(fin_path)
    gclass = pd.read_pickle(gclass_path)
    FNAME = pd.read_pickle(fname_path)
    return FIN, gclass, FNAME

def process_metabolites(FIN, threshold=1.5):
    """Process metabolites based on thresholds."""
    D = []
    for df in FIN:
        ratio = df.iloc[:, 1].astype(float)
        index1 = ratio >= threshold
        index2 = ratio < (1 / threshold)
        A = df.loc[index1, df.columns[0]].tolist()
        B = df.loc[index2, df.columns[0]].tolist()
        D.append((A, B))
    return D

def calculate_odds_ratios(D, FNAME):
    """Calculate odds ratios and perform chi-square tests."""
    results = []
    for i, (A1, B1) in enumerate(D):
        for j, (A2, B2) in enumerate(D):
            a = len(set(A1) & set(A2))
            b = len(set(A1) & set(B2))
            c = len(set(B1) & set(A2))
            d = len(set(B1) & set(B2))

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
                "a": a,
                "b": b,
                "c": c,
                "d": d,
                "i_name": FNAME[i],
                "j_name": FNAME[j]
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
