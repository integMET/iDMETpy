import pandas as pd
import numpy as np
import os

# Clear the workspace (not needed in Python)

# Step 1
# Metabolite name dictionary upload
file = "./metabodic.csv"
K = pd.read_csv(file).iloc[:, [0, 1]]

# List of csv files (metabolomic data)
# path = "./data/csv_set1"       # data set 1
path = "./data/csv_set2"         # data set 2

# Matchings
ALL = []
L = os.listdir(path)  # Assuming L is a list of filenames in the directory

for i in range(len(L)):
    # csv files
    file = os.path.join(path, L[i])
    D = pd.read_csv(file, header=0).iloc[1:, :]  # Skip the first row
    
    # Metabolite name
    M = D.iloc[:, 0].astype(str)
    
    # Matching with metabolite dictionary
    MT = np.full(len(M), np.nan)
    for k in range(len(M)):  # Metabolite name of differential data
        for j in range(K.shape[0]):  # Metabolite name of dictionary
            a = K.iloc[j, 1]
            b = a.split(';')  # Each splitted metabolite name of dictionary
            for l in range(len(b)):
                if M[k].lower() == b[l].lower():
                    MT[k] = K.iloc[j, 0]
    
    index = ~np.isnan(MT)
    all_data = pd.concat([pd.Series(MT, name='MT'), D], axis=1).loc[index]
    # All of differential table
    ALL.append(all_data)

ALL_dict = dict(zip(L, ALL))

# Data save
import pickle
with open('ALL.pickle', 'wb') as handle:
    pickle.dump(ALL_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
