import h5py

# Open file
with h5py.File('ZI-HTH_GABA_modisco_results.h5', 'r') as f:
    # List all groups
    print("Keys:", list(f.keys()))
    # Access a dataset
    data = f['neg_patterns']["pattern_0"]
    data.keys()  # List keys in the dataset
    data_save = data["contrib_scores"][:]
    data2 = f['pos_patterns']

df = pd.read_hdf('ZI-HTH_GABA_modisco_results.h5')