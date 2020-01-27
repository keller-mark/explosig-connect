import sys, os, numpy as np, pandas as pd

from ardnmf.ardnmf import ARDNMF, PRIORS, EXP_PRIOR

from scipy.stats import iqr
from tqdm.notebook import tqdm
from collections import defaultdict

# Adapted from https://github.com/cmsc828p-f18/Reproducing-Kim2016
def extract_signatures_Kim2016(sbs96_df):
    
    # Initialize ARDNMF hyperparameters
    n_init = 50 # default: 50 (following Kim, et al. paper)
    a = 10 # default: 10 (following SignatureAnalyzer v1.1 source code)
    beta = 1 # default: 1 (KL-divergence)
    tolerance = 1e-7 # default: 1e-7 (following Kim, et al. paper)
    max_iter = 100000 # default: 100000 (following SignatureAnalyzer v1.1 source code)
    prior = EXP_PRIOR # default: exponential (following Kim, et al. paper)
    tau = 0.001 # default: 0.001 (following SignatureAnalyzer v1.1 source code)
    n_permutations = 25000 # default: 25000
    random_seed = 1 # default: 1
    
    # Split hypermutators, which are defined to be samples where:
    # N_SNV > median(N_SNV) + 1.5*IQR (interquartile range)
    def split_hypermutators(X):
        # Repeat the hypermutator splitting until it's over
        sample_indices = list(range(X.shape[0]))
        hypermuts = True
        X_star = np.array(X)
        while hypermuts:
            # Row sums give N_SNV, then use median + 1.5 IQR to identify
            # hypermutated indices
            n_snvs = X_star.sum(axis=1)
            hypermut_thresh = np.median(n_snvs) + 1.5*iqr(n_snvs)
            hypermut_idx = set(np.where(n_snvs > hypermut_thresh)[0])
            hypermuts = len(hypermut_idx) > 0

            # Split the rows at hypermutated indices in two
            rows = []
            new_sample_idx = []
            for i, row in enumerate(X_star):
                if i in hypermut_idx:
                    rows.append(row/2.)
                    rows.append(row/2.)
                    new_sample_idx.append(sample_indices[i])
                    new_sample_idx.append(sample_indices[i])
                else:
                    rows.append(row)
                    new_sample_idx.append(sample_indices[i])

            # Create new X_star
            X_star = np.array(rows)
            sample_indices = new_sample_idx

        return X_star, sample_indices
    
    X = sbs96_df.values
    samples = sbs96_df.index
    N = len(samples)
    categories = list(sbs96_df.columns)
    cat_index  = dict(zip(categories, range(len(categories))))
    
    # Create mutation count matrix and split "hypermutators"
    X_star, sample_indices = split_hypermutators(X)
    
    np.random.seed(random_seed)
    
    models = []
    K_count = defaultdict(int)
    with tqdm(range(n_init), desc='Initializations') as t:
        for ini in t:
            # Create and fit model
            random_state = np.random.randint(2**32-1)
            model = ARDNMF(a=a, max_iter=max_iter, tol=tolerance,
                           norm_H=True, prior=prior, beta=beta,
                           random_state=random_state, tau=tau)
            H = model.fit_transform(X_star.T) # transpose because samples are on columns

            # Record and report
            models.append(model)
            K_count[model.k_eff_] += 1
            K_eff = max(K_count.keys(), key=lambda k: K_count[k])
      
    
    # Get the distribution of Ks, and choose the K that occurs the most often.
    # Then choose the model with LOWEST posterior probability for that K.
    # This goes against what they say in their paper but matches the SignatureAnalyzer
    # version 1 code.
    model = sorted([ m for m in models if m.k_eff_ == K_eff ], key=lambda m: m.obj_[-1])[0] 

    # Merge the exposures per sample together
    H_star = model.H_.T # transpose so samples are back on rows again
    H = np.zeros((N, K_eff))
    for i, s in enumerate(sample_indices):
        H[s] += H_star[i]

    # Reorder signatures by the maximum exposure to any category (seems to be
    # the convention followed by Kim, et al.)
    ordered_sigs = sorted(range(K_eff), key=lambda i: max(H[:, i].sum()*model.W_[:, i]), reverse=True)
    H = H[:, ordered_sigs]
    W = model.W_.T[ordered_sigs]
    
    sig_names = ["Extracted Signature %d" % i for i in range(W.shape[0])]
    exps_df = pd.DataFrame(index=sbs96_df.index.values.tolist(), columns=sig_names, data=H)
    sigs_df = pd.DataFrame(index=sig_names, columns=sbs96_df.columns.values.tolist(), data=W)
    
    exps_df.to_csv('Kim2016_reproduced_exposures.tsv', sep='\t')
    sigs_df.to_csv('Kim2016_reproduced_signatures.tsv', sep='\t')
    
    return sigs_df, exps_df