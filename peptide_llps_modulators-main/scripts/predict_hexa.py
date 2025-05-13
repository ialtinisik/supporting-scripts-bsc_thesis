from pathlib import Path
import pickle
import numpy as np
import pandas as pd
from sklearn.neural_network import MLPRegressor
from sklearn.svm import SVR
from sklearn.linear_model import Ridge
from sklearn.ensemble import RandomForestRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from train_model import featurize


def load_filtered_seqs():
    all_seqs = []
    for file in list(Path('./data/seqs').glob('*.txt')):
        with open(file, 'r') as f:
            seqs = f.readlines()
            seqs = [seq.strip() for seq in seqs]
            all_seqs.extend(seqs)
    return all_seqs


def predict_topk(descriptor: str, modelname: str, conc: str, optim: str):
    batch_size = 1000
    n_batch = 250
    k = 192
    assert conc in ['0.1', '1.0']
    assert optim in ['minimize', 'maximize']
    n_splits = 10

    all_seqs = load_filtered_seqs()

    # optionally filter subset to
    np.random.seed(9000)
    subset_indices = np.random.randint(0, len(all_seqs), size=n_batch*batch_size)
    print(subset_indices)
    all_seqs = [all_seqs[ix] for ix in subset_indices]

    models = []
    for fold in range(n_splits):
        with open(f'./out/models/{descriptor}-{modelname}_fold{fold}_conc{conc}.pkl', "rb") as f:
            model = pickle.load(f)
            models.append(model)

    topk_vals = np.ones((k,)) * np.inf
    topk_vals = topk_vals * -1 if optim == 'maximize' else topk_vals
    topk_seqs = [None] * k

    for index in range(0, len(all_seqs), batch_size):
        print('predicting seqs', index, 'conc', conc, '-', index+batch_size,
              '\t', 'best:', topk_seqs[np.argmax(topk_vals)], round(np.max(topk_vals), 4))
        # select next `batch_size` seqs, featurize, predict w/ model
        seqs = [all_seqs[ix] for ix in range(index, index+batch_size)]
        feats = featurize(seqs, descriptor=descriptor)

        # ensemble
        # preds = model.predict(feats)

        preds = np.mean([model.predict(feats) for model in models], axis=0)
        # print('preds', preds.shape)

        # concat topk with current batch (vals, seqs), partition (or via torch.topk)
        # top_k through np.partition is a more efficient way than to sort
        vals_topk_batch = np.concat([topk_vals, preds])
        seqs_topk_batch = topk_seqs + seqs

        if optim == 'maximize':
            topk_idx = np.argpartition(vals_topk_batch, -k)[-k:]
        elif optim == 'minimize':
            # bottom_k = -top_k (minimize)
            topk_idx = np.argpartition(vals_topk_batch, k)[:k]
        assert len(vals_topk_batch) == len(seqs_topk_batch)

        # select top_k values, keep track of vals & its sequences
        topk_vals = vals_topk_batch[topk_idx]
        topk_seqs = [seqs_topk_batch[ix] for ix in topk_idx]

        # check if this breaks?
        temp_df = pd.DataFrame({'sequences': topk_seqs, 'property': topk_vals})
        temp_df.to_csv(f'./out/intermittent_topk_seqs_preds_{conc}.csv')

    # sort final output since np.partition does keep topk sort itself
    if optim == 'maximize':
        sort_idx = np.argsort(-topk_vals)
    elif optim == 'minimize':
        sort_idx = np.argsort(topk_vals)

    topk_vals = topk_vals[sort_idx]
    topk_seqs = [topk_seqs[ix] for ix in sort_idx]
    [print(seq, round(val, 2))
     for seq, val in list(zip(topk_seqs[:32], topk_vals[:32]))]

    return topk_vals, topk_seqs


if __name__ == "__main__":
    # select best model from train / evaluation
    descriptor = 'peptides'
    modelname = 'svr'

    for optim in ['minimize', 'maximize']:
        for conc in ['0.1', '1.0']:
            vals, seqs = predict_topk(
                descriptor=descriptor, modelname=modelname, 
                conc=conc, optim=optim
            )

            # save output to csv
            df = pd.DataFrame({'sequences': seqs, 'property': vals})
            df.to_csv(f'./out/topk_seqs_{descriptor}-{modelname}_{conc}_{optim}.csv')
