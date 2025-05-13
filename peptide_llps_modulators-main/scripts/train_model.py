import pickle

import numpy as np
import pandas as pd
import peptides
from rdkit import Chem
from rdkit.Chem import AllChem, rdFingerprintGenerator
from rdkit.Chem.Descriptors import CalcMolDescriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.linear_model import Lasso, LinearRegression, Ridge
from sklearn.metrics import (mean_absolute_error, r2_score,
                             root_mean_squared_error)
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR


def calc_error_metrics(preds, labels):
    mae = mean_absolute_error(preds, labels)
    rmse = root_mean_squared_error(preds, labels)
    r2 = r2_score(preds, labels)
    return mae, rmse, r2


def featurize(sequences: list[str], descriptor: str):
    # most of these are composition-only descriptors,
    # and do not explicitly take sequence into account
    if descriptor == "peptides":
        # m x n == (n_peptides, 102)
        features = np.array([
            list(peptides.Peptide(s).descriptors().values()) for s in sequences
        ])
        return features

    elif descriptor == "rdkit":
        # m x n == (n_peptides, 200)
        features = np.array([
            list(CalcMolDescriptors(Chem.MolFromFASTA(s)).values()) for s in sequences

        ])
        return features

    elif descriptor == "ecfp":
        # missing convert sequence to smiles
        featurizer = rdFingerprintGenerator.GetMorganGenerator(
            fpSize=2048, radius=2
        )
        features = np.array([
            featurizer.GetFingerprint(Chem.MolFromFASTA(s)) for s in sequences
        ])
        return features

    else:
        # TODO implement sequence-based descriptor
        raise NotImplementedError


def train_model(data: pd.DataFrame, propname: str, descriptor: str, modelname: str, conc: str):
    assert 'peptide' in data.columns and propname in data.columns
    sequences = data['peptide']
    properties = data[propname]
    n_splits = 10

    features = featurize(sequences, descriptor=descriptor)
    train_maes, train_rmses = [], []
    test_maes, test_rmses = [], []

    for fold in range(n_splits):
        # train test split ## (maybe validation)
        X_train, X_test, y_train, y_test = train_test_split(
            features, properties, test_size=0.1, random_state=9000+fold
        )

        # scaling (may not be necessary)
        # scaler = StandardScaler()
        # scaler = scaler.fit(X_train)
        # X_train = scaler.transform(X_train)
        # X_test = scaler.transform(X_test)

        print("train: X", X_train.shape, '\t', 'y', y_train.shape)
        print("test: X", X_test.shape, '\t', 'y', y_test.shape)
        # print("test sequences\n", y_test)

        # models #
        if modelname == "rf":
            model = RandomForestRegressor(
                # n_estimators=100, min_samples_split=2,
                # min_samples_leaf=1, min_samples_split=2,
                # random_state=42
            )
        elif modelname == "linear":
            model = LinearRegression()
        elif modelname == "ridge":
            model = Ridge()
        elif modelname == "lasso":
            model = Lasso()
        elif modelname == "svr":
            model = SVR()
        else:
            raise NotImplementedError

        ## hyperparameter tuning ##
            # hyperparameters = {}
            # GridSearchCV

        ## model selection ##
            # pick best model from GridSearchCV
            # train again on full data
            # model = ...

        ## model training ##
        model = model.fit(X_train, y_train)
        print('trained', modelname, 'w/ descriptor:', descriptor)

        # optional: evaluate metrics on train
        preds_train = model.predict(X_train)
        mae, rmse, r2 = calc_error_metrics(preds_train, y_train)
        # print('TRAIN', 'mae', mae:.2f, 'rmse', rmse:.2f, 'r2', r2:.2f)
        print(f'TRAIN mae {mae:.3f} rmse {rmse:.3f} r2 {r2:.3f}')
        train_maes.append(mae)
        train_rmses.append(rmse)

        # predict for test set
        preds = model.predict(X_test)

        # evaluate test set (mae/rmse/r2)
        mae, rmse, r2 = calc_error_metrics(preds, y_test)
        print(f'TEST  mae {mae:.3f} rmse {rmse:.3f} r2 {r2:.3f} \n')

        test_rmses.append(rmse)
        test_maes.append(mae)
        # save model #
        with open(f'./out/models/{descriptor}-{modelname}_fold{fold}_conc{conc}.pkl', "wb") as f:
            pickle.dump(model, f, protocol=5)

    # predict on whole dataset in separate script #
    # with open(f'./out/models/{descriptor}-{modelname}.pkl', "rb") as f:
    #     model = pickle.load(f)

    avg_train_mae = np.round(np.mean(train_maes), 3)
    avg_train_rmse = np.round(np.mean(train_rmses), 3)
    avg_test_mae = np.round(np.mean(test_maes), 3)
    avg_test_rmse = np.round(np.mean(test_rmses), 3)
    # maybe to table & save results
    return {'modelname': modelname, 'descriptor': descriptor,
            'avg_train_mae': avg_train_mae, 'avg_test_mae': avg_test_mae,
            'avg_train_rmse': avg_train_rmse, 'avg_test_rmse': avg_test_rmse,
            'conc': conc}


if __name__ == "__main__":
    # filepaths = [
    #     './out/csc_integrals_sequences_SYN18.csv',
    #     './out/csc_integrals_sequences_SYN24.csv',
    #     './out/csc_integrals_sequences_SYN42.csv',
    # ]
    # sequence_prop_df = pd.concat([
    #     [pd.read_csv(path) for path in filepaths]
    # ], dim=0)
    # column name of dataframe - expects 'sequence'
    # propname = 'CSC' # or 'csc'

    # df generated from src/toy_targets.py w/ filtered enumerator.py seqs
    # sequence_prop_df = pd.read_csv('./out/seqs_toy_targets.csv')

    for conc in ['0.1', '1.0']:
        # sequence_prop_df = pd.read_csv('./data/PSM0010_output_1.0.csv')
        sequence_prop_df = pd.read_csv(f'./data/PSM0010_output_{conc}.csv')

        # average across replicates
        sequence_prop_df = sequence_prop_df.groupby('peptide').agg(
            normintegral_avg=('normintegral', 'mean'))
        sequence_prop_df = sequence_prop_df.reset_index()

        print(sequence_prop_df.columns)
        # sequence_prop_df['normintegral_avg'] = list(sequence_prop_df.groupby('peptide')['normintegral'].agg('mean')
        print(sequence_prop_df.head(5))
        print(sequence_prop_df.sample(5))

        # sequence_prop_df = sequence_prop_df.sample(1000)  # sample 100 for testing
        propname = 'normintegral_avg'
        # propname = 'y_seq' # y_comp
        print(sequence_prop_df.sample(10))

        results_list = []
        for descriptor in ["peptides", "rdkit", "ecfp"]:
            print('\n\n', 'training models w/ descriptor:',
                  descriptor, 'conc', conc)

            for modelname in ["linear", "ridge", "lasso", "svr", "rf"]:
                print('training model', modelname,
                      'w/ descriptor:', descriptor)

                result = train_model(data=sequence_prop_df, propname=propname,
                                     descriptor=descriptor, modelname=modelname,
                                     conc=conc)

                results_list.append(result)

        results_df = pd.DataFrame(results_list)
        print(results_df)
        results_df.to_csv(f'./out/model_results_{conc}.csv', index=False)
