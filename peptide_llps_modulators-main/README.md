# peptide_llps_modulators
Irem's project - LLPS modulation peptide discovery 

## repository structure
```
├── src
│     ├── design_heuristics.py          # filters sequence based on synthesizability heuristics
│     ├── enumerator.py                 # enumerates all sequences up to `len` recursively (len=6)
│     ├── toy_targets.py                # creates toy target 'labels' to test model pipeline
│     └── todo.py
│
├── scripts
│     ├── CSCanalysisscript.py          # 
│     ├── train_model.py                # trains all [descriptors] * [models] on 
│     ├── todo.py                       # 
│
├── data
│   ├── csc_integrals_sequences_df.csv  # processed df: seq, csc
│   ├── seqs                            # enumerated and filtered hexapeptides (~16M)
```


## conda environment
```
conda create -n ml python=3.12
conda activate ml
pip install -r requirements.txt
export PYTHONPATH="./peptide_llps_modulators:$PYTHONPATH"
```

# pipeline
## descriptors
- peptides package
- rdkit Descriptors (rdkit.Chem.Descriptors.CalcMolDescriptors)
- structural fingerprint (ECFP/Morgan)

## models
- SupportVectorRegressor
- RandomForestRegressor
- Ridge/Lasso

## CSC 
measure turbditiy at increasing salt concentration using nephelometry 
calculate critical salt concentration from this curve


# workflow
### experimental
- synthesis (biotage)
- purification (ether precipitation)
- formulation (OT-2 flex)
- nephelometry (w/ salt titration)

### modeling
- calculate turbidity to CSC
- match measurement to input sequence --> df
- fit + evaluate models
- select best model+hyperparameters and fit full data
- filter all possible hexapeptides for synthesizability heuristics [`src.design_heuristics.py`]
- predict CSC for all sequences, select top-k ~ 192
- synthesize those and repeat!
