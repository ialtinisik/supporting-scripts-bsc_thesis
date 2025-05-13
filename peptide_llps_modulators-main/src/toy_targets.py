import itertools
import re


def toy_target_composition(seq):
    ''' composition-based toy target:
        +1 for each AA in ['S', 'T', 'E', 'F', 'A', 'N'] '''
    counts = [sum([1 if seq[i] in max_aa else 0 for i, aa in enumerate(seq)])
              for max_aa in ['S', 'T', 'E', 'F', 'A', 'N']]
    return sum(counts)


def toy_target_sequence(seq):
    ''' sequence-based toy target: 1 if seq contains any permutation of
        'AT/NL/ES/IR' as substring ~ 1.2% of seqs '''
    # pattern = '|'.join([''.join(p) for p in itertools.permutations('AT')])
    pattern = '|'.join(['AT', 'TA', 'NL', 'LN', 'ES', 'SE', 'IR', 'RI'])
    return 1 if re.search(pattern, seq) else 0


if __name__ == "__main__":
    from pathlib import Path
    import pandas as pd

    seqs = ['AMRIEN', 'RMREI', 'WGMEIR', 'AAMERI', 'IREMLA',
            'QRAKGK', 'AARKGK', 'ARADGK', 'AVIKGK', 'APRPSK']
    [print(s, toy_target_composition(s), toy_target_sequence(s)) for s in seqs]

    all_seqs = []
    for file in list(Path('./data/seqs').glob('*.txt')):
        with open(file, 'r') as f:
            seqs = f.readlines()
            seqs = [seq.strip() for seq in seqs]
            all_seqs.extend(seqs)

    sample_seqs = pd.DataFrame(all_seqs, columns=['sequence']).sample(10000)
    sample_seqs['y_comp'] = sample_seqs['sequence'].apply(
        toy_target_composition)
    sample_seqs['y_seq'] = sample_seqs['sequence'].apply(
        toy_target_sequence)
    print('frac y_seq', sum(sample_seqs['y_seq']) / len(sample_seqs['y_seq']))
    sample_seqs.to_csv('./out/seqs_toy_targets.csv', index=False)
    print(sample_seqs.reset_index(drop=True).head(5))
