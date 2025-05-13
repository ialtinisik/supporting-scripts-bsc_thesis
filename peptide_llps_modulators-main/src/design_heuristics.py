def filter_synthesizable(seq):
    """ heuristics from literature to filter out non-synthesizble peptides
        each rule returns `True` if filter passes
        if any(*rule) == False --> not synthesizable
        adapted from doi:todo

        with rule: fraction_charged > 0.16 (~1/6):
        yields 16M (16265921) out of 64M hexapeptides ~ 25.41%
        {A: 1092427, C: 857686, D: 900855, E: 1070496, F: 826132, G: 1031610, H: 1070496,
         I: 826132, K: 1070496, L: 826132, M: 857686, N: 0, P: 832014, Q: 0,
         R: 1070496, S: 704418, T: 655143, V: 826132, W: 1092427, Y: 655143}
        """

    def no_forbidden_terminals(seq):
        ''' N-terminal Q is unstable due to cyclization under acidic conditions,
        N-terminal N should be avoided as its protecting group is hard to remove during cleavage. ''' 
        return False if seq[0] in ['Q', 'N'] else True

    def no_repeats(seq):
        ''' avoid consecutive repeats of same AA.
            repeat=True if any neighbors are same aa '''
        repeat = [True if seq[i] == seq[i+1] else False for i in range(len(seq)-1)]
        return False if any(repeat) else True

    def no_bad_aspartic_neighbors(seq):
        ''' Aspartic acid ('D') can undergo hydrolysis and cause peptide cleavage under
            acidic conditions when paired with glycine, proline or serine. '''
        # scan from 0 to len(-1) to avoid wrapping index with i+1
        bad_neighbor = [True if (seq[i] in ['G', 'P', 'S'] and seq[i+1] == 'D')
                             or (seq[i] == 'D' and seq[i+1] in ['G', 'P', 'S']) 
                        else False for i in range(len(seq)-1)]
        return False if any(bad_neighbor) else True

    def no_beta_sheets(seq):
        ''' A series of AAs "VIYFLQT" can form β-sheets, leading to incomplete solvation 
            during synthesis and causing deletions. '''
        beta = ['V', 'I', 'Y', 'F', 'L', 'Q', 'T']
        beta_neighbors = [True if (seq[i] in beta and seq[i+1] in beta) else False
                          for i in range(len(seq)-1)]
        return False if any(beta_neighbors) else True

    def no_more_than_one(seq):
        ''' Multiple of AAs "PSCM" residues in a sequence can cause
            significant deletions during synthesis '''
        # max_one_aa = ['P', 'S', 'C', 'M']
        # could also solve this with peptide.count(next_aa)
        counts = [sum([1 if seq[i] in max_aa else 0 for i, aa in enumerate(seq)])
                  for max_aa in ['P', 'S', 'C', 'M']]
        return False if any([c > 1 for c in counts]) else True

    def no_tryptophan_terminal(seq):
        ''' From our experience, peptides ending with tryptophan ('W')
            are often very low in (or have no) purity '''
        return False if seq[-1] == "W" else True

    def min_charged_fraction(seq):
        ''' A rule of thumb in designing soluble peptides is to ensure that at least
            1/5 AA residues is charged '''
        charged_count = sum([1 for aa in seq if aa in ['D', 'E', 'R', 'H', 'K']])
        return False if charged_count <= len(seq) * 0.16 else True

    def no_excess_crosslinking(seq):
        ''' Peptides with over 75% of AAs "DEHKNQRSTY" can form intermolecular hydrogen
            bonds (crosslinks), leading to gel formation in aqueous solutions. '''
        limited_count = sum([1 for aa in seq if aa in set("DEHKNQRSTY")])
        return False if limited_count >= len(seq) * 0.75 else True

    filters = [
        no_forbidden_terminals(seq),
        no_repeats(seq),
        no_bad_aspartic_neighbors(seq),
        no_beta_sheets(seq),
        no_more_than_one(seq),
        no_tryptophan_terminal(seq),
        min_charged_fraction(seq),
        no_excess_crosslinking(seq),
    ]
    # print(seq, 'passed' if all(filters) else 'filtered', '\t', filters)
    return seq if all(filters) else None

if __name__ == "__main__":
    valid_seqs = ['STEFANH', 'ARAKGK']
    assert all([filter_synthesizable(seq)for seq in valid_seqs])
    invalid_seqs = [
        'QRAKGK',  # Fails no_forbidden_terminals (starts with Q)
        'AARKGK',  # Fails no_repeats (two A's adjacent)
        'ARADGK',  # Fails no_bad_aspartic_neighbors (D is adjacent to G)
        'AVIKGK',  # Fails no_beta_sheets (V and I are adjacent, which are in beta-sheet forming set)
        'APRPSK',  # Fails no_more_than_one (P is repeated twice)
        'ARAKGW',  # Fails no_tryptophan_terminal (ends with W)
        'AFCDFG',  # Fails min_charged_fraction (only one charged residue, D)
        'DEKRNA',  # Fails no_excess_crosslinking (5/6 residues are in the crosslinking set)
    ]
    assert not all([filter_synthesizable(seq) for seq in invalid_seqs])
