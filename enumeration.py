from collections import defaultdict
from heapq import heapify, heappop, heappush

import data as data

VALID_PAIRS = set(['AU', 'UA', 'CG', 'GC'])

def extract_base_pairs(struc, one_base_indexed=False):
    stack = []
    base_pairs = []
    for i, base in enumerate(struc):
        if base == '(':
            stack.append(i)
        elif base == ')':
            j = stack.pop()
            if one_base_indexed:
                base_pairs.append((j + 1, i + 1))
            else:
                base_pairs.append((j, i))
    return base_pairs


def kbest(s, k):  # using the division in total (no need to check duplicate), not in best
    left = lambda m: opt[i + 1, m - 1]  # i ( xxx ) m ___ j
    right = lambda m: opt[m, j]  # i ( ___ ) m xxx j
    single = lambda: opt[i + 1, j]  # i . xxxxxxxxxxx j

    def tryadd(split, index1, index2, func=heappush):  # i paired with split-1
        if index1 < len(left(split)) and index2 < len(right(split)) and not (split, index1, index2) in used:
            used.add((split, index1, index2))
            func(heap, (-left(split)[index1] - right(split)[index2] - 1, split, index1, index2))  # pair

    def tryadd1(index, func=heappush):
        if index < len(single()):
            func(heap, (-single()[index], index))  # min-heap, i unpaired

    def solution(i, j, index):
        if i == j:
            return ""  # empty
        if i == j - 1:
            return "."  # singleton
        split, index1, index2 = back[i, j][index]
        if split is None:
            return ".%s" % solution(i + 1, j, index1)
        return "(%s)%s" % (solution(i + 1, split - 1, index1), solution(split, j, index2))

    n = len(s)
    opt = defaultdict(list)  # list of values
    back = defaultdict(list)
    for i in range(n + 1):  # N.B. opt[n,n]=0
        opt[i, i + 1] = [0]  # singleton
        opt[i, i] = [0]  # empty

    for span in range(2, n + 1):  # span length: 2..n
        for i in range(n - span + 1):  # left boundary-: 0..(n-span)
            j = i + span  # right boundary

            heap, used = [], set()

            tryadd1(0, func=list.append)  # i unpaired  # don't heappush
            for m in range(i + 1, j + 1):  # i pair with some k
                if s[i] + s[m - 1] in VALID_PAIRS:  # i( ___ )k ___ j
                    tryadd(m, 0, 0, func=list.append)  # dont' heappush

            heapify(heap)

            while heap != [] and len(opt[i, j]) < k:
                item = heappop(heap)
                if len(item) == 2:  # i unpaired
                    value, index = item
                    opt[i, j].append(-value)
                    back[i, j].append((None, index, None))
                    tryadd1(index + 1)
                else:  # i paired, split
                    value, split, index1, index2 = item
                    opt[i, j].append(-value)
                    back[i, j].append((split, index1, index2))
                    tryadd(split, index1 + 1, index2)
                    tryadd(split, index1, index2 + 1)

    return [(opt[0, n][p], solution(0, n, p)) for p in range(min(k, len(opt[0, n])))]

# Time complexity: O(4^(n-p)), where p is the number of paired bases and n is the length of the sequence
def get_seqs(structure):
    length, base_pairs = len(structure), extract_base_pairs(structure)
    wc_pairs = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
    valid_seqs = []

    def is_valid(seq, base_pairs):
        for (i, j) in base_pairs:
            if seq[i] != wc_pairs.get(seq[j], ''):
                return False
        return True

    def backtrack(seq, position, base_pairs):
        if position == length:
            if is_valid(seq, base_pairs):
                valid_seqs.append(''.join(seq))
            else:
                # throw error
                raise ValueError('Invalid seq')
            return

        if seq[position] == '':
            paired_pos = None
            for (i, j) in base_pairs:
                if i == position:
                    paired_pos = j
                    break
                elif j == position:
                    paired_pos = i
                    break
            
            if paired_pos is not None and seq[paired_pos] == '':
                for base, pair in wc_pairs.items():
                    seq[position] = base
                    seq[paired_pos] = pair
                    backtrack(seq, position + 1, base_pairs)
                    seq[position] = ''
                    seq[paired_pos] = ''
            else:
                for base in 'AUCG':
                    seq[position] = base
                    backtrack(seq, position + 1, base_pairs)
                    seq[position] = ''
        else:
            backtrack(seq, position + 1, base_pairs)

    seq = [''] * length
    backtrack(seq, 0, base_pairs)
    return valid_seqs

def analyze_design(struc, seq):
    base_pairs = set(extract_base_pairs(struc))
    pidx = {}   # pair index of each base if it is paired
    paired_pos = set()
    unpaired_pos = set(range(len(seq)))
    for (i, j) in base_pairs:
        pidx[i] = j
        pidx[j] = i
        paired_pos.add(i)
        paired_pos.add(j)
        unpaired_pos.remove(i)
        unpaired_pos.remove(j)

    def val_pair(i, j):
        return seq[i] + seq[j] in VALID_PAIRS
    
    psuedo_pairs = defaultdict(list)
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if val_pair(i, j) and (i, j) not in base_pairs:
                for k in range(i + 1, j):
                    if struc[k] != '.' and (pidx[k] < i or pidx[k] > j):
                        psuedo_pairs[i, j].append((min(k, pidx[k]), max(k, pidx[k])))

    def get_broken_nuc_pos(i, j):
        i, j = min(i, j), max(i, j)
        broken_pos, inside, outside = set(), [], []
        broken_pairs = psuedo_pairs[i, j]
        for (p, q) in set(broken_pairs):
            broken_pos.add(p)
            broken_pos.add(q)
        if struc[i] != '.':
            broken_pos.add(i)
            broken_pos.add(pidx[i])
        if struc[j] != '.':
            broken_pos.add(j)
            broken_pos.add(pidx[j])
        for bpos in broken_pos:
            if bpos >= i and bpos <= j:
                inside.append(bpos) 
            else:
                outside.append(bpos) 
        return inside, outside
    
    def get_impssbl_pos(positions):
        no_pair_pos = []
        for i in positions:
            pair_found = False
            for j in positions:
                if val_pair(i, j):
                    pair_found = True
                    break
            if not pair_found:
                no_pair_pos.append(i)
        return no_pair_pos
    
    
    total_distance = len(seq)
    wrong_pairs = []
    for i in range(len(seq)):
        if struc[i] == '.':  # unpaired
            for j in range(0, len(seq)):
                if val_pair(i, j) and (i, j) not in base_pairs:
                    in_brk_pos, out_brk_pos = get_broken_nuc_pos(i, j)
                    in_unp_pos = set(in_brk_pos) - {i, j}
                    out_unp_pos = set(out_brk_pos)

                    min_break_count = 2
                    if struc[j] == '.':
                        min_break_count += 1

                    for q in unpaired_pos - {i, j}:
                        if q > min(i, j) and q < max(i, j):     # q inside the pair
                            in_unp_pos.add(q)
                            min_break_count += 1
                        else:   # q outside the pair
                            for p in out_brk_pos:
                                if len(set(get_broken_nuc_pos(min(p, q), max(p, q))[0]) - set(in_brk_pos)) <= 1:
                                    out_unp_pos.add(q)
                                    min_break_count += 1
                                    break

                    if len(get_impssbl_pos(in_unp_pos)) + len(get_impssbl_pos(out_unp_pos)) < min_break_count:
                        total_distance -= 1
                        wrong_pairs.append((i + 1, j + 1))
                
        else:  # paired
            for j in range(i + 1, len(seq)):
                if val_pair(i, j) and (i, j) not in base_pairs:
                    valid = False
                    for x in range(i + 1, j):
                        if struc[x] != '.' and (pidx[x] <= i or pidx[x] >= j):
                            valid = True
                    if not valid:
                        total_distance -= 1
                        wrong_pairs.append((i + 1, j + 1))


    return total_distance, wrong_pairs


def get_designs(structure, verbose=False):
    if verbose:
        print('Designing structure %s' % structure)
    designs = []
    for seq in get_seqs(structure):
        assert len(seq) == len(structure), 'Incompatible sequence %s for structure %s' % (seq, structure)
        if analyze_design(structure, seq)[0] == len(seq):
            designs.append(seq)
            if verbose:
                print(seq)
    if len(designs) == 0 and verbose:
        print('No valid designs found!')
    else:
        print('Found %d designs for structure %s' % (len(designs), structure))
    return designs

# used to verify the solution
def verify(structure, ref_seqs = [], verbose=False):
    success = True
    invalid_designs = []
    seqs = get_designs(structure, verbose=False)
    for seq in seqs:
        two_best_structs = kbest(seq, 2)  # get top 2 structures
        if len(two_best_structs) < 2:
            continue
        if two_best_structs[0][0] == two_best_structs[1][0]:  # if top 2 structures have the same energy, then its 
                                                              # not unique, therefore invalid design
            invalid_designs.append((seq, two_best_structs))
            success = False

    if verbose:
        print('Verifying designs...')
        if success:
            print('[SUCCESS] Designing Successful! All folded structures are unique! :)')
            if len(seqs) == 0:
                print('[INFO] No valid designs found!')
        else:
            print('[FAILURE] Designing Failed! Some folded structures are not unique! :(')
            for design, structs in invalid_designs:
                print('%s %s' % (design, structs))
        if ref_seqs:
            for ref_seq in ref_seqs:
                if ref_seq in seqs:
                    print("[SUCCESS] Reference sequence '%s' found in the designs" % ref_seq)
                else:
                    print("[FAILURE] Reference sequence '%s' not found in the designs" % ref_seq)
                    if ref_seq in get_seqs(structure):
                        print("[INFO] Reference sequence is in the enumeration of compatible sequences")
                    else:
                        print("[INFO] Reference sequence is not in the enumeration of compatible sequences")

    if len(seqs) == 0:
        print('[INFO] No valid designs found!')
        seqs = get_seqs(structure)
        num_pairs = len(extract_base_pairs(structure))
        for seq in seqs:
            out = kbest(seq, 10)            
            while out and out[0][0] > num_pairs:
                out.pop(0)
            if len(out) > 1 and out[0][0] != out[1][0]:
                print('heres the deign:', seq, out)

    return success, invalid_designs
    

if __name__ == '__main__':
    for tidx in range(0, len(data.structs)):
        verify(data.structs[tidx], data.ref_seqs[tidx], verbose=True)
        print('\n', end='')

    # for struc in data.generate_structures(10):
    #     verify(struc, verbose=False)
    #     print('\n', end='')
    
    # for debugging ------------------------------------------------
    # print(analyze_design('..((.)(.))', 'AACACUUGAG'))
    # print(analyze_design('()(..(().))', 'AUUCCUUAGAA'))
    # print(analyze_design('().()', 'AUUUA'))
