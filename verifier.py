from collections import defaultdict
from heapq import heappush, heappop, heapify

from utility import extract_base_pairs, VALID_PAIRS

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


# used to verify the solution
def verify(structure, seqs, ref_seqs = [], verbose=False):
    print('Verifying %s designs for structure %s' % (len(seqs), structure))
    # if verbose:
    #     for seq in seqs:
    #         print(''.join(seq))
    success = True
    invalid_designs = []
    for seq in seqs:
        two_best_structs = kbest(seq, 2)  # get top 2 structures
        if len(two_best_structs) < 2:
            continue
        if two_best_structs[0][0] == two_best_structs[1][0]:  # if top 2 structures have the same energy, then its 
                                                              # not unique, therefore invalid design
            invalid_designs.append((seq, two_best_structs))
            success = False

    if verbose:
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

    # if verbose and len(seqs) == 0:
    #     print('[INFO] No valid designs found!')
    #     seqs = get_seqs(structure)
    #     num_pairs = len(extract_base_pairs(structure))
    #     for seq in seqs:
    #         out = kbest(seq, 10)            
    #         while out and out[0][0] > num_pairs:
    #             out.pop(0)
    #         if len(out) > 1 and out[0][0] != out[1][0]:
    #             print('heres the deign:', seq, out)

    return success, invalid_designs