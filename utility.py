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
