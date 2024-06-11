ref_seqs = {
    0: ['GUAC', 'CCGG'],
    1: ['GCCG'],
    2: ['ACAGU'],
    3: ['GCACG'],
    4: ['CCCGGG'],
    5: [],
    6: [],
    7: ['UAACGGACUCU'],
    8: ['AAUCGGACUCU'],
    9: ['AAUCGGAUCU'],
    10: [],
    11: [],
    12: ['UUCAGAA'],
    13: [],
    14: [],
    15: [],
    16: [],
    17: [],
    # 18: [],
    # 19: [],
}

structs = {
    0: '(())',
    1: '()()',
    2: '((.))',
    3: '().()',
    4: '((()))',
    5: '(()())',
    6: '(())()',
    7: '()(()((.)))',
    8: '(()()((.)))',
    9: '(()()(()))',
    10: '()(.(()()))()',
    11: '()(..(().))',
    12: '(((.)))',
    13: '(.()()(.))',
    14: '.((...))',
    15: '(()((.)))',
    16: '(((()(.)))).',
    17: '..((.)(.))',
    # 18: '(((()()((()(.))))((.))))',
    # 19: '.(()())...((((()()))).())',
}

def is_valid_structure(structure):
    stack = []
    for char in structure:
        if char == '(':
            stack.append(char)
        elif char == ')':
            if not stack:
                return False
            stack.pop()
    return not stack

def generate_structures(n, pos=0, structure='', open_brackets=0):
    if pos == n:
        if is_valid_structure(structure):
            return [structure]
        return []

    structures = []
    if open_brackets > 0:
        structures += generate_structures(n, pos + 1, structure + ')', open_brackets - 1)
    structures += generate_structures(n, pos + 1, structure + '.', open_brackets)
    if pos < n - open_brackets:
        structures += generate_structures(n, pos + 1, structure + '(', open_brackets + 1)
    return structures

if __name__ == '__main__':
    n = 5
    structures = generate_structures(n)
    for s in structures:
        print(s)

    
