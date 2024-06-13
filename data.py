import random

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

def generate_random_structure(n, m):
    if m > n or (n - m) % 2 != 0:
        raise ValueError("Invalid parameters: k should be less than or equal to n and n-k must be even")

    def is_valid_structure(structure):
        balance = 0
        for char in structure:
            if char == '(':
                balance += 1
            elif char == ')':
                balance -= 1
            if balance < 0:
                return False
        return balance == 0

    num_pairs = (n - m) // 2
    structure = ['.'] * m + ['('] * num_pairs + [')'] * num_pairs

    while True:
        random.shuffle(structure)
        if is_valid_structure(structure):
            return ''.join(structure)

if __name__ == '__main__':
    n = 10
    structures = generate_structures(n)
    print(len(structures))
    # for s in structures:
    #     print(s)

    # for i in range(4, 25, 2):
    #     structure = generate_random_structure(i,  2)
    #     print(structure)

    
