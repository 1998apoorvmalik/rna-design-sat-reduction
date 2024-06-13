from pysat.solvers import Solver
from pysat.formula import CNF

from itertools import combinations
from collections import defaultdict


from utility import extract_base_pairs
from verifier import verify

import data as data
import random
import time
import matplotlib.pyplot as plt

def encode_rna_design(struc):
    target_pairs = set(extract_base_pairs(struc))

    pidx = {}   # pair index of each base if it is paired
    paired_pos = set()
    unpaired_pos = set(range(len(struc)))
    for (i, j) in target_pairs:
        pidx[i] = j
        pidx[j] = i
        paired_pos.add(i)
        paired_pos.add(j)
        unpaired_pos.remove(i)
        unpaired_pos.remove(j)

    psuedo_pairs = defaultdict(list)
    for i in range(len(struc)):
        for j in range(i + 1, len(struc)):
            if (i, j) not in target_pairs:
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
        return set(inside) - {i, j}, set(outside)
    
    cnf = CNF()
    
    # Define base variables
    def b(i, base):
        bases = {'A': 0, 'U': 1, 'G': 2, 'C': 3}
        return i * 4 + bases[base] + 1
    
    def add_non_pairing_constraint(i, j):
        cnf.append([-b(i, 'A'), -b(j, 'U')])
        cnf.append([-b(i, 'U'), -b(j, 'A')])
        cnf.append([-b(i, 'G'), -b(j, 'C')])
        cnf.append([-b(i, 'C'), -b(j, 'G')])

    # Add constraints to ensure each position has exactly one base
    for i in range(len(struc)):
        cnf.append([b(i, 'A'), b(i, 'U'), b(i, 'G'), b(i, 'C')])  # At least one base
        cnf.append([-b(i, 'A'), -b(i, 'U')])
        cnf.append([-b(i, 'A'), -b(i, 'G')])
        cnf.append([-b(i, 'A'), -b(i, 'C')])
        cnf.append([-b(i, 'U'), -b(i, 'G')])
        cnf.append([-b(i, 'U'), -b(i, 'C')])
        cnf.append([-b(i, 'G'), -b(i, 'C')])

    # Add Watson-Crick pairing constraints for target pairs
    for (i, j) in target_pairs:
        cnf.append([-b(i, 'A'), b(j, 'U')])
        cnf.append([-b(i, 'U'), b(j, 'A')])
        cnf.append([-b(i, 'G'), b(j, 'C')])
        cnf.append([-b(i, 'C'), b(j, 'G')])
    
    for i in range(len(struc)):
        for j in range(i + 1, len(struc)):
            if (i, j) not in target_pairs:
                        
                in_bpos, out_bpos = get_broken_nuc_pos(i, j)

                min_brk_count = 1 + int(struc[i] == '.' and struc[j] == '.')
                
                if len(set(in_bpos)) < min_brk_count or len(set(out_bpos)) < min_brk_count:
                    add_non_pairing_constraint(i, j)
                elif struc[i] == '.' or struc[j] == '.':

                    # for p, q in combinations(in_bpos, 2):
                    #     if (p, q) not in target_pairs:
                    #         add_non_pairing_constraint(p, q)

                    # for p, q in combinations(out_bpos, 2):
                    #     if (p, q) not in target_pairs:
                    #         add_non_pairing_constraint(p, q)


                    for q in unpaired_pos:
                        if q == i or q == j:
                            continue

                        brk_pos = in_bpos if q > i and q < j else out_bpos
                        for p in brk_pos:
                            tmp = set()
                            for m in [i, j, p, q]:
                                if struc[m] != '.':
                                    tmp.add(m)
                                    tmp.add(pidx[m])

                                    
                            if (q > i and q < j) and len(set(in_bpos) - set(get_broken_nuc_pos(p, q)[0]) - tmp) == 0:
                                add_non_pairing_constraint(p, q)

                            if (q < i or q > j) and len(set(get_broken_nuc_pos(p, q)[0]) - set(in_bpos) - tmp) == 0:
                                add_non_pairing_constraint(p, q)
    return cnf

def solve_rna_design(struc, max_seq=float('inf')):
    cnf = encode_rna_design(struc)
    solutions = []

    with Solver(bootstrap_with=cnf) as solver:
        while solver.solve() and len(solutions) < max_seq:
            model = solver.get_model()
            base_assignment = [''] * len(struc)
            bases = ['A', 'U', 'G', 'C']
            
            for i in range(len(struc)):
                for b_idx, b in enumerate(bases):
                    base_var = i * 4 + b_idx + 1
                    if base_var in model:
                        base_assignment[i] = b
                        break

            solutions.append(''.join(base_assignment))
            
            # Add a blocking clause to prevent the solver from finding the same solution again
            blocking_clause = [-v for v in model if v > 0]
            solver.add_clause(blocking_clause)

    return solutions


if __name__ == '__main__':
    # Example usage
    # struc = '()(()())((())(())()())'
    # designs = solve_rna_design(struc)
    # verify(struc, designs, verbose=True)
    # print("All RNA Sequences:")
    # for solution in all_solutions:
    #     print("".join(solution))

    ## test case 2:
    for tidx in range(0, len(data.structs)):
        struc = data.structs[tidx]
        ref_seqs = data.ref_seqs[tidx]
        designs = solve_rna_design(struc)
        verify(struc, designs, ref_seqs, verbose=True)
        print('\n', end='')



    # Initialize the dictionary to store designable structures
    # designable_structures = defaultdict(list)
    # lengths = []
    # dots = []
    # times = []
    # data_size = 10

    # # Generate structures and collect data
    # for m in range(0, 7, 2):
    #     while len(designable_structures[m]) < data_size:        
    #         n = random.randint(10, 201)
    #         if n % 2:
    #             n += 1
    #         start_time = time.time()
    #         structure = data.generate_random_structure(n, m)
    #         designs = solve_rna_design(structure)
    #         end_time = time.time()
    #         if len(designs):
    #             if verify(structure, designs, verbose=False):
    #                 designable_structures[m].append((structure, designs[0]))
    #                 lengths.append(n)
    #                 dots.append(m)
    #                 times.append(end_time - start_time)
    #             else:
    #                 print("[FATAL ERROR]: Design verification failed!")

    # import numpy as np
    # import matplotlib.pyplot as plt
    # from scipy.optimize import curve_fit

    # # set font size
    # plt.rcParams.update({'font.size': 14})

    # # Function to fit the data
    # def fit_function(x, a, b):
    #     return a * np.power(x, b)

    # # Function to format the fitted curve as Big O notation
    # def big_o_notation(a, b):
    #     if b == 1:
    #         return r"$O(n)$"
    #     elif b.is_integer():
    #         return fr"$O(n^{{{int(b)}}})$"
    #     else:
    #         return fr"$O(n^{{{b:.2f}}})$"
        

    # # Aggregate all data points for length vs time
    # x_data = np.array(lengths)
    # y_data = np.array([time * 1000 for time in times])  # convert to milliseconds

    # # Fit the curve
    # popt, _ = curve_fit(fit_function, x_data, y_data)

    # # Generate curve data for plotting
    # x_fit = np.linspace(min(x_data), max(x_data), 200)
    # y_fit = fit_function(x_fit, *popt)

    # # Format the Big O notation
    # big_o = big_o_notation(*popt)

    # plt.figure(figsize=(12, 8))
    # plt.scatter(x_data, y_data, label='$y$ (structure)')
    # plt.plot(x_fit, y_fit, color='red', label=f'{big_o}')
    # plt.title('RNA Design: Time vs Sequence Length')
    # plt.xlabel('Sequence Length ($n$)')
    # plt.ylabel('Time ($ms$)')
    # plt.grid(True)
    # plt.legend()

    # # Save the plot
    # plt.savefig('time_vs_length.png')
    # plt.close()

    # # Plotting time vs length (all structures have same dots)
    # unique_dots = set(dots)
    # for dot in unique_dots:
    #     dot_indices = [i for i, d in enumerate(dots) if d == dot]
    #     x_data = np.array([lengths[i] for i in dot_indices])
    #     y_data = np.array([times[i] * 1000 for i in dot_indices])  # convert to milliseconds
        
    #     # Fit the curve
    #     popt, _ = curve_fit(fit_function, x_data, y_data)
        
    #     # Generate curve data for plotting
    #     x_fit = np.linspace(min(x_data), max(x_data), 100)
    #     y_fit = fit_function(x_fit, *popt)
        
    #     # Format the Big O notation
    #     big_o = big_o_notation(*popt)
        
    #     plt.figure(figsize=(12, 8))
    #     plt.scatter(x_data, y_data, label='$y$ (structure)')
    #     plt.plot(x_fit, y_fit, color='red', label=f'{big_o}')
    #     plt.title(f'RNA Design: Time vs Structure Length ({dot} Unpaired positions)')
    #     plt.xlabel('Length ($n$)')
    #     plt.ylabel('Time ($ms$)')
    #     plt.grid(True)
    #     plt.legend()

    #     # save the plot
    #     plt.savefig(f'time_vs_length_{dot}_dots.png')
        
    #     # plt.show()
