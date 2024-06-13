import os
import sys
from flask import Flask, render_template, request

# Add the parent directory to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sat import solve_rna_design
from verifier import kbest
from data import structs

app = Flask(__name__)

def get_rna_design(structure, num_sequences, verify_design):
    designs = solve_rna_design(structure, num_sequences)
    output = []
    for seq in designs:
        two_best_structs = kbest(seq, 2)  # get top 2 structures
        if verify_design:
            output.append([seq, not two_best_structs[0][0] == two_best_structs[1][0]])
        else:
            output.append([seq])
    return output

@app.route('/', methods=['GET', 'POST'])
def index():
    samples = [
        structs[10],
        structs[16],
        structs[8],
        "((..((....))..))",
        '(((()()((()(.))))((.))))',
        # "....((..))....",
    ]

    designed_sequences = None
    dot_bracket = ""
    num_sequences = 1
    verify_design = False
    if request.method == 'POST':
        dot_bracket = request.form['dot_bracket']
        num_sequences = int(request.form['num_sequences'])
        verify_design = 'verify_design' in request.form
        designed_sequences = get_rna_design(dot_bracket, num_sequences, verify_design)
        print(designed_sequences)
    return render_template('index.html', designed_sequences=designed_sequences,
                           dot_bracket=dot_bracket, num_sequences=num_sequences,
                           verify_design=verify_design, samples=samples)

if __name__ == '__main__':
    app.run(debug=True)
