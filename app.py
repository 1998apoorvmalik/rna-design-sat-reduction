from flask import Flask, render_template, request
import sat as Solver
from verifier import kbest

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html', designs=[])

@app.route('/submit_form', methods=['POST'])
def submit_form():
    input_struc = request.form['input']
    designs = Solver.solve_rna_design(input_struc)
    design_results = []

    for des in designs:
        two_best_structs = kbest(des, 2)
        if len(two_best_structs) < 2:
            design_results.append({'design': des, 'result': 'Valid Design'})
        elif two_best_structs[0][0] == two_best_structs[1][0]:
            design_results.append({'design': des, 'result': 'Invalid Design'})
        else:
            design_results.append({'design': des, 'result': 'Valid Design'})

    return render_template('index.html', designs=design_results)

if __name__ == '__main__':
    app.run(debug=True)
