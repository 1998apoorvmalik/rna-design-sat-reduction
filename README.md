# RNA Design to SAT Reduction

This repository contains our project work for reducing the RNA design problem to SAT. It uses `pysat` (SAT solver) to find the boolean variables assignment which corresponds to the solution of the RNA design problem, if it exists. The solution to the RNA design problem is a sequence of nucleotides which uniquely folds into the input secondary structure. Any alternative structure that the sequence can fold into must contain strictly fewer base pairs than the input secondary structure.

The file `enumeration.py` contains the code for solving the RNA design problem using the naive enumeration method. \
The file `sat.py` contains the code for reducing the RNA design problem to SAT and then solving it.

## Web Server

We have hosted our RNA Design SAT Solver at [https://rna-design-web-server.uw.r.appspot.com](https://rna-design-web-server.uw.r.appspot.com). You can use this web server to design your own RNA sequences for a target secondary structure.

Alternatively, you can run the code on your local machine by following the instructions below.

## Running the Web Server Locally

### Prerequisites

- Python 3.6+
- `pip` (Python package installer)
- `pysat` library
- `python-sat` library
- `Flask` web framework

You also need to have your pysat configured, to be ready to use.

### Installation

1. Clone the repository to your local machine:

   ```sh
   git clone https://github.com/1998apoorvmalik/rna-design-sat-reduction
   cd rna-design-sat-reduction
   ```

2. [Optional] Create and activate a virtual environment:

   ```sh
   python3 -m venv venv
   source venv/bin/activate
   ```

3. Install the required packages:

   ```sh
   pip install -r requirements.txt
   ```

### Running the Application

1. Run the Flask development server:

   ```sh
   cd web-app
   python app.py
   ```

2. Open your web browser and navigate to:

   ```
   http://127.0.0.1:5000/
   ```

## Contributing

We welcome contributions to improve the project. Please submit a pull request or open an issue to discuss changes.

## Contact

If you have any questions or need further assistance, feel free to contact us at [malikap@oregonstate.edu], [ameyass@oregonstate.edu].
