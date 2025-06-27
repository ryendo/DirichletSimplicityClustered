Of course. I will create a new version of the `README.md` file, incorporating your corrections and adding more technical detail and mathematical formulas from the paper to provide a deeper, more accurate description of the project.

Here is the revised `README.md`.

-----

# Computer-Assisted Proof for Dirichlet Eigenvalue Simplicity

[cite\_start]This project provides the source code and computational framework for the computer-assisted proof presented in the paper "Rigorous estimation for the difference quotients of multiple eigenvalues" (arXiv:2305.14063v5)[cite: 1]. [cite\_start]The primary goal is to rigorously validate the simplicity of the second Dirichlet eigenvalue for nearly equilateral triangles, offering a partial solution to Henrot's Conjecture 6.47[cite: 4, 15].

## Background

[cite\_start]Determining the eigenvalue multiplicity of the Laplace operator is a challenging problem, especially when eigenvalues are nearly degenerate, as is the case for the second and third Dirichlet eigenvalues on equilateral triangles where $\\lambda\_{2}=\\lambda\_{3}$[cite: 10]. Standard numerical methods struggle to separate these tightly clustered eigenvalues with mathematical rigor.

This work introduces a novel method to overcome this limitation. [cite\_start]Instead of relying on traditional shape derivatives, which are difficult to analyze over a neighborhood with a given radius[cite: 37], this project analyzes the **difference quotient of eigenvalues**, defined as:

[cite\_start]$$D_{t}\lambda_{i}:=\frac{\lambda_{i}(\Omega_{t})-\lambda(\Omega_{0})}{t}$$ [cite: 1]

By applying guaranteed computation techniques based on the Finite Element Method (FEM) and interval arithmetic, we can obtain rigorous bounds on this quotient.

The overall strategy is a two-pronged attack, dividing the parameter space of triangles $\\Omega$ into distinct regions:

  * [cite\_start]**For nearly equilateral triangles ($\\Omega\_{up}$):** The difference quotient method is used to prove that $D\\lambda\_{2}(p\_{0},p)\<D\\lambda\_{3}(p\_0, p)$ for any perturbation $p$ from the equilateral vertex $p\_0$[cite: 152]. [cite\_start]Since $\\lambda\_{2}^{p\_{0}}=\\lambda\_{3}^{p\_{0}}$[cite: 153], this rigorously implies $\\lambda\_{2}^{p}\<\\lambda\_{3}^{p}$.
  * [cite\_start]**For other triangles ($\\Omega\_{down}^{(1)}, \\Omega\_{down}^{(2)}$):** High-precision upper and lower bounds are computed directly to show a definitive gap, $\\overline{\\lambda}*{2}^{p}\<\\underline{\\lambda}*{3}^{p}$[cite: 153].

## Project Structure

The repository is organized as follows:

```
.
├── Each_Process/         # Core MATLAB functions, libraries, logs, and temporary results
│   ├── functions/        # Helper MATLAB functions for computation
│   ├── HighOrderFEM_CGYOU_2016/ # FEM library
│   ├── Intlab_Group/     # INTLAB library for guaranteed computation
│   └── log/              # Log files from parallel execution
├── prep/                 # Scripts for preparing the computation
├── results/              # Final stored results from computations
├── main_algo1.m          # Main MATLAB script for Algorithm 1
├── main_algo2.sh         # Main execution script for Algorithm 2
├── main_algo3.sh         # Execution script for Algorithm 3
├── prep.sh               # Master script to run all preparation steps
└── README.md             # This file
```

## Setup and Prerequisites

### Prerequisites

1.  **MATLAB:** A working installation of MATLAB is required.
2.  **INTLAB Library:** The INTLAB toolbox for interval arithmetic is essential for the guaranteed computations.
3.  **Unix-like Environment:** A shell like `bash` or `zsh` and standard command-line tools (`grep`, `awk`, `sed`, etc.) are needed to run the execution scripts.

### Initial Configuration

1.  **Clone the Repository:**

    ```bash
    git clone <your-repository-url>
    cd DirichletSimplicityClustered
    ```

2.  **Configure MATLAB Command-Line Access:**
    The shell scripts need to be able to call `matlab` from the command line. You must add the MATLAB binary directory to your shell's `PATH` environment variable.

      * First, find the path to your MATLAB executable using `which matlab` or `mdfind -name matlab | grep '/bin/matlab$'`. The path will be similar to `/Applications/MATLAB_R2024a.app/bin`.
      * Add this directory to your shell's configuration file.
          * For **bash**, add the following line to `~/.bash_profile`:
            ```bash
            export PATH="/path/to/your/matlab/bin:$PATH"
            ```
          * For **zsh**, add the same line to `~/.zshrc` or `~/.zprofile`.

3.  **Place the INTLAB Library:**
    Download and place the INTLAB library folder (named `Intlab_V12`) into the following directory:
    `Each_Process/Intlab_Group/`

## Execution Workflow

### 1\. Preparation

First, run the master preparation script from the project root directory.

```bash
./prep.sh
```

This command will:

  * Execute `prep/make_indices_algo2.m` to generate the task list `prep/list_j.csv`.
  * Execute `prep/prepare_parallel.sh` to create isolated INTLAB environments for each parallel worker.

### 2\. Main Computations

The proof is executed by running the different main scripts, each corresponding to a specific algorithm and region from the paper.

  * **Algorithm 1 (`main_algo1.m`):** Run from within MATLAB to analyze the nearly equilateral region $\\Omega\_{up}$.
    ```matlab
    >> main_algo1
    ```
  * **Algorithm 2 (`main_algo2.sh`):** Run from the shell to analyze the $\\Omega\_{down}^{(1)}$ region.
    ```bash
    bash ./main_algo2.sh
    ```
  * **Algorithm 3 (`main_algo3.sh`):** Run from the shell to analyze the $\\Omega\_{down}^{(2)}$ region.
    ```bash
    bash ./main_algo3.sh
    ```

### 3\. Monitoring the Process

For the parallel scripts (`main_algo2.sh`, `main_algo3.sh`), you can monitor the status:

  * **Check for running processes:** `ps aux | grep -i matlab`
  * **Watch log files in real-time:** `tail -f Each_Process/log/process_no1.log`
  * **Check for results:** Final results are written to CSV files in the `results/` directory.

## Code and Script Roles

### Root Directory

  * [cite\_start]`main_algo1.m`: Implements **Algorithm 1** from the paper[cite: 192]. This script analyzes the nearly equilateral region $\\Omega\_{up}$ by rigorously computing the difference quotients of $\\lambda\_2$ and $\\lambda\_3$. [cite\_start]It solves the generalized matrix eigenvalue problem $M\_{t}\\sigma=\\mu N\_{t}\\sigma$ [cite: 91] using interval arithmetic to obtain guaranteed bounds on the difference quotients.
  * [cite\_start]`main_algo2.sh`: The main parallel execution engine for **Algorithm 2**[cite: 208]. It analyzes the $\\Omega\_{down}^{(1)}$ region by partitioning it into small rectangles $R\_{ij}$ and applying perturbation estimates. It uses a robust worker pool model to distribute the computation of the eigenvalue bounds, which are governed by the relation:
    $$
    $$$$m\_{x}((x\_{i+1},y\_{j}),(x\_{i},y\_{j}))\\underline{\\lambda}*{k}^{(x*{i+1},y\_{j})}\\le\\lambda\_{k}^{p}\\le M\_{x}((x\_{i},y\_{j+1}),(x\_{i+1},y\_{j+1}))\\overline{\\lambda}*{k}^{(x*{i},y\_{j+1})}
    [cite\_start]$$ [cite: 207]
  * [cite\_start]`main_algo3.sh`: A parallel runner for **Algorithm 3**[cite: 221]. It analyzes the $\\Omega\_{down}^{(2)}$ region, where triangles are more degenerate. This algorithm leverages the domain monotonicity property of Dirichlet eigenvalues to establish bounds, using the relation:
    $$
    $$$$\\lambda\_{k}^{(r\_{i+1},h\_{j+1})}\\le\\lambda\_{k}^{(r,h)}\\le\\lambda\_{k}^{(r\_{i},h\_{j})}
    [cite\_start]$$ [cite: 220]

### Preparation Scripts (`prep/`)

  * [cite\_start]`make_indices_algo2.m`: Generates the master task list `prep/list_j.csv` for `main_algo2.sh`, containing the `j` indices (1-1220) for the computation over the $\\Omega\_{down}^{(1)}$ region[cite: 206].
  * `prepare_parallel.sh`: Sets up the execution environment by creating isolated copies of the INTLAB library for each parallel worker.
  * `run_matlab.sh`: A helper script called by `main_algo2.sh` to execute a single MATLAB computation task.

### Core MATLAB Code (`Each_Process/`)

  * [cite\_start]`func_algo2.m`: Implements the core logic for **Algorithm 2**[cite: 208]. It takes a list of `j` indices, iterates through the `i` index, and computes rigorous eigenvalue bounds for each subdomain `R_ij`.
  * [cite\_start]`func_algo3.m`: Implements the core logic for **Algorithm 3**[cite: 221], using domain monotonicity.
  * `my_intlab_mode_config.m`: A crucial function that initializes the INTLAB guaranteed computation environment for each specific worker process.
  * **`functions/` directory:** Contains various helper functions for the main algorithms.
      * [cite\_start]`calc_eigen_bounds_any_order.m`: A core function that calculates high-precision eigenvalue bounds, utilizing the Lehmann-Goerisch method as described in Lemma 4.3 [cite: 166] to achieve the required accuracy.
      * `Lagrange_...` functions: Related to the Lagrange Finite Element Method (FEM) spaces used for approximation.
      * `get_mesh_...` functions: Responsible for generating the specific triangulations (`\mathcal{T}^{h}`) of the domains needed for FEM computation.