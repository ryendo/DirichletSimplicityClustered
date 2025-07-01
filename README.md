# Computer-Assisted Proof for Dirichlet Eigenvalue Simplicity

This project provides the source code and computational framework for the computer-assisted proof presented in the paper "Rigorous estimation for the difference quotients of multiple eigenvalues". The primary goal is to rigorously validate the simplicity of the second Dirichlet eigenvalue for nearly equilateral triangles, offering a partial solution to a conjecture posed by R. Laugesen and B. Siudeja (discussed as Conjecture 6.47 in A. Henrot's "Shape Optimization and Spectral Theory").

---

## Background

Determining the eigenvalue multiplicity of the Laplace operator is a challenging problem, especially when eigenvalues are nearly degenerate, as is the case for the second and third Dirichlet eigenvalues on equilateral triangles where $\lambda_{2}=\lambda_{3}$. Standard numerical methods struggle to separate these tightly clustered eigenvalues with mathematical rigor.

This work introduces a novel method to overcome this limitation. Instead of relying on traditional shape derivatives, which are difficult to analyze over a neighborhood with a given radius, this project analyzes the **difference quotient of eigenvalues**, defined as:

$$
D_{t}\lambda_{i}:=\frac{\lambda_{i}(\Omega_{t})-\lambda(\Omega_{0})}{t}
$$

By applying guaranteed computation techniques based on the Finite Element Method (FEM) and interval arithmetic, we can obtain rigorous bounds on this quotient.

The overall strategy is a two-pronged attack, dividing the parameter space of triangles $\Omega$ into distinct regions:

* For nearly equilateral triangles ($\Omega_{up}$), the difference quotient method is used to prove that $D\lambda_{2}(p_{0},p)<D\lambda_{3}(p_0, p)$ for any perturbation $p$ from the equilateral vertex $p_0$. Since $\lambda_{2}^{p_{0}}=\lambda_{3}^{p_{0}}$, this rigorously implies $\lambda_{2}^{p}<\lambda_{3}^{p}$.
* For other triangles ($\Omega_{down}^{(1)}, \Omega_{down}^{(2)}$), high-precision upper and lower bounds are computed directly to show a definitive gap, $\overline{\lambda}*{2}^{p}<\underline{\lambda}*{3}^{p}$.

---

## Project Structure

The repository is organized as follows:

```
.
├── Each_Process/                   # Core MATLAB functions, logs, and temporary results for parallel jobs
│   ├── FEM_Functions/              # Helper MATLAB functions for FEM computation
│   ├── Intlab_Group/               # INTLAB library for guaranteed computation
│   ├── log/                        # Log files from parallel execution
│   └── ...
├── HighOrderFEM_CGYOU_2016/        # FEM library
├── HighOrderFEM_CGYOU_2016_Dirichlet/ # FEM library for Dirichlet problems
├── prep/                           # Scripts for preparing the computation
├── results/                        # Final stored results from computations
├── verified_eig_estimation/        # Scripts for eigenvalue verification
├── main_algo1.m                    # Main MATLAB script for Algorithm 1
├── main_algo2.sh                   # Main execution script for Algorithm 2
├── main_algo3.sh                   # Execution script for Algorithm 3
├── my_intlab_config_alone.m        # Standalone INTLAB configuration
├── prep.sh                         # Master script to run all preparation steps
└── README.md                       # This file
```

---

## Setup and Prerequisites

### Prerequisites

1.  **MATLAB:** A working installation of MATLAB is required.
2.  **INTLAB Library:** The INTLAB toolbox for interval arithmetic is essential for the guaranteed computations.
3.  **Unix-like Environment:** A shell like `bash` or `zsh` and standard command-line tools are needed to run the parallel execution scripts.

### Initial Configuration

1.  **Clone the Repository:**
    ```bash
    git clone https://github.com/ryendo/DirichletSimplicityClustered
    cd DirichletSimplicityClustered
    ```

2.  **Configure MATLAB Command-Line Access:**
    The shell scripts need to be able to call `matlab` from the command line. You must add the MATLAB binary directory to your shell's `PATH` environment variable.

    * First, find the path to your MATLAB executable (e.g., using `which matlab`). The path will be similar to `/Applications/MATLAB_R2024a.app/bin`.
    * Add this directory to your shell's configuration file (`~/.bash_profile` for bash, `~/.zshrc` or `~/.zprofile` for zsh).
        ```bash
        export PATH="/path/to/your/matlab/bin:$PATH"
        ```

3.  **Place the INTLAB Library:**
    Download and place the INTLAB library folder (named `Intlab_V12`) into the following directory:
    `Each_Process/Intlab_Group/`

---

## Execution Workflow

The proof is executed by running the different main scripts, each corresponding to a specific algorithm and region from the paper.

### Algorithm 1 (Difference Quotient Analysis)

This algorithm is run as a standard MATLAB script and does not use the parallel execution framework.

1.  **Initialize INTLAB:** First, run the standalone configuration script from the MATLAB command window.
    ```matlab
    >> my_intlab_config_alone
    ```
2.  **Run the Algorithm:** After the setup is complete, run the main script.
    ```matlab
    >> main_algo1
    ```

### Algorithms 2 & 3 (Parallel Eigenvalue Bounds)

These algorithms are run from the shell and utilize a parallel framework.

1.  **Configuration (One-time setup):**

    * **Set Job Count:** Before running, you must configure the number of parallel jobs. Open `main_algo2.sh`, `main_algo3.sh`, and `prep/prepare_parallel.sh` and set the `MAX_JOBS` variable at the top of each file to your desired number of cores. For example:
        ```bash
        MAX_JOBS=60
        ```
    * **Set Permissions:** Grant execute permission to the MATLAB runner script. This only needs to be done once.
        ```bash
        chmod +x prep/run_matlab.sh
        ```

2.  **Preparation:** Run the master preparation script from the project root directory. This script uses the `MAX_JOBS` setting from `prep/prepare_parallel.sh`.

    ```bash
    ./prep.sh
    ```

    This command will:

    * Execute `prep/make_indices_algo2.m` to generate the task list `prep/list_j.csv`.
    * Execute `prep/prepare_parallel.sh` to create isolated INTLAB environments for each parallel worker.

3.  **Main Computations:**

    * To analyze the $\Omega_{down}^{(1)}$ region (Algorithm 2):
        ```bash
        bash ./main_algo2.sh
        ```
    * To analyze the $\Omega_{down}^{(2)}$ region (Algorithm 3):
        ```bash
        bash ./main_algo3.sh
        ```

4.  **Monitoring the Process:**
    For the parallel scripts, you can monitor the status:

    * **Check for running processes:** `ps aux | grep -i matlab`
    * **Watch log files in real-time:** `tail -f Each_Process/log/process_no1.log`
    * **Check for results:** Final results are written to CSV files in the `results/` directory.

---

## Code and Script Roles

### Root Directory

* `main_algo1.m`: Implements **Algorithm 1** from the paper. This script analyzes the nearly equilateral region $\Omega_{up}$ by rigorously computing the difference quotients of $\lambda_2$ and $\lambda_3$. It solves the generalized matrix eigenvalue problem $M_{t}\sigma=\mu N_{t}\sigma$ using interval arithmetic to obtain guaranteed bounds on the difference quotients.
* `main_algo2.sh`: The main parallel execution engine for **Algorithm 2**. It analyzes the $\Omega_{down}^{(1)}$ region by partitioning it into small rectangles $R_{ij}$ and applying perturbation estimates. It uses a robust worker pool model to distribute the computation of the eigenvalue bounds.
* `main_algo3.sh`: A parallel runner for **Algorithm 3**. It analyzes the $\Omega_{down}^{(2)}$ region, where triangles are more degenerate. This algorithm leverages the domain monotonicity property of Dirichlet eigenvalues to establish bounds.
* `my_intlab_config_alone.m`: A setup script to initialize the INTLAB guaranteed computation environment for a single, non-parallel MATLAB session. It must be run in MATLAB before executing `main_algo1.m`.
* `HighOrderFEM_CGYOU_2016/`, `HighOrderFEM_CGYOU_2016_Dirichlet/`: Finite Element Method (FEM) libraries used for the numerical computations.

### Preparation Scripts (`prep/`)

* `make_indices_algo2.m`: Generates the master task list `prep/list_j.csv` for `main_algo2.sh`, containing the `j` indices (1-1220).
* `prepare_parallel.sh`: Sets up the execution environment by creating isolated copies of the INTLAB library for each parallel worker.
* `run_matlab.sh`: A helper script called by `main_algo2.sh` and `main_algo3.sh` to execute a single MATLAB computation task.

### Core MATLAB Code (`Each_Process/`)

* `func_algo2.m`: Implements the core logic for **Algorithm 2**. It computes rigorous eigenvalue bounds for each subdomain `R_ij`.
* `func_algo3.m`: Implements the core logic for **Algorithm 3**, using domain monotonicity.
* `my_intlab_mode_config.m`: Configures the INTLAB environment for each specific parallel worker process.
* **`FEM_Functions/` directory:** Contains various helper functions for finite element calculations, including the Lehmann-Goerisch method (`calc_eigen_bounds_any_order.m`) and Lagrange element functions (`Lagrange_...`).
* **`mesh/` directory:** Contains helper functions responsible for generating the computational triangulations ($\mathcal{T}^{h}$) of the domains (`get_mesh_...`).
