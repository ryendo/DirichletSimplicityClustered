# Computer-Assisted Proof for Dirichlet Eigenvalue Simplicity

This project provides the source code and computational framework for the computer-assisted proof presented in the paper “Rigorous estimation for the difference quotients of multiple eigenvalues”. The primary goal is to rigorously validate the simplicity of the second Dirichlet eigenvalue for nearly equilateral triangles, offering a partial solution to a conjecture posed by R. Laugesen and B. Siudeja (discussed as Conjecture 6.47 in A. Henrot’s *Shape Optimization and Spectral Theory*).

## Background

Determining the eigenvalue multiplicity of the Laplace operator is challenging, especially when eigenvalues are nearly degenerate—as for the second and third Dirichlet eigenvalues on equilateral triangles where $\lambda_{2}=\lambda_{3}$. Standard numerical methods struggle to separate these tightly clustered eigenvalues with mathematical rigor.

This work analyzes the **difference quotient of eigenvalues** instead of classical shape derivatives:

$$
D_{t}\lambda_{i} := \frac{\lambda_{i}(\Omega_{t})-\lambda_{i}(\Omega_{0})}{t}.
$$

By applying guaranteed computation techniques based on the Finite Element Method (FEM) and interval arithmetic, we obtain rigorous bounds on this quotient.

We split the parameter space of triangles $\Omega$ into two regions:

* Nearly equilateral $\Omega_{\text{up}}$: prove $D\lambda_{2}(p_{0},p) < D\lambda_{3}(p_{0},p)$ for any perturbation $p$ from the equilateral vertex $p_{0}$. Since $\lambda_{2}^{p_{0}}=\lambda_{3}^{p_{0}}$, this implies $\lambda_{2}^{p}<\lambda_{3}^{p}$.
* The complement $\Omega_{\text{down}}$: compute high-precision bounds directly to show a definitive gap, $\overline{\lambda_2}^{p}<\underline{\lambda_3}^{p}$.

## Project Structure

```
.
├── Each_Process/                 # Core MATLAB functions, libraries, logs, and temporary results
│   ├── functions/                # Helper MATLAB functions for computation
│   ├── HighOrderFEM_CGYOU_2016/  # FEM library
│   ├── Intlab_Group/             # INTLAB library for guaranteed computation
│   └── log/                      # Log files from parallel execution
├── prep/                         # Scripts for preparing the computation
├── results/                      # Final stored results from computations
├── ProofRunner.m                 # Orchestrator class (progress, ETA, selective runs)
├── main_algo1.m                  # Main MATLAB script for Algorithm 1 (unchanged)
├── main_algo2.sh                 # Main execution script for Algorithm 2 (unchanged)
├── my_intlab_config_alone.m      # Standalone INTLAB configuration
├── prep.sh                       # Master script to run all preparation steps
└── README.md                     # This file
```

## Setup and Prerequisites

### Prerequisites

1. **MATLAB:** A working installation of MATLAB.
2. **INTLAB Library:** Required for interval arithmetic and verified numerics.
3. **Unix-like Environment:** `bash`/`zsh` and common CLI tools to run shell scripts.

### Initial Configuration

1. **Clone the Repository**

   ```bash
   git clone https://github.com/ryendo/DirichletSimplicityClustered
   cd DirichletSimplicityClustered
   ```

2. **Enable `matlab` on the Command Line**

   ```bash
   # Add MATLAB bin to your PATH (adjust the path to your installation)
   export PATH="/path/to/your/matlab/bin:$PATH"
   ```

3. **Place INTLAB**
   Put `Intlab_V12` under:

   ```
   Each_Process/Intlab_Group/
   ```

## Quick Start with `ProofRunner` (Recommended)

```matlab
% QUICK START
%   s = ProofRunner;
%   s.setupAll();                      % INTLAB + prep.sh
%   s.runAlgo1All();                   % full Ω_up sweep with ETA/progress
%   s.summarizeAlgo1CSV();             % check sup(mu1) < inf(mu2)
%
%   s.runAlgo1Direction(I_pi/12);        % single direction near equilateral
%   s.runAlgo1Interval([I_pi/20,I_pi/20+intval('1e-5')]); % small interval
%
%   s.boundsAtPoint(intval('1')/2,sqrt(intval('3'))/2-intval('1e-3'));       % eigenvalue bounds at (s,t)
%   s.boundsOnBox(intval('1')/2,intval('1')/2,sqrt(intval('3'))/2-2*intval('1e-3'),sqrt(intval('3'))/2-intval('1e-3'));     %[s_inf,s_sup,t_inf,t_sup]
%
%   s.prepareAlgo2();                  % runs ./prep.sh
%   s.runAlgo2('detached',true);       % launch main_algo2.sh (Unix)
%   s.monitorAlgo2();                  % live progress from CSV
```

`ProofRunner` is a thin orchestration layer; it **does not modify** `main_algo1.m` or `main_algo2.sh`.

### Algorithm 1 (Ω_up): Difference-Quotient Analysis

Run the full sweep $\delta \in [0,\pi/3]$ with progress and ETA:

```matlab
>> s = ProofRunner;
>> s.setupAll();
>> s.runAlgo1All();
```

Single direction:

```matlab
>> s.runAlgo1Direction(I_pi/12);
```

Direction interval with a chosen number of bins:

```matlab
>> s.runAlgo1Interval(hull(I_pi/20, I_pi/15), 25);
```

Summarize results from a CSV:

```matlab
>> s.summarizeAlgo1CSV();                           % latest file
>> s.summarizeAlgo1CSV('results/quotients_....csv');% specific file
```

**Output CSV header**

```
idx,inf_mu1,sup_mu1,inf_mu2,sup_mu2
```

The proof condition per row is `sup_mu1 < inf_mu2`.

### Algorithm 2 (Ω_down): Parallel Eigenvalue Bounds

Prepare tasks:

```matlab
>> s.prepareAlgo2();     % runs ./prep.sh
```

Run and monitor:

```matlab
>> s.runAlgo2();                 % foreground
>> s.runAlgo2('detached', true); % background
>> s.monitorAlgo2();             % reads prep/algo2_list_j.csv and shows ETA
```

### Partial Evidence Near Equilateral

Bounds at a specific vertex $(s,t)$ offset from $(\tfrac12,\tfrac{\sqrt{3}}{2})$:

```matlab
>> intval('1')/2,sqrt(intval('3'))/2-intval('1e-3');
% prints: λ2 ≤ up2,  λ3 ≥ lo3,  gap ≥ lo3 - up2
```

Bounds on a box $[a_{\inf},a_{\sup}] \times [t_{\inf},t_{\sup}]$:

```matlab
>> s.boundsOnBox(intval('1')/2,intval('1')/2,sqrt(intval('3'))/2-2*intval('1e-3'),sqrt(intval('3'))/2-intval('1e-3'));
```

---

## Running the Original Scripts (Baseline)

You may also use the legacy entry points.

### Algorithm 1

```matlab
>> my_intlab_config_alone
>> main_algo1
```

### Algorithm 2

```bash
./prep.sh
bash ./main_algo2.sh
```

Monitor with:

* `ps aux | grep -i matlab`
* `tail -f Each_Process/log/process_no1.log`
* Results appear in `results/` (separate CSVs per task).

---

## Verifying the Results

### Algorithm 1 (CSV example)

```csv
idx,inf_mu1,sup_mu1,inf_mu2,sup_mu2
1,98.3723341163222,109.514164604184,173.118847437632,186.288746858983
2,98.3656467084821,109.520542044462,173.110819954873,186.296463918525
3,98.3588074629973,109.526760266504,173.10264018102,186.304022212805
...
1000,26.3925588844467,39.7994719938128,100.131634868792,117.580812296542
```

Use:

```matlab
>> s.summarizeAlgo1CSV('results/quotients_...csv');
```

### Algorithm 2

The main script produces one CSV per task (e.g., `results/algo2_j_1.csv`). A region is verified when `upper(λ₂) < lower(λ₃)` for that task. The union of successful files constitutes the proof over (\Omega_{\text{down}}). Progress and ETA can be displayed via:

```matlab
>> s.monitorAlgo2();
```

---

## Code and Script Roles

### Root Directory

* `ProofRunner.m`: Orchestration class (progress, ETA, selective runs). Calls `prep.sh` and `main_algo2.sh` without modifying them.
* `main_algo1.m`: Algorithm 1 (unchanged).
* `main_algo2.sh`: Algorithm 2 (unchanged).
* `my_intlab_config_alone.m`: INTLAB init for single MATLAB sessions.

### `prep/`

* `make_indices_algo2.m`: Generates `prep/algo2_list_j.csv`.
* `prepare_parallel.sh`: Sets up per-worker environments (if used).
* `run_matlab.sh`: Helper invoked by `main_algo2.sh`.

### `Each_Process/`

* `func_algo2.m`: Core of Algorithm 2 (verified bounds per subdomain).
* `my_intlab_mode_config.m`: Per-worker INTLAB config.
* `functions/`: FEM helpers (e.g., `calc_eigen_bounds_any_order.m`, `Lagrange_...`, `get_mesh_...`).
