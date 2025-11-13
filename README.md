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

* Nearly equilateral $\Omega_{\text{up}}$: prove $D\lambda_{2}(p_{0},p) < D\lambda_{3}(p_{0},p)$ for any perturbation $p$ from the equilateral vertex $p_{0}$. Since $\lambda_{2}^{p_{0}}=\lambda_{3}^{p_{0}}$, this implies $\lambda_{2}^{p}<\lambda_{3}^{p} $.
* The complement $\Omega_{\text{down}}$: compute high-precision bounds directly to show a definitive gap,
  $\overline{\lambda}_2^{p}<\underline{\lambda}_3^{p}$.

## Method Overview

**Algorithm 1 ($\Omega_{\text{up}}$) — Difference-Quotient near equilateral.**  
This algorithm rigorously estimates the difference quotients $D_t\lambda_2$ and $D_t\lambda_3$ for triangles obtained by perturbing the equilateral triangle along a direction $e$ (or a small direction interval). By sweeping directions $\delta \in [0,\pi/3]$ and enclosing the generalized eigenvalues $\mu_1,\mu_2$ of the quotient problem with interval arithmetic, the proof is established when $\sup(\mu_1) < \inf(\mu_2)$ holds for all sectors.

**Algorithm 2 ($\Omega_{\text{down}}$) — Verified eigenvalue gap away from equilateral.**  
This algorithm tiles the remaining parameter space into rectangles $R_{ij}$ and, for each subregion, computes certified upper/lower bounds $\overline{\lambda}_2$ and $\underline{\lambda}_3$ using high-order FEM and interval arithmetic (Lehmann–Goerisch-based routines). A subregion is certified when $\overline{\lambda}_2 < \underline{\lambda}_3$; collecting all certified subregions yields the proof on $\Omega_{\text{down}}$.


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

2. **MATLAB on the command line (only if you run shell scripts):**
   If you will run `prep.sh` / `main_algo2.sh`, make sure the shell can invoke MATLAB.
    Put the PATH tweak inside `prep.sh`:

   ```bash
   # --- at the very top of prep.sh (optional, if 'matlab' is not already on PATH) ---
   export PATH="/path/to/your/matlab/bin:$PATH"   # e.g., /Applications/MATLAB_R2024a.app/bin
   ```

3. **Place INTLAB:**
   Open `my_intlab_config_alone.m` and replace this line:
   ```matlab
   addpath("/path/to/your/INTLAB")
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

This algorithm estimates the difference-quotient for the 2nd and 3rd Dirichlet eigenvalues near the equilateral triangle. It perturbs the domain along a chosen direction $e$ (or a small direction interval), assembles a verified $2\times2$ generalized eigenproblem for each angular sector $\delta$, and records rigorous enclosures $(\mu_1,\mu_2)$ to CSV. The desired ordering $D\lambda_2 < D\lambda_3$ is concluded when $\sup(\mu_1) < \inf(\mu_2)$ for all processed sectors.


**Run the full sweep $\delta \in [0,\pi/3]$ with progress and ETA:**

```matlab
>> s = ProofRunner;
>> s.setupAll();
>> s.runAlgo1All();
```

**Single direction:**

```matlab
>> s.runAlgo1Direction(I_pi/12);
```
(Result)
```matlab
--- Running Algorithm 1 on 1 bin(s) ---
δ=[0.261799,0.261799]: 
 Dtλ2 ∈ [95.2543287913,102.967388819] 
 Dtλ3 ∈ [170.315416349,179.427292308],
 gap ≥ 67.3480275294
```

**Direction interval with a chosen number of bins:**

```matlab
>> s.runAlgo1Interval([I_pi/20, I_pi/20+intval('1e-4')], 5);
```
(Result)
```matlab--- Algorithm 1: rigorous difference-quotient bounds (Ω_up) ---
What is computed:
  - Verified interval bounds for the difference quotients of λ₂ and λ₃
    over a direction interval (nearly equilateral triangle).
Parameters:
  mesh_N=32, ord=5, ep=1e-05, omega_N=1000
Planned bins (K): 5 
Output CSV: results/quotients_2025-11-12_15-47-42.csv
Planned steps: 5 | warm-up per-step ≈ 0.34s | ETA ≈ 00h:00m:02s
--- Running Algorithm 1 on 5 bin(s) ---
δ=[0.15708,0.1571]: 
 Dtλ2 ∈ [98.6215102961,105.772579022] 
 Dtλ3 ∈ [173.747700975,182.167431973],
 gap ≥ 67.9751219528
[======..........................] 1/5 (20.0%) ETA 00h:00m:01sδ=[0.1571,0.15712]: 
 Dtλ2 ∈ [98.6210089801,105.772192796] 
 Dtλ3 ∈ [173.747187234,182.167058161],
 gap ≥ 67.9749944389
[=============...................] 2/5 (40.0%) ETA 00h:00m:01sδ=[0.15712,0.15714]: 
 Dtλ2 ∈ [98.6205076095,105.771806512] 
 Dtλ3 ∈ [173.746673439,182.166684292],
 gap ≥ 67.9748669274
[===================.............] 3/5 (60.0%) ETA 00h:00m:00sδ=[0.15714,0.15716]: 
 Dtλ2 ∈ [98.6200061841,105.771420171] 
 Dtλ3 ∈ [173.74615959,182.166310366],
 gap ≥ 67.9747394186
[==========================......] 4/5 (80.0%) ETA 00h:00m:00sδ=[0.15716,0.15718]: 
 Dtλ2 ∈ [98.6195047041,105.771033773] 
 Dtλ3 ∈ [173.745645685,182.165936382],
 gap ≥ 67.9746119125
[================================] 5/5 (100.0%) ETA 00h:00m:00s
Done. Results → results/quotients_2025-11-12_15-47-42.csv
```

**Summarize results from a CSV:**

```matlab
>> s.summarizeAlgo1CSV();                           % latest file
>> s.summarizeAlgo1CSV('results/quotients_....csv');% specific file
```
(Result)
```matlab
>> s.summarizeAlgo1CSV();
File: results/quotients_2025-11-12_15-47-42.csv
Bins: 5 | OK: 5 | NG: 0
>>> PROOF SUCCESSFUL: sup(mu1) < inf(mu2) holds for all bins.
```

**Output CSV header**

```
idx,inf_mu1,sup_mu1,inf_mu2,sup_mu2
```

The proof condition per row is `sup_mu1 < inf_mu2`.

### Algorithm 2 (Ω_down): Parallel Eigenvalue Bounds

This algorithm evaluates certified eigenvalue bounds for triangles away from the equilateral case. The parameter space is partitioned into subregions $R_{ij}$; for each, the code bounds $\lambda_2$ from above and $\lambda_3$ from below with high-order FEM and interval arithmetic. A subregion is certified when $\overline{\lambda}_2 < \underline{\lambda}_3$. The provided shell scripts launch and monitor these jobs in parallel and can safely resume after interruptions.

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

### Pointwise Eigenvalue Bounds

**Bounds at a specific vertex ((s,t)):**

```matlab
>> s.boundsAtPoint(intval('1')/2,sqrt(intval('3'))/2-intval('1e-3'));
```
(Result)
```matlab
Point (a=[0.5,0.5], b=[0.865026,0.865026]): 
 λ2 ∈ [122.92567091,122.926307183] 
 λ3 ∈ [123.00155729,123.002192357],
 gap ≥ 0.0752501072466
```

**Bounds on a box $[s_{\inf},s_{\sup}] \times [t_{\inf},t_{\sup}]$:**

```matlab
>> s.boundsOnBox(intval('1')/2,intval('1')/2,sqrt(intval('3'))/2-intval('1e-3')-intval('1e-10'),sqrt(intval('3'))/2-intval('1e-3'));
```
(Result)
```matlab
Box [0.5,0.5]×[0.865026,0.865026]:  sup λ2 ≤ 122.926307244,  inf λ3 ≥ 123.00155729,  min gap ≥ 0.0752500465989
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

* `ProofRunner.m`: Orchestration class (progress, ETA, selective runs). Calls `prep.sh` and `main_algo2.sh`.
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
