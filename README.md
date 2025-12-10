# Computer-Assisted Proof for Dirichlet Eigenvalue Simplicity

This project provides the source code and computational framework for the computer-assisted proof presented in the paper:

> **[Paper B] The Second Dirichlet Eigenvalue is Simple on Every Non-equilateral Triangle, Part II: Nearly Equilateral Triangles** \> *(to appear in Numerische Mathematik)*

The primary goal is to rigorously validate the simplicity of the second Dirichlet eigenvalue for nearly equilateral triangles, offering a complete solution to a conjecture posed by R. Laugesen and B. Siudeja.

The computer-assisted proof for nearly degenerate triangles is presented in the following paper:

> **[Paper A] The second Dirichlet eigenvalue is simple on every non-equilateral triangle, Part I: Nearly degenerate triangles** \> *Journal of Differential Equations 447, 113629*
> [https://doi.org/10.1016/j.jde.2025.113629](https://doi.org/10.1016/j.jde.2025.113629)

## Background

Determining the eigenvalue multiplicity of the Laplace operator is challenging, especially when eigenvalues are nearly degenerate—as for the second and third Dirichlet eigenvalues on equilateral triangles where $\lambda_{2}=\lambda_{3}$. Standard numerical methods struggle to separate these tightly clustered eigenvalues with mathematical rigor.

We split the parameter space of triangles $\Omega$ into three regions:

  * **Region** $\Omega_{\text{up}}$ (nearly equilateral triangles): Using **Algorithm 1**, we prove the difference quotient satisfies $D\lambda_{2} < D\lambda_{3}$, implying separation.
  * **Region** $\Omega_{\text{down}}$: Using **Algorithm 2**, we compute high-precision bounds directly to show $\overline{\lambda}_2 < \underline{\lambda}_3$.
  * **Region** $\Omega_{\text{rest}}$ (nearly degenerate triangles): This case is handled in **[Paper A]**.

For the definition of each region, see [Paper B].

## Core Libraries & Dependencies

This project relies on specialized libraries for verified numerical computation:

1.  **INTLAB**: The fundamental toolbox for rigorous interval arithmetic in MATLAB.
      * **Source:** [http://www.tuhh.de/ti3/intlab/](http://www.tuhh.de/ti3/intlab/)
2.  Revised version of **VFEM2D**: Used for rigorous finite element matrix assembly and high-precision eigenvalue bounds (Lehmann–Goerisch method).
      * **Source:** [https://github.com/xfliu/VFEM2D](https://github.com/xfliu/VFEM2D)
3.  **veigs**: Used for solving generalized matrix eigenvalue problems with rigorous error bounds with the information of indices.
      * **Source:** [https://github.com/yuuka-math/veigs](https://github.com/yuuka-math/veigs)

## Project Structure

```text
.
├── ProofRunner.m                 # Orchestrator class (Main entry point)
├── verification_step_1.m         # Core logic for Algorithm 1 (Difference Quotients)
├── verification_step_2.m         # Core logic for Algorithm 2 (Domain Monotonicity)
├── inputs/
│   └── cell_def.csv              # Input grid definition for Algorithm 2
├── FEM_Functions/                # Helper functions for FEM operations
├── mesh/                         # Mesh generation routines
├── mode_swith_interface/         # Interface for switching calculation modes
├── VFEM2D/                       # Original VFEM2D library
├── VFEM2D_revised/               # Revised VFEM2D library for this proof
├── results/                      # Output directory
│   ├── results_algo1.csv         # Output from Algo 1
│   └── results_algo2.csv               # Output from Algo 2
└── my_intlab_config.m            # INTLAB configuration
````

## Installation & Setup

Before running the code, ensure you have **MATLAB** (R2020b or later).

### 1\. Clone the Repository

```bash
git clone [https://github.com/ryendo/DirichletSimplicityClustered](https://github.com/ryendo/DirichletSimplicityClustered)
cd DirichletSimplicityClustered
```

### 2\. Configure INTLAB

This project requires the **INTLAB** (Interval Laboratory) toolbox. You must configure the path to your local INTLAB installation.

Open `my_intlab_config.m` and edit the `addpath` line:

```matlab
% Open my_intlab_config.m
addpath('/path/to/your/INTLAB_directory'); 
% e.g., addpath('/Applications/Intlab_V12');
```

-----

## Usage: The `ProofRunner` Class

The `ProofRunner` class is the central controller for this project. It handles setup, execution, and provides utility methods for debugging or verifying specific cases.

### 1\. Initialization

Load the configuration and prepare the environment within MATLAB.

```matlab
s = ProofRunner;
s.setupAll();  % Initializes INTLAB
```

### 2\. Algorithm 1: Difference Quotients ($\Omega_{\text{up}}$)

This algorithm verifies that the difference quotients of eigenvalues separate as we move away from the equilateral shape.

#### Run Full Verification

Sweeps the entire angular range $\delta \in [0, \pi/3]$. Results are saved to `results/results_algo1.csv`.

```matlab
s.runAlgo1All();
```

**Output Sample:**

```text
--- Starting Algorithm 1 Loop ---
[=================...............] 532/1000 (53.2%) ETA 00h:12m:30s
...
```

#### Check a Specific Interval

```matlab
% Verify the interval [0.1, 0.2] by dividing it into 5 sub-bins.
s.runAlgo1Interval([intval('0.1'), intval('0.2')], 5);
```

#### Verify Results (`summarizeAlgo1CSV`)

Analyzes the generated CSV to confirm the mathematical proof condition ($\sup \mu_1 < \inf \mu_2$).

```matlab
s.summarizeAlgo1CSV(); 
```

**Output Sample:**

```text
--- Algo 1 Summary ---
File: results/results_algo1.csv
Total: 1000, Passed: 1000, Failed: 0
>>> PROOF SUCCESSFUL (Algo 1)
```

### 3\. Algorithm 2: Global Verification ($\Omega_{\text{down}}$)

This algorithm verifies the eigenvalue gap on the remaining domain using domain monotonicity. It works in a batch mode by processing a grid definition file.

**Prerequisite:** Ensure `inputs/cell_def.csv` exists. This file defines the grid cells ($x, \theta$) and the FEM parameters required for verification.

#### Execution (`runAlgo2All`)

This command reads `inputs/cell_def.csv`, computes rigorous bounds for each cell, and appends the results to `results/results.csv`.

*Note: This process initializes a new results file and processes the entire batch defined in the CSV.*

```matlab
% Run using default input path 'inputs/cell_def.csv'
s.runAlgo2All();
```

**Output Sample:**

```text
--- Algorithm 2: Running Batch Verification ---
Input: inputs/cell_def.csv
Output: results/results.csv
Initialized results.csv with headers.
Processing cell i=1 ...
Processing cell i=2 ...
...
Algorithm 2 Batch Completed. Duration: 01h:20m:15s
```

#### Verify Results (`summarizeAlgo2CSV`)

Reads the output file (`results/results.csv`) and verifies the eigenvalue gap condition ($\sup \lambda_2 < \inf \lambda_3$) for every processed cell.

```matlab
s.summarizeAlgo2CSV();
```

**Output Sample:**

```text
--- Verifying Algorithm 2 Results: results/results.csv ---
>>> PROOF SUCCESSFUL (Algo 2): All 12514 cells verified.
```

### 4\. Pointwise & Box Verification Tools

These methods calculate rigorous eigenvalue bounds ($\lambda_2, \lambda_3$) for specific triangle shapes defined by the top vertex coordinate $(a, b)$ or box regions.

#### Single Point Check (`boundsAtPoint`)

Computes eigenvalues for a specific triangle vertex $(a, b)$ assuming base vertices at $(0,0)$ and $(1,0)$.

```matlab
% Check triangle with top vertex at (0.5, 0.8)
[lam2, lam3] = s.boundsAtPoint(intval('0.5'), intval('0.8'));
```

**Output Sample:**

```text
Point Bounds:
 lam2 in [42.1532, 42.1545]
 lam3 in [48.9012, 48.9025]
```

#### Box Check (`boundsOnBox`)

Verifies the gap condition over a rectangular region in the $(x, \theta)$ parameter space, corresponding to **Algorithm 4** in the paper.

The domain is defined by a triangle with vertices $A(0,0)$, $B(1,0)$, and $C(x,y)$. We use the coordinate system:

$$
x = x, \quad \theta = \arctan\left(\frac{y}{x}\right) \implies y = x \tan(\theta)
$$

Based on the geometric construction, for a cell defined by $x \in [x_{\min}, x_{\max}]$ and $\theta \in [\theta_{\min}, \theta_{\max}]$, the domain inclusion holds:

$$
T^{(x_{\min}, \theta_{\min})} \subset T^{(x, \theta)} \subset T^{(x_{\max}, \theta_{\max})}
$$

Due to the monotonicity of Dirichlet eigenvalues with respect to domain inclusion, we obtain:

  * **$\sup(\lambda_2)$** is computed using the "smallest" geometry defined by $(x_{\min}, \theta_{\min})$.
  * **$\inf(\lambda_3)$** is computed using the "largest" geometry defined by $(x_{\max}, \theta_{\max})$.

<!-- end list -->

```matlab
% Check box x=[0.5, 0.51], theta=[0.9, 0.91]
[up2, lo3] = s.boundsOnBox(intval('0.5'), intval('0.51'), ...
                           intval('0.9'), intval('0.91'));
```
