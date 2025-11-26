Here is the updated **README.md**. The sections regarding **Project Structure**, **Installation**, and **Algorithm 2 Usage** have been modified to reflect the switch from the shell script orchestration to the direct MATLAB function call.

-----

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
  * **Region** $\Omega_{\text{rest}}$ (nearly degenerate triangles): This case is handeled in **[Paper A]**.

For the definition of each region, see [Papaer B].

## Core Libraries & Dependencies

This project relies on specialized libraries for verified numerical computation:

1.  **INTLAB**: The fundamental toolbox for rigorous interval arithmetic in MATLAB.

      * **Source:** [http://www.tuhh.de/ti3/intlab/](http://www.tuhh.de/ti3/intlab/)

2.  Revised version of **VFEM2D**: Used for rigorous finite element matrix assembly and high-precision eigenvalue bounds (Lehmann–Goerisch method).

      * **Source:** [https://github.com/xfliu/VFEM2D](https://github.com/xfliu/VFEM2D)

3.  **veigs**: Used for solving generalized matrix eigenvalue problems with rigorous error bounds with the information of indices.

      * **Source:** [https://github.com/yuuka-math/veigs](https://github.com/yuuka-math/veigs)

## Project Structure

```
.
├── ProofRunner.m                 # Orchestrator class (Main entry point)
├── verification_step_1.m         # Core logic for Algorithm 1 (Difference Quotients)
├── verification_step_2.m         # Core logic for Algorithm 2 (Domain Monotonicity)
├── routines/                     # Core MATLAB functions and libraries
├── veigs/                        # veig library
├── VFEM2D_revised/               # Revised version of VFEM2D
├── results/                      # Output directory
└── my_intlab_config_alone.m      # INTLAB configuration
```

## Installation & Setup

Before running the code, ensure you have **MATLAB** (R2020b or later).

### 1\. Clone the Repository

```bash
git clone https://github.com/ryendo/DirichletSimplicityClustered
cd DirichletSimplicityClustered
```

### 2\. Configure INTLAB

This project requires the **INTLAB** (Interval Laboratory) toolbox. You must configure the path to your local INTLAB installation.

Open `my_intlab_config_alone.m` and edit the `addpath` line:

```matlab
% Open my_intlab_config_alone.m
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

Sweeps the entire angular range $\delta \in [0, \pi/3]$.

```matlab
s.runAlgo1All();
```

**Output Sample:**

```text
--- Running Algorithm 1 on 1000 bin(s) ---
[=================...............] 532/1000 (53.2%) ETA 00h:12m:30s
...
Done. Results -> results/quotients_2025-11-25_12-00-00.csv
```

#### Check a Single Direction (`runAlgo1Direction`)

Useful for debugging a specific perturbation direction $\delta$ (angle in radians).

```matlab
% Check the direction delta = pi/12
s.runAlgo1Direction(intval('pi')/12);
```

**Output Sample:**

```text
--- Algorithm 1: Rigorous difference-quotients ---
Parameters: mesh_N=32, ord=5, ep=1e-5, omega_N=1
Bins: 1
...
δ=[0.2617,0.2618]: gap >= 65.39
Done. Results -> results/quotients_....csv
```

#### Check a Specific Interval (`runAlgo1Interval`)

Verifies a custom sub-interval of angles.

```matlab
% Check interval [0.1, 0.2] divided into 5 bins
s.runAlgo1Interval([intval('0.1'), intval('0.2')], 5);
```

#### Verify Results (`summarizeAlgo1CSV`)

Analyzes the generated CSV to confirm the mathematical proof condition ($\sup \mu_1 < \inf \mu_2$).

```matlab
s.summarizeAlgo1CSV(); 
% Or specify a file: s.summarizeAlgo1CSV('results/quotients_xxx.csv');
```

**Output Sample:**

```text
Bins: 1000 | OK: 1000 | NG: 0
>>> PROOF SUCCESSFUL: sup(mu1) < inf(mu2) for all checked bins.
```

### 3\. Algorithm 2: Global Verification ($\Omega_{\text{down}}$)

This algorithm verifies the eigenvalue gap on the remaining domain using domain monotonicity. It iterates over a grid of angular indices ($j$) and radial indices ($i$).

#### Execution

You can run the full verification or process specific batches of indices.

```matlab
% Run the full verification (all j indices from 1 to N)
s.runAlgo2();

% Run specific indices (e.g., j from 1 to 10) for testing or partial verification
s.runAlgo2(1:10);
```

**Output Sample:**

```text
--- Algorithm 2: Domain Monotonicity (Direct Call) ---
Grid: M=40 (x), N=200 (theta)
Processing 200 j-indices. Output prefix: results/step2_bounds
...
Processing j = 1 / 200 ...
Processing j = 2 / 200 ...
```

#### Verify Results (`summarizeAlgo2Results`)

Reads all output CSV files (`step2_bounds_*.csv`) generated by the algorithm and verifies the eigenvalue gap condition ($\sup \lambda_2 < \inf \lambda_3$).

```matlab
s.summarizeAlgo2Results();
```

**Output Sample:**

```text
--- Verifying Algorithm 2 Results in results ---
Total subregions checked: 8000
>>> PROOF SUCCESSFUL (Algo 2): sup(lam2) < inf(lam3) verified for all loaded data.
```

### 4\. Pointwise & Box Verification Tools

These methods calculate rigorous eigenvalue bounds ($\lambda_2, \lambda_3$) for specific triangle shapes defined by the top vertex coordinate $(a, b)$.

#### Single Point Check (`boundsAtPoint`)

Computes eigenvalues for a specific triangle.

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

Verifies the gap condition over a rectangular region in the $(r, h)$ parameter space.

The domain is defined by a triangle with vertices $A(0,0)$, $B(1,0)$, and $C(x,y)$. We use the coordinate transformation:

$$
r = \sqrt{x^2+y^2}, \quad h = y
$$

Due to the monotonicity of Dirichlet eigenvalues with respect to domain inclusion, the eigenvalues are monotonically decreasing in both $r$ and $h$.
Therefore, the uniform upper bound $\sup(\lambda_2)$ for $\lambda_2$ and the uniform lower bound $\inf(\lambda_3)$ for $\lambda_3$ in the region $R:=[r_{\min},r_{\max}]\times[h_{\min},h_{\max}]$ are obtained in the following way:

  * **$\sup(\lambda_2)$** is computed using the "smallest" geometry $(r_{\min}, h_{\min})$.
  * **$\inf(\lambda_3)$** is computed using the "largest" geometry $(r_{\max}, h_{\max})$.

<!-- end list -->

```matlab
% Check box r=[1.0, 1.01], h=[0.80, 0.81]
[up2, lo3] = s.boundsOnBox(intval('1.0'), intval('1.01'), ...
                           intval('0.80'), intval('0.81'));
```

**Output Sample:**

```text
Box Bounds [(r,h)]: sup(lam2) <= 128.47474683161985, inf(lam3) >= 136.91708321477296
```