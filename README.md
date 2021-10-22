# Codes of "Rigorous numerics for nonlinear heat equations in the complex plane of time"

This repository contains the MATLAB codes associated with the paper:
"Rigorous numerics for nonlinear heat equations in the complex plane of time"
by A. Takayasu, J.-P. Lessard, J. Jaquette, and H. Okamoto. ([arXiv:1910.12472 [math.DS]](https://arxiv.org/abs/1910.12472))

**Abstract**: In this paper, we introduce a method for computing rigorous local inclusions of solutions of Cauchy problems for nonlinear heat equations for complex time values. The proof is constructive and provides explicit bounds for the inclusion of the solution of the Cauchy problem, which is rewritten as a zero-finding problem on a certain Banach space. Using a solution map operator, we construct a simplified Newton operator and show that it has a unique fixed point. The fixed point together with its rigorous bounds provides the local inclusion of the solution of the Cauchy problem. The local inclusion technique is then applied iteratively to compute solutions over long time intervals. This technique is used to prove the existence of a branching singularity in the nonlinear heat equation. Finally, we introduce an approach based on the Lyapunov-Perron method to  calculate part of a center-stable manifold and prove that an open set of solutions of the Cauchy problem converge to zero, hence yielding the global existence of the solutions in the complex plane of time.

These codes require *MATLAB* with [*INTLAB* - INTerval LABoratory](http://www.ti3.tu-harburg.de/rump/intlab/) (MATLAB toolbox for interval arithmetic) version 11 and [*Chebfun* - numerical computing with functions](https://www.chebfun.org/) version 5.7.0.

All computations are carried out on Windows 10, Intel(R) Core(TM) i7-6700K CPU @ 4.00GHz, and MATLAB 2019a with INTLAB - INTerval LABoratory version 11 and Chebfun - numerical computing with functions version 5.7.0.

_Note that this code was originally published at [this link](https://www.risk.tsukuba.ac.jp/~takitoshi/codes/RNcnheq.zip), and has been moved to this repository._

---

Before running the codes one should addpath of Chebfun toolbox and INTLAB, e.g., `addpath('chebfun-master/') ` and `addpath('Intlab/')`

Firstly, one add the following paths:

```
addpath('verify_solution/')
addpath('verify_defect/')
addpath('variational_problem/')
```

### Proof of Theorem 1.1

```
% script_proof_upper_bound_blowup % uncomment when you execute the proof of Theorem 1.1
% script_proof_lower_bound_blowup % uncomment when you execute the proof of Theorem 1.1
script_plot_fig3
script_plot_fig4
```

### Proof of Theorem 1.2

```
% script_proof_of_GE_60 % uncomment when you execute the proof of Global existence (GE) in the case of \theta = 60
% script_proof_of_GE_45 % uncomment when you execute the proof of Global existence (GE) in the case of \theta = 45
% script_proof_of_GE_30 % uncomment when you execute the proof of Global existence (GE) in the case of \theta = 30
% script_proof_of_GE_15 % uncomment when you execute the proof of Global existence (GE) in the case of \theta = 15
script_plot_fig5
```
---

Copyright (C) 2019  A. Takayasu, J.-P. Lessard, J. Jaquette, and H. Okamoto.
