# SpatialRegret

[![License: 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](http://creativecommons.org/licenses/by/4.0/)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2021a+-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![arXiv](https://img.shields.io/badge/arXiv-2511.14280-b31b1b.svg)](https://arxiv.org/abs/2511.14280)

**Distributed Controller Synthesis using Spatial Regret Optimization**

This repository contains the MATLAB code accompanying the paper ["A graph-informed regret metric for optimal distributed control"](https://arxiv.org/abs/2511.14280) (Martinelli et al., 2025). The code implements distributed controller synthesis for networked dynamical systems using the spatial regret framework, with both L1-optimal and SDP-based (H2/H∞) synthesis methods, supporting System Level Synthesis (SLS) and Youla parameterization approaches.

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Repository Structure](#repository-structure)
- [Main Scripts](#main-scripts)
- [Core Algorithms](#core-algorithms)
- [Examples](#examples)
- [Requirements](#requirements)
- [Citation](#citation)
- [Contributing](#contributing)
- [License](#license)

## Overview

**Spatial Regret** is a performance metric for distributed control that measures the worst-case performance degradation compared to an oracle controller with enhanced information-sharing capabilities. This framework enables the design of distributed controllers that:

- Respect communication constraints (sparsity patterns)
- Minimize performance loss relative to idealized architectures
- Scale to large networked systems
<!-- - Provide provable performance guarantees -->

The repository implements synthesis methods for:
- **H2/H∞/L1-optimal control** using System Level Synthesis (SLS)
- **H2/H∞ control** using Youla parameterization
- **Spatial regret minimization** for both approaches
- **Distributed optimization** using ADMM (Alternating Direction Method of Multipliers)

## Key Features

**Multiple Synthesis Methods**
- System Level Synthesis (SLS) with FIR approximations
- Youla parameterization with coprime factorizations
- Transfer function sampling on the unit circle

**Flexible Communication Topologies**
- Grid networks (arbitrary dimensions)
- Custom adjacency matrices
- Distance-based sparsity constraints (Input/Output delay constraints for the resulting controller)

**Performance Metrics**
- L1 norm (peak performance)
- H2 norm (average performance)
- H∞ norm (worst-case frequency response)
- Spatial regret norm (oracle comparison)

**Scalable Optimization**
- Centralized convex optimization
- Distributed ADMM with row/column decomposition
- Warm-start capabilities
- Multiple solver support (Gurobi, Mosek, SEDUMI)

## Installation

### Prerequisites

1. **MATLAB** (R2021a or later recommended)
   - Control System Toolbox
   - Optimization Toolbox

2. **YALMIP** - MATLAB optimization modeling toolbox
   ```bash
   # Download from: https://yalmip.github.io/
   # Add to MATLAB path
   ```

3. **Convex Optimization Solver** (at least one):
   - [Gurobi](https://www.gurobi.com/) (recommended for L1 problems)
   - [MOSEK](https://www.mosek.com/) (recommended for SDP problems)
   - SEDUMI (free, included with YALMIP)

### Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/DecodEPFL/SpatialRegret.git
   cd SpatialRegret
   ```

2. Add to MATLAB path:
   ```matlab
   addpath(genpath('./Functions_SpRegret'));
   ```

3. Verify installation:
   ```matlab
   % Run the example script
   EXAMPLE_FOR_USE
   ```

## Quick Start

### Synthesizing a Spatial Regret Optimal Controller

```matlab
% Define the LFT system structure (see EXAMPLE_FOR_USE.m for details)
% from the ``plant'' we aim to control
sys.A = plant.A;
sys.B2 = plant.B;
sys.C2 = plant.C;
sys.D22 = plant.D;
sys.plant = plant;
Adjacency_matrix = ...
sys.n = size(sys.A, 1);
sys.m = size(sys.B2, 2);
sys.p = size(sys.C2, 1);

% Define performance weights and disturbance structure
sys.C1 = [sqrtm(eye(sys.n)); zeros(sys.m, sys.n)];
sys.D12 = [zeros(sys.n, sys.m); sqrtm(eye(sys.m))];
sys.n_z = size(sys.C1, 1);
sys.n_w = n_agents;
sys.B1 = kron(eye(n_agents), [0; 1]);
sys.D11 = zeros(sys.n_z, sys.n_w);
sys.D21 = eye(sys.p, sys.n_w);

% Compute coprime factorization for Youla parameterization
[F, L] = calculus_F_and_L(sys, Adjacency_matrix); %Not required if plant is already stable
sys.F = F; sys.L = L; 
[P11, P12, P21] = coprime_factorization(sys);
sys.P11 = P11; sys.P12 = P12; sys.P21 = P21;

% Define communication delays and oracle structure
Graph = digraph(Adjacency_matrix ~= 0);
delays_matrix = distances(Graph);
oracle_delays = delays_matrix;
oracle_delays(:, end) = 0;  % All agents share with last agent

% Synthesize oracle controller
options = get_default_options('N_tf', 20, 'method', 'youla');
[K_oracle, Q_oracle, ~] = calculus_distributed(sys, 'hinf', oracle_delays, options);
lft_oracle = sys.P11 - sys.P12 * Q_oracle * sys.P21;

% Synthesize spatial regret controller
options_spreg = get_default_options('number_points', 2000, 'method', 'tf_sampled');
[K_spreg, ~, spreg_cost] = calculus_spatial_regret(sys, lft_oracle, ...
                                                    delays_matrix, options_spreg);

fprintf('Spatial regret cost: %.4f\n', spreg_cost);
```

### Results of Section IV.A (SDP formulation)

```matlab
% Run the 5-agent chain example with SDP methods
main_SDP_GRID  % Synthesizes H2, Hinf, Oracle, and Spatial Regret controllers
```

### Results of Section IV.B (L1 formulation)

```matlab
% Run the 16-agent grid example
main_L1_GRID  % Synthesizes L1, Oracle, and Spatial Regret controllers
```



## Repository Structure

```
SpatialRegret/
├── main_L1_GRID.m           # Main script for 16-agent grid (L1 synthesis)
├── main_SDP_GRID.m          # Main script for 5-agent chain (SDP synthesis)
├── EXAMPLE_FOR_USE.m        # Tutorial example script
├── Functions_SpRegret/      # Core algorithms and utilities
│   ├── calculus_spatial_regret.m        # Spatial regret synthesis (SDP)
│   ├── calculus_spatial_regret_L1.m     # Spatial regret synthesis (L1)
│   ├── calculus_distributed.m           # Distributed controller synthesis
│   ├── spregnorm.m                      # Spatial regret norm computation
│   ├── generate_plant_homogeneous.m     # Plant model generation
│   ├── coprime_factorization.m          # Youla parameterization setup
│   ├── sls_achievability_constraints.m  # SLS constraints
│   ├── plots_for_L1.m                   # Visualization for L1 results
│   ├── plots_for_SDP.m                  # Visualization for SDP results
│   └── ... (40+ utility functions)
├── figures/                 # Where generated plots are stored
├── results/                 # Where simulation results are stored
├── README.md               # This file
├── LICENSE                 # MIT License
└── CONTRIBUTING.md         # Contribution guidelines
```

## Main Scripts

### `main_L1_GRID.m`
Synthesizes L1-optimal controllers for a 4×4 grid network (16 agents).

**Key steps:**
1. Build 16-node grid topology
2. Synthesize Oracle controller (enhanced communication)
3. Synthesize Spatial Regret controller (L1, SLS)
4. Synthesize baseline L1 controller
5. Compare L1 norms and generate plots

**Output:** Controller norms, frequency plots, impulse response simulations

### `main_SDP_GRID.m`
Synthesizes H2/H∞ controllers for a 5-agent linear chain using Youla parameterization.

**Key steps:**
1. Build 5-node chain topology
2. Synthesize Oracle controller (H∞, enhanced information)
3. Synthesize Spatial Regret controller
4. Synthesize H2 and H∞ baseline controllers
5. Compare all three norms (H2, H∞, Spatial Regret)

**Output:** Performance comparison table, frequency responses, time-domain simulations

### `EXAMPLE_FOR_USE.m`
Tutorial script demonstrating the complete workflow with detailed comments.

## Core Algorithms

### Spatial Regret Synthesis

**`calculus_spatial_regret.m`** - SDP-based spatial regret minimization
```matlab
[K, Q, cost] = calculus_spatial_regret(sys, oracle_LFT, delays_matrix, options)
```
- Minimizes: $\sup_{\omega} \sigma_{\max}(T_{zw}(e^{j\omega}) - T_{zw}^{\text{oracle}}(e^{j\omega}))$
- Methods: `'sls'`, `'sampled_youla'`, `'tf_sampled'`

**`calculus_spatial_regret_L1.m`** - L1-based spatial regret minimization
```matlab
[K, Phis, cost] = calculus_spatial_regret_L1(sys, oracle_Phis, delays_matrix, options)
```
- Minimizes: $\|T_{zw} - T_{zw}^{\text{oracle}}\|_{\ell_1}$
- Supports centralized and distributed (ADMM) optimization

### Distributed Controller Synthesis

**`calculus_distributed.m`** - General distributed synthesis
```matlab
[K, Q_or_Phis, objective] = calculus_distributed(sys, problem_type, delays, options)
```
- `problem_type`: `'h2'`, `'hinf'`, `'l1'`
- `options.method`: `'sls'`, `'youla'`, `'sampled_youla'`

### Performance Evaluation

**`spregnorm.m`** - Compute spatial regret norm
```matlab
lambda = spregnorm(system_LFT, oracle_LFT, number_points)
```
Evaluates: $\inf \{\lambda : T_{zw}^*(e^{j\omega}) T_{zw}(e^{j\omega}) \preceq \lambda I + T_{zw}^{\text{oracle},*}(e^{j\omega}) T_{zw}^{\text{oracle}}(e^{j\omega}), \forall \omega\}$

## Examples

The repository includes several examples of increasing complexity:

1. **`EXAMPLE_FOR_USE.m`**: Basic 4-agent system with detailed explanations
2. **`main_SDP_GRID.m`**: 5-agent chain with H2/H∞ synthesis
3. **`main_L1_GRID.m`**: 16-agent grid with L1 synthesis and distributed optimization

### Running Examples

```matlab
% Example 1: Tutorial (recommended starting point)
EXAMPLE_FOR_USE

% Example 2: Small-scale SDP synthesis
main_SDP_GRID

% Example 3: Large-scale L1 synthesis
main_L1_GRID
```

### Expected Runtime
- `EXAMPLE_FOR_USE`: ~1-2 minutes
- `main_SDP_GRID`: ~5-10 minutes
- `main_L1_GRID`: ~1-5 minutes (depending on solver and hardware)

## Requirements

### MATLAB Toolboxes
- **Control System Toolbox** (required)
- **Optimization Toolbox** (required)
- **Symbolic Math Toolbox** (optional, for some utilities)

### External Dependencies
- **YALMIP** (required) - [Download](https://yalmip.github.io/)
- **Optimization Solver** (at least one):
  - Gurobi (recommended for large-scale L1 problems)
  - MOSEK (recommended for SDP problems)
  - SEDUMI (free alternative, slower)

### System Requirements
- MATLAB
- 8GB+ RAM recommended for 16-agent examples
- Multi-core CPU beneficial for distributed optimization

## Citation

If you use this code in your research, please cite:

```bibtex
@misc{martinelli2025graphinformedregretmetricoptimal,
      title={A graph-informed regret metric for optimal distributed control}, 
      author={Daniele Martinelli and Andrea Martin and Giancarlo Ferrari-Trecate and Luca Furieri},
      year={2025},
      eprint={2511.14280},
      archivePrefix={arXiv},
      primaryClass={eess.SY},
      url={https://arxiv.org/abs/2511.14280}
}
```

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on:
- Reporting bugs
- Suggesting enhancements
- Code style conventions
- Pull request process

Quick contribution checklist:
- ✅ Follow MATLAB coding standards
- ✅ Add documentation to new functions
- ✅ Test on small examples
- ✅ Update README if adding features

## License
This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by] 

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

## Acknowledgments

This work was developed at [EPFL DECODE Lab](https://www.epfl.ch/labs/decode/).

## Contact

For questions or issues:
- Open an issue on [GitHub](https://github.com/DecodEPFL/SpatialRegret/issues)
- Contact the maintainers via the repository
- Contact the maintainers via the institutional mail.

<!-- ## Troubleshooting

### Common Issues

**Problem:** "Undefined function or variable 'sdpvar'"
- **Solution:** Install YALMIP and add to MATLAB path

**Problem:** "No solver found"
- **Solution:** Install at least one optimization solver (Gurobi, MOSEK, or SEDUMI)

**Problem:** Slow L1 optimization
- **Solution:** Use Gurobi solver, reduce `N_tf`, or enable distributed optimization

**Problem:** Out of memory errors
- **Solution:** Reduce system size, decrease `N_tf`, or use a machine with more RAM -->

### Getting Help

1. Check the example scripts for usage patterns
2. Review function documentation (type `help function_name`)
3. Open an issue with a minimal reproducible example


---

**Last Updated:** November 2025
