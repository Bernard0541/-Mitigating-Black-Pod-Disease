# -Mitigating-Black-Pod-Disease
Mathematical modeling and analysis of Black Pod Disease (BPD) transmission dynamics in cocoa plants using compartmental modeling, bifurcation theory, optimal control, sensitivity analysis, uncertainty quantification, and cost-effectiveness analysis.

# Mitigating Black Pod Disease in Cocoa Plants: An Optimal Control and Cost-Effectiveness Analysis

## Overview

Black Pod Disease (BPD), caused primarily by *Phytophthora palmivora* and *Phytophthora megakarya*, remains one of the major threats to global cocoa production. This repository contains mathematical and computational implementations for studying BPD transmission dynamics and evaluating intervention strategies using optimal control theory and cost-effectiveness analysis.

This work develops a nonlinear compartmental model that incorporates disease transmission through primary and secondary infections, investigates disease persistence through bifurcation analysis, and evaluates intervention strategies under resource limitations.

---

## Objectives

The project aims to:

- Develop a mathematical model for Black Pod Disease transmission
- Compute the basic reproduction number
- Study local and global stability properties
- Investigate backward bifurcation behavior
- Perform uncertainty and sensitivity analyses
- Evaluate optimal intervention strategies
- Perform cost-effectiveness analysis
- Identify efficient and sustainable disease management strategies

---

## Model Compartments

The model consists of:

- Sc : Cherelles pods
- Sp : Young and mature pods
- Sr : Ripe cocoa pods
- E : Exposed cocoa pods
- I : Infectious cocoa pods
- R : Removed cocoa pods
- Is : Secondary infection spores
- Ip : Primary infection spores

Total population:

```latex
N(t)=Sc+Sp+Sr+E+I+R+Is+Ip
```

---

## Intervention Strategies

Time-dependent controls include:

### u1 — Cultural practices

- Weed removal
- Sanitation
- Farm maintenance

### u2 — Fungicide application

- Chemical treatment of infected plants

### u3 — Rouging strategy

- Removal of infected pods and infected plant tissues

---

## Methods Implemented

### Mathematical Analysis

- Positivity and boundedness
- Disease-free equilibrium
- Endemic equilibrium
- Basic reproduction number
- Local stability analysis
- Global stability analysis
- Backward bifurcation analysis

### Sensitivity and Uncertainty Analysis

- Normalized Forward Sensitivity Index
- Latin Hypercube Sampling (LHS)
- Partial Rank Correlation Coefficient (PRCC)
- Robustness analysis

### Optimal Control

- Pontryagin Maximum Principle
- Adjoint equations
- Forward–Backward Sweep Method

### Cost-effectiveness Analysis

- Incremental Cost-Effectiveness Ratio (ICER)
- Comparison of intervention strategies

---

## Main Results

Key findings include:

- Black Pod Disease exhibits backward bifurcation behavior.
- Disease elimination may require reducing the reproduction threshold below one.
- Sensitivity analysis identifies critical transmission parameters influencing disease spread.
- Combining intervention strategies substantially reduces disease burden.
- Cultural practices combined with rouging provide the most cost-effective control strategy.

## Citation

If you use this repository, please cite:

Afful, B.A.

*Mitigating Black Pod Disease in Cocoa Plants: An Optimal Control and Cost-Effectiveness Analysis.*

```bibtex
@misc{afful2026BPD,
author={Afful, Bernard Asamoah},
title={Mitigating Black Pod Disease in Cocoa Plants: An Optimal Control and Cost-Effectiveness Analysis},
year={2026}
}
```

---

## Author

Bernard Asamoah Afful
Department of Mathematics and Statistics  
Utah State University  
Logan, Utah, USA
