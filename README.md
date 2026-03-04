# VCG-CAFNCB

Variational Curvature Governance Framework
Engineering Release (v1.0)

---

## Overview

VCG-CAFNCB is a multi-generation optimization framework for stabilizing high-dimensional systems, evolving from PNN (Gen-1) macro computations to full VCG (Gen-4) with geometric governance.

Core mechanisms:
- Macro-micro unification (NCB Gen-2)
- Cross-region fusion (CAFNCB Gen-3)
- Curvature-triggered governance (kappa signal)
- One-step lag smoothing for stability
- Willow rollback for deviation control
- Density-aware scaling with hit gating
- Structural consistency enforcement (SEC)

Focus: Numerical stability, reproducibility, and auditability in iterative processes.

---

## Architecture

PNN (Gen-1) → NCB (Gen-2) → CAFNCB (Gen-3) → VCG-CAFNCB (Gen-4)

VCG introduces:
1. Curvature governance for non-linear control
2. Rollback stabilization via Willow
3. Density gating at state/increment levels
4. Smoothed signals for robust updates

---

## Installation

Requires MATLAB (tested on R2023b+). No external dependencies beyond base toolboxes.

Clone repo:
git clone https://huggingface.co/<your-username>/vcg-cafncb

---

## How to Run

Navigate to matlab/ folder and run:
results = Copy_of_compare_all();

View metrics:
disp(results);

Expected structure (approximate):
NCB_MSE, CAF_MSE, VCG_MSE
NCB_qecr, CAF_qecr, VCG_qecr
VCG_willow_w, VCG_willow_rho

---

## Stability Features (v1.0)

- Cold-start guards prevent invalid audits
- Denominator regularization (1e-3 min for division stability)
- Kappa smoothing with lag (configurable eta_smooth)
- Density governance with iterative scaling (shrink/expand factors)
- Hit-based activation to avoid over-governance
- Rollback gamma (0.2 default) for controlled mixing

---

## License

MIT (see LICENSE file)

---

## Citation

If using this framework, cite as:

@software{vcg_cafncb_2026,
  author = {Your Name},
  title = {VCG-CAFNCB: Variational Curvature Governance Framework},
  version = {1.0},
  year = {2026},
  url = {https://huggingface.co/<your-username>/vcg-cafncb}
}