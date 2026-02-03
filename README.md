# ✈️ Robust Fault Detection in Jet Engine Systems
### *Adaptive Learning, Estimation and Supervision of Dynamical Systems - 2025*

## 📄 Overview

This project focuses on the implementation and validation of **Unknown Input Observers (UIO)** combined with **Fault Detection Filters (BFDF)**. The work reproduces the methodology presented in the paper *"Design of Unknown Input Observers and Robust Fault Detection Filters"* (Chen, Patton, Zhang), applying it to a non-linear **Jet Engine model**.

The primary objective is to generate **robust directional residuals** that allow for fault isolation even in the presence of modelling uncertainties (linearization errors) treated as unknown disturbances.

## 🎯 Project Goals

* **Simulate** the non-linear dynamics of a Jet Engine system.
* **Design** a full-order Unknown Input Observer (UIO) to decouple disturbances.
* **Implement** a robust residual generator with directional properties for fault isolation.
* **Reproduce** the results of Table 1 (Eigenvalues) and Table 2 (Fault Isolation Logic) from the reference paper.

## 🧠 Methodology & Theory

The project follows a rigorous control theory pipeline derived from the paper's main theorems:

### 1. Disturbance Decoupling (UIO)
We model the system with unknown inputs $d(t)$ (representing linearization errors):
$$\dot{x}(t) = Ax(t) + Bu(t) + Ed(t)$$
$$y(t) = Cx(t)$$

The observer is designed to decouple $d(t)$ by satisfying the existence conditions:
* **Rank Condition:** $\text{rank}(CE) = \text{rank}(E)$
* **Detectability:** the pair $(C, A_1)$ must be detectable, where $A_1 = A - E[(CE)^T CE]^{-1}(CE)^T CA$.

### 2. Directional Residuals (BFDF)
Once robustness against disturbances is achieved, the remaining design freedom in the gain matrix $K$ is exploited to enforce **Directional Residuals**.
* The residual vector $r(t)$ is designed to lie in a fixed, fault-specific subspace.
* Faults are isolated by comparing $r(t)$ against pre-computed fault signatures.

## 💻 Implementation Details

The solution is implemented in **Python** using `numpy` and `scipy.signal`. The pipeline consists of:

1.  **System Linearization:** approximating the non-linear Jet Engine model at the operating point.
2.  **Matrix Synthesis:**
    * Computation of $H^* = E[(CE)^T CE]^{-1}(CE)^T$ (Eq. 12).
    * Stabilization of $F = A_1 - K_1 C$.
3.  **Simulation:** numerical integration of the system dynamics with injected sensor/actuator faults.

## 🛠️ Tech Stack
* **Language:** MATLAB
* **Reference Paper:** Chen, J., Patton, R. J., & Zhang, H. Y. *Design of Unknown Input Observers and Robust Fault Detection Filters*.

## 👨‍💻 Team Members

| Nome | StudentID | GitHub |
|------|-----------|--------|
| Daniele Gotti | 1079011 | [@DanieleGotti](https://github.com/DanieleGotti) |
| Filippo Antonio Bolis | 1079493 | [@BolisFilippo](https://github.com/BolisFilippo) |
