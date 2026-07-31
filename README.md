# Mitigating Overconfidence in Nonlinear Kalman Filters via Covariance Recalibration

This is the code for the research paper with the same title, accepted for publication in *Automatica* in 2026. A preprint is available at <https://arxiv.org/abs/2407.05717>.

In the paper, we introduce a new covariance-recalibrated framework that can reduce the state estimation errors of various types of nonlinear Kalman filters by more than an order of magnitude. The nonlinear Kalman filters investigated in the paper include the extended Kalman filter, second-order extended Kalman filter, unscented Kalman filter, and cubature Kalman filter.

Running these codes allows you to reproduce all the figures in the paper. The five MATLAB scripts in the root directory correspond to the five applications studied in the paper:

| Script | Application |
| --- | --- |
| `Target_tracking.m` | 3D target tracking |
| `Terrain_referenced_navigation.m` | Terrain-referenced navigation |
| `Synchronous_generator_state_estimation.m` | Synchronous generator state estimation |
| `Pendulum_state_estimation.m` | Pendulum state estimation |
| `Battery_state_estimation.m` | Battery state estimation |

`Supplementary_Material.pdf` contains the detailed algorithm listings for the EKF, EKF2, CKF, UKF, and the conventional IEKF benchmark, together with additional results for the "back out" step.

## Notes on the code

The iterated EKF2 (IEKF2) is implemented and simulated in the scripts, but its curves are not drawn, because the figures in the paper do not include them. To display the IEKF2 results, uncomment the corresponding plotting lines (the `h6_1` handle) and the associated runtime `disp` lines.

## Updates

**2025-02-04:** In UKF and CKF, `sqrtm(Variance)` is now replaced by `chol(Variance).'`, making the algorithm faster.

**2025-11-25:** We conducted additional experiments to validate the need for the "back out" step. The related code can be found in the folder `Necessity_of_back_out`. Additionally, the newly uploaded `Target_tracking_ANEE.m` compares the accuracy of covariance estimation between the old and new frameworks.

**2026-07-31:** The paper has been accepted by *Automatica*.

## How to cite

```bibtex
@article{jiang2026mitigating,
  title={Mitigating Overconfidence in Nonlinear {K}alman Filters via Covariance Recalibration},
  author={Jiang, Shida and Shi, Junzhe and Moura, Scott},
  journal={Automatica},
  year={2026}
}
```
