# SlyB binder pipeline status

This repository now contains the orchestration layer for continuing
the SlyB binder-design campaign on VUB Hydra.

The runtime data remain under the project scratch directory and are
not copied into Git.

## Workflow

1. `doctor`
2. `inventory`
3. `geometry`
4. review viable parent geometry
5. `rfd3-plan`
6. install/activate RFdiffusion3 on a GPU node
7. authorize a bounded RFD3 generation run
8. sequence design
9. RF3/Boltz validation
10. full C11 geometry and ranking

## Primary recovered candidates

- d2_s0
- d1_s3
- d6_s6

Reserve / challenge controls:

- d2_s1
- d6_s3
- d6_s7

## Existing parent backbone families

- design0
- design1
- design2
- design6

New RFdiffusion3 generation is intentionally gated until the existing
C11 geometry has been screened and target hotspots are reviewed.
