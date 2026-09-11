# SlyB binder-design research log

## Current validation state

RFD3 single-binder penta backbone generation produced geometry-passing candidates.

LigandMPNN sequence design produced the current candidate set.

RF3 multichain iPTM validation was rejected as an optimisation metric because
the native SlyB positive control itself scored poorly.

## Current next step

Calibrate Boltz2 using the native SlyB A/B positive control before evaluating
the candidate and negative-control binder sequences.

No Boltz2 candidate ranking should be interpreted until the native control
demonstrates that the predictor can recognise the SlyB interface.

## Boltz2 validation variants

Two Boltz2 co-folding protocols are retained separately for provenance.

### Free co-folding baseline

Directory:

    configs/boltz2/calibration_inputs/

This is the earlier sequence-only experiment using the 96-aa soluble SlyB
construct. Chains A, B and binder C are predicted without a structural template.

### 7OJG template-conditioned co-folding

Directory:

    configs/boltz2/template_inputs/

This experiment uses the deposited 155-aa SlyB sequence from PDB 7OJG for
chains A and B.

Chains A and B are structurally conditioned on experimental 7OJG chains A/B:

    force: true
    threshold: 1.0

Binder chain C has no structural template and must find a compatible pose
against the template-conditioned SlyB pair.

The two input sets are intentionally retained separately and should not be
overwritten because they represent different validation experiments.

## vertical60-mpnn-pilot-20260911

Four 60-aa binder backbones were generated with classic RFdiffusion using the SlyB K/A/B/C structural context and balanced A/B hotspot conditioning. ProteinMPNN generated four sequences per backbone; one sequence per backbone was selected for an initial Boltz-2 cofolding pilot. During cofolding, the two native 138-aa SlyB chains are strongly template-conditioned while the 60-aa binder is predicted freely.
