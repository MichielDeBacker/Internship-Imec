# Prosapia-style SlyB binder pipeline

This directory is an **adapter/workbench layer inside `Internship-Imec`**. It does not edit, patch, vendor, or overwrite Prosapia source code.

The campaign keeps all 150 existing RFD3 backbones and runs two major downstream stages:

1. sequence design with ProteinMPNN, designing only the binder chain;
2. Boltz cofolding with the two SlyB target chains supplied as structural templates.

The canonical chain convention is fixed across every stage:

- `A` = SlyB subunit 1, fixed during sequence design, templated in Boltz
- `B` = SlyB subunit 2, fixed during sequence design, templated in Boltz
- `C` = binder, designed by ProteinMPNN, never templated in Boltz

## Important interpretation

Boltz templates condition the target structure. Even with `force: true` and inference-time potentials, this is not a mathematical coordinate freeze. The pipeline therefore calls A/B "templated" rather than claiming they are exactly immobile.

## Layout

```text
prosapia_binder_pipeline/
├── config/slyb_hydra.yaml
├── prosapia_binder_pipeline/common.py
├── scripts/
│   ├── 00_check_environment.py
│   ├── 01_prepare_backbones.py
│   ├── 02_sequence_one.py
│   ├── 03_collect_sequences.py
│   ├── 04_prepare_boltz_inputs.py
│   ├── 05_boltz_one.py
│   ├── 06_collect_boltz_scores.py
│   ├── 07_status.py
│   ├── submit_sequence.sh
│   └── submit_boltz.sh
└── slurm/
    ├── sequence_array.sbatch
    └── boltz_array.sbatch
```

Runtime artifacts are written outside this source folder to:

```text
work/prosapia_binder_pipeline/slyb_vertical96/
├── 01_backbones/
│   ├── backbones.tsv
│   ├── canonical_pdb/
│   └── target_templates/
├── 02_sequence/
│   ├── runs/
│   └── sequences.tsv
├── 03_boltz/
│   ├── inputs/
│   ├── boltz_inputs.tsv
│   └── runs/
└── 04_results/
    ├── boltz_scores_all_models.csv
    ├── boltz_scores_best.csv
    └── cofold_summary.csv
```

This keeps source code separate from campaign state and makes every stage restartable from a manifest.

A thin `workbench.sh` front end is included, so the same stages can be called as:

```bash
./workbench.sh check
./workbench.sh prepare
./workbench.sh submit-seq
./workbench.sh collect-seq
./workbench.sh prepare-boltz
./workbench.sh submit-boltz
./workbench.sh collect-boltz
./workbench.sh status
```

## 0. Copy into the repository

Place this whole directory at the root of:

```text
/scratch/brussel/vo/000/bvo00014/vsc39230/repos/Internship-Imec/
```

Then:

```bash
cd /scratch/brussel/vo/000/bvo00014/vsc39230/repos/Internship-Imec/prosapia_binder_pipeline
```

The only light Python dependencies used by the adapter scripts are:

```bash
pip install -r requirements.txt
```

If your cluster environments are already provisioned, it is better to point the config at their existing Python/executables rather than reinstalling model packages.

## 1. Edit external tool paths once

Open:

```text
config/slyb_hydra.yaml
```

At minimum verify:

```yaml
external:
  prosapia_root: /path/to/prosapia
  proteinmpnn_root: /path/to/ProteinMPNN
  proteinmpnn_python: /path/to/proteinmpnn/python
  boltz_executable: /path/to/boltz
  boltz_cache: /scratch/.../boltz_cache
```

`prosapia_root` is informational/read-only in this implementation. Nothing in this repository writes into it.

Check the environment:

```bash
python3 scripts/00_check_environment.py --config config/slyb_hydra.yaml
```

## 2. Register all 150 RFD3 backbones

```bash
python3 scripts/01_prepare_backbones.py --config config/slyb_hydra.yaml
```

This stage deliberately requires exactly 150 matches. It refuses to silently continue with a partial set or accidental extra smoke outputs.

Each RFD3 complex is rewritten into a canonical PDB:

```text
A = target partner 1
B = target partner 2
C = binder
```

It also writes a target-only `AB_template.pdb` for each backbone. Original `.cif.gz` files are never modified.

Check:

```bash
wc -l /scratch/brussel/vo/000/bvo00014/vsc39230/repos/Internship-Imec/work/prosapia_binder_pipeline/slyb_vertical96/01_backbones/backbones.tsv
```

Expected: `151` lines = header + 150 backbones.

## 3. Sequence all 150 backbones

The default configuration generates 8 ProteinMPNN sequences per backbone at temperature 0.15. Only chain C is designed. Chains A and B remain fixed context.

Submit:

```bash
THROTTLE=20 bash scripts/submit_sequence.sh
```

The submission helper derives the array size from `backbones.tsv`; the Slurm worker processes exactly one backbone per array task. A failed backbone can therefore be resubmitted independently.

Monitor:

```bash
squeue -u "$USER"
python3 scripts/07_status.py --config config/slyb_hydra.yaml
```

After all 150 sequencing tasks complete:

```bash
python3 scripts/03_collect_sequences.py --config config/slyb_hydra.yaml
```

Default expected result:

```text
150 backbones x 8 sequences = 1200 binder sequences
```

The sequence manifest contains the ProteinMPNN score and the exact backbone/template provenance for every binder sequence.

## 4. Create Boltz cofolding inputs

```bash
python3 scripts/04_prepare_boltz_inputs.py --config config/slyb_hydra.yaml
```

Each YAML is conceptually:

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: TARGET_A
      msa: empty
  - protein:
      id: B
      sequence: TARGET_B
      msa: empty
  - protein:
      id: C
      sequence: DESIGNED_BINDER
      msa: empty

templates:
  - pdb: /absolute/path/to/backbone__AB_template.pdb
    chain_id: [A, B]
    template_id: [A1, B1]
    force: true
    threshold: 1.0
```

A and B are templated. C is intentionally absent from `templates` so the cofold is a test of the designed binder sequence rather than a forced reproduction of its RFD3 backbone.

The default config uses single-sequence mode (`msa: empty`) for reproducibility and to avoid 1200 network calls. To use the Boltz MSA server, change the relevant `*_msa` config values to `server`. Existing MSA files can also be given as paths.

### Number of Boltz jobs

By default:

```yaml
boltz:
  sequences_per_backbone: all
```

so all 1200 generated sequences are prepared for Boltz while preserving all 150 backbones.

If that is too expensive for the first screen, set an integer such as:

```yaml
boltz:
  sequences_per_backbone: 2
```

The preparation stage then keeps the lowest ProteinMPNN-score sequences per backbone, but still requires all 150 backbone identities to remain represented.

## 5. Run Boltz

Submit all prepared inputs:

```bash
THROTTLE=12 bash scripts/submit_boltz.sh
```

The worker uses scratch-backed runtime caches and writes one independent output directory per binder sequence.

By default the Boltz configuration uses:

```yaml
recycling_steps: 3
sampling_steps: 200
diffusion_samples: 1
use_potentials: true
template_force: true
template_threshold_angstrom: 1.0
```

Increase `diffusion_samples` only when you actually want multiple stochastic structure samples per sequence; the collector supports multiple models.

## 6. Collect the cofolding ipTM table

After all Boltz array tasks finish:

```bash
python3 scripts/06_collect_boltz_scores.py --config config/slyb_hydra.yaml
```

The compact output is:

```text
work/prosapia_binder_pipeline/slyb_vertical96/04_results/cofold_summary.csv
```

with the requested columns:

```text
name
role
iptm
AB_pair_iptm
AC_pair_iptm
BC_pair_iptm
binder_mean
binder_min
binder_delta
```

where:

```text
AB = SlyB subunit 1 vs SlyB subunit 2
AC = SlyB subunit 1 vs binder
BC = SlyB subunit 2 vs binder
binder_mean  = mean(AC, BC)
binder_min   = min(AC, BC)
binder_delta = abs(AC - BC)
```

Boltz reports `pair_chains_iptm` directionally, so `A->C` and `C->A` may differ. The pipeline therefore retains every directional value in `boltz_scores_best.csv`. The compact `AC_pair_iptm`, `BC_pair_iptm`, and `AB_pair_iptm` use `scoring.pair_summary_mode` from the config. Default is `max`.

Available modes are:

```text
max
mean
min
forward
```

This makes the convention explicit instead of hiding Boltz's asymmetric raw values.

## 7. Safe restart behavior

Both expensive worker stages are idempotent at the campaign level:

- sequencing skips a backbone when the expected FASTA already contains the requested number of samples;
- Boltz skips an input when the expected confidence JSONs already exist;
- original RFD3 `.cif.gz` files are never modified;
- Prosapia source is never modified.

To intentionally rerun one Boltz input:

```bash
python3 scripts/05_boltz_one.py \
  --config config/slyb_hydra.yaml \
  --index 17 \
  --override
```

## 8. Git workflow

Only commit the source adapter folder, not the large work outputs:

```bash
cd /scratch/brussel/vo/000/bvo00014/vsc39230/repos/Internship-Imec
git status
git add prosapia_binder_pipeline
git commit -m "Add Prosapia-style SlyB sequencing and Boltz cofold pipeline"
git push
```

The `work/` outputs should remain scratch/campaign data unless you explicitly decide to version small summary tables.

## Upstream interfaces used

- ProteinMPNN: `protein_mpnn_run.py --pdb_path ... --pdb_path_chains C` makes C the designed chain while the other chains are fixed.
- Boltz: YAML supports protein sequences and structural templates with explicit `chain_id` / `template_id` mapping; prediction confidence JSON contains `iptm` and `pair_chains_iptm`.

This adapter intentionally calls those public interfaces instead of modifying upstream Prosapia, ProteinMPNN, or Boltz source code.
