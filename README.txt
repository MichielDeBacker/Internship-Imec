SlyB vertical-groove RFD3 overnight package

Files:
  pipeline/lib/build_vertical_groove_spec.py
  pipeline/lib/screen_vertical_rfd3.py
  pipeline/slurm/rfd3_vertical_overnight.sbatch
  pipeline/slurm/submit_vertical.sh

Install on Hydra from the Internship-Imec repo root by extracting this archive over the repo.
Then run:
  bash pipeline/slurm/submit_vertical.sh

The batch job generates up to 192 new RFD3 backbones (3 input arms x 16 batches x 4 samples), then screens them for:
  - vertical principal binder axis <=20 degrees from native C11 symmetry axis
  - target fit <=2.0 A
  - symmetry placement residual <=1.0 A
  - >=3 distinct binder residues contacting each adjacent partner at 5 A for every copy
  - binder-ring backbone minimum >=1.8 A
  - all 55 binder-binder copy pairs minimum >=1.8 A

Passing structures are classified NEEDS_VISUAL_REVIEW, never automatically PASS.
