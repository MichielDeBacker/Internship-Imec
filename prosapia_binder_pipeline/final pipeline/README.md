# Final pipeline - SlyB interface binder design

This directory contains the retained candidate sequences and the files
needed to continue the final structural-selection stage.

## Goal

Design a 58-62 aa de novo binder that spans the interface between two
adjacent SlyB subunits on the outside/soluble side of the C11 ring.

## RFD3 target representation

Each target subunit contains 96 residues:

- native residues 18-60: 43 aa
- native residues 103-155: 53 aa

Native residues 61-102 are deliberately excluded from the design target.

RFD3 binder contig:

`58-62,/0,B18-60,B103-155,/0,A18-60,A103-155`

Two hotspot scenarios were generated:

- A110 / B132 / B155
- A27 / B132 / B155

ProteinMPNN sampled 8 binder sequences per backbone, giving 800 sequences.

## Final Boltz context

The production cofold uses four neighbouring soluble target fragments:

- Boltz A = native K96
- Boltz B = native A96
- Boltz C = native B96
- Boltz D = native C96
- Boltz E = binder

The intended central interface is native A/B.
Native K/C provide neighbouring steric context.

## Production result

800/800 predictions completed successfully.

16/800 had both:

- binder->native A pairwise iPTM > 0.50
- binder->native B pairwise iPTM > 0.50

14/800 also had the binder centroid on the radial outside of the ring.

These 14 form the BROAD shortlist in this directory.

They are not all validated dual-interface binders.

Several high-iPTM models have zero physical contact with one of the two
central target subunits. Pairwise iPTM must therefore never be used alone.

## Files

### sequences/

`broad14_binder_sequences.fasta`
- exact binder sequences for all 14 broad-shortlist candidates.

`priority5_binder_sequences.fasta`
- five candidates currently considered particularly useful for close
  structural inspection.

### metadata/

`broad14_shortlist.tsv`
- sequence, pairwise iPTM, contact counts, outside-ring metric, alignment
  RMSD and full-ring severe-clash count.

`priority5.txt`
- current visual-inspection priority set.

`full_production_manifest_800.tsv`
- provenance mapping for the complete 800-design production campaign.

### inputs/boltz_4x96/

Exact Boltz YAML inputs for the 14 retained candidates.

### results/confidence/

Original Boltz confidence JSON files.

### results/structures/

Original Boltz model_0 CIF structures.

### results/full_ring_viewers/

Full-native-C11 viewers, when generated.

Viewer convention:

- native chains A-K = unchanged complete C11 ring
- chain Z = transformed Boltz binder

### reference/

`target_ring_C11.pdb`
- complete native reference ring.

`native_KABC_as_ABCD_96aa.pdb`
- validated 4x96 K/A/B/C template used for Boltz conditioning.

### parameters/

Exact known Boltz production parameters and chain mappings.

## Current broad shortlist

The strongest broad-shortlist candidates by minAB include:

- A110_bb044_s01
- A27_bb010_s04
- A27_bb010_s03
- A27_bb001_s02
- A27_bb047_s02
- A27_bb015_s02
- A27_bb027_s02
- A27_bb004_s03
- A27_bb031_s04
- A27_bb015_s06
- A110_bb033_s08
- A110_bb032_s08
- A27_bb015_s01
- A27_bb010_s08

Candidates currently deserving particularly close geometric inspection:

- A27_bb047_s02
- A110_bb032_s08
- A27_bb001_s02
- A110_bb044_s01
- A27_bb010_s04

## Important remaining work

Before experimental selection:

1. visually inspect all broad-shortlist full-ring models;
2. require meaningful physical contact with both native A and B;
3. reject obvious K/C binding or pore-facing poses;
4. inspect every <2 A clash atom-by-atom;
5. distinguish backbone clashes from side-chain rotamer clashes;
6. calculate scenario-specific hotspot distances;
7. evaluate interface PAE / confidence;
8. calculate buried interface area / SASA;
9. structurally relax the surviving complexes;
10. perform an energetic interface assessment.

Do not select candidates from pairwise iPTM alone.
