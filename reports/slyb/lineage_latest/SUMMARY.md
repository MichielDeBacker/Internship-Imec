# SlyB lineage recovery

Native SlyB chain A: 138 residues
Historical RFD3 structures scanned: 117
Cofold samples indexed: 30

## Parent representations

- design0: chains [50, 480] -> **ambiguous_two_chain**
- design1: chains [37, 480] -> **ambiguous_two_chain**
- design2: chains [51, 480] -> **ambiguous_two_chain**
- design6: chains [38, 480] -> **ambiguous_two_chain**

## RF3 monomer -> parent matches

- d2_s0 -> design2 chain A: CA RMSD=0.41077 Å, seq identity=0.45098
- d2_s1 -> design2 chain A: CA RMSD=0.36538 Å, seq identity=0.529412
- d1_s3 -> design1 chain A: CA RMSD=0.47842 Å, seq identity=0.378378
- d6_s6 -> design6 chain A: CA RMSD=0.67426 Å, seq identity=0.447368
- d6_s3 -> design6 chain A: CA RMSD=11.29454 Å, seq identity=0.368421
- d6_s7 -> design6 chain A: CA RMSD=9.48374 Å, seq identity=0.447368

## Best cofold samples

- d2_s0: ranking=0.1693, iPTM=0.13800372183322906, chains=96,96,51, binder=C, clash=False
- d2_s1: ranking=0.1722, iPTM=0.1458863615989685, chains=96,96,51, binder=C, clash=False
- d1_s3: ranking=0.1457, iPTM=0.12190351635217667, chains=96,96,37, binder=C, clash=False
- d6_s3: ranking=0.1709, iPTM=0.1561616212129593, chains=96,96,38, binder=C, clash=False
- d6_s7: ranking=0.2125, iPTM=0.16949540376663208, chains=96,96,38, binder=C, clash=False

## Next decision

Use historical_rfd3_matches.tsv to restore the original target/input specification before authorizing new RFD3.
