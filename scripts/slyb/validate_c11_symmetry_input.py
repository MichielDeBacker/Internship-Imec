import json
import numpy as np
from rfd3.inference.input_parsing import DesignInputSpecification

path = "inputs/rfd3_c11_symmetry_smoke.json"
with open(path) as fh:
    data = json.load(fh)["c11_50_single"]
spec = DesignInputSpecification(**data)
aa, meta = spec.build(return_metadata=True)
print("metadata", meta["extra"])
print("chains", list(np.unique(aa.chain_id)))
for c in np.unique(aa.chain_id):
    m = aa.chain_id == c
    residues = len(set(map(int, aa.res_id[m])))
    fixed = int(aa.is_motif_atom_with_fixed_coord[m].astype(bool).sum())
    hotspots = int(aa.is_atom_level_hotspot[m].astype(bool).sum())
    tids = sorted(set(map(int, aa.sym_transform_id[m])))
    eids = sorted(set(map(int, aa.sym_entity_id[m])))
    print(c, "atoms", int(m.sum()), "res", residues, "fixed_atoms", fixed, "hotspot_atoms", hotspots, "tids", tids, "eids", eids)
print("transform_ids", {int(x): int((aa.sym_transform_id == x).sum()) for x in np.unique(aa.sym_transform_id)})
print("entity_ids", {int(x): int((aa.sym_entity_id == x).sum()) for x in np.unique(aa.sym_entity_id)})
print("asu_atoms", int(aa.is_sym_asu.astype(bool).sum()))
print("fixed_atoms_total", int(aa.is_motif_atom_with_fixed_coord.astype(bool).sum()))
print("hotspot_atoms_total", int(aa.is_atom_level_hotspot.astype(bool).sum()))
print("src_components_hotspot", sorted(set(aa.src_component[aa.is_atom_level_hotspot.astype(bool)])))
