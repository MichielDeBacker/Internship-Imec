#!/usr/bin/env python3
from pathlib import Path
import argparse, csv, gzip, json
from collections import Counter
import numpy as np
from scipy.spatial import cKDTree
import check_rfd3_symmetry as old

BB = np.array(["N","CA","C","O"])
HARD = 2.2
CONTACT = 5.0
MIN_CONTACTS = 20
MAX_ALIGN = 0.5

def bb(a): return a[np.isin(a.atom_name, BB)]
def mdist(x,y): return float(cKDTree(y).query(x)[0].min())
def contacts(x,y): return int((cKDTree(y).query(x)[0] < CONTACT).sum())
def ca_map(a,ch):
    q=a[(a.chain_id==ch)&(a.atom_name=="CA")]
    return {int(r):x for r,x in zip(q.res_id,q.coord)}
def fit_map(moving,fixed):
    ks=sorted(set(moving)&set(fixed)); P=np.array([moving[k] for k in ks]); Q=np.array([fixed[k] for k in ks])
    pm=P.mean(0); qm=Q.mean(0); U,S,Vt=np.linalg.svd((P-pm).T@(Q-qm)); R=U@Vt
    if np.linalg.det(R)<0: U[:,-1]*=-1; R=U@Vt
    return R,pm,qm
def xform(x,t):
    R,pm,qm=t; return (x-pm)@R+qm

def one(cif, ring, penta, tnext, tprev):
    meta=json.load(open(str(cif)[:-7]+".json"))
    st=old.load_structure(cif)
    tc,R,mc,fc,ar,nm=old.build_alignment(st,meta,penta)
    binds=[str(c) for c in np.unique(st.chain_id) if c!=tc]
    if len(binds)!=2: raise ValueError(f"expected 2 binder chains, got {binds}")
    xyz={}; xyzbb={}; lens={}
    for b in binds:
        a=st[st.chain_id==b]; lens[b]=len(set(map(int,a.res_id)))
        xyz[b]=old.apply_transform(a.coord,R,mc,fc)
        bba=bb(a); xyzbb[b]=old.apply_transform(bba.coord,R,mc,fc)
    # assignment on all-atom contacts to penta chains K/A/B
    cc={(b,c):contacts(xyz[b],penta[penta.chain_id==c].coord) for b in binds for c in ["K","A","B"]}
    b0,b1=binds
    options=[{b0:("K","A"),b1:("A","B")},{b0:("A","B"),b1:("K","A")}]
    assignment=max(options,key=lambda a:sum(cc[(b,c)] for b,p in a.items() for c in p))
    bridge={b:all(cc[(b,c)]>=MIN_CONTACTS for c in assignment[b]) for b in binds}
    ringbb=bb(ring).coord
    targetmin={b:mdist(xyzbb[b],ringbb) for b in binds}
    pairmin=mdist(xyzbb[b0],xyzbb[b1])
    selfmin={}
    for b in binds:
        selfmin[b]=min(mdist(xyzbb[b],xform(xyzbb[b],tnext)),mdist(xyzbb[b],xform(xyzbb[b],tprev)))
    # compactness diagnostic around ring axis
    centre,axis=old.ring_frame(ring); z=axis/np.linalg.norm(axis); trial=np.array([1.,0,0])
    if abs(trial@z)>0.9: trial=np.array([0.,1,0])
    ex=trial-(trial@z)*z; ex/=np.linalg.norm(ex); ey=np.cross(z,ex)
    footprints={}
    for b in binds:
        ca=st[(st.chain_id==b)&(st.atom_name=="CA")]; q=old.apply_transform(ca.coord,R,mc,fc)-centre
        a=np.sort(np.mod(np.arctan2(q@ey,q@ex),2*np.pi)); gaps=np.diff(np.r_[a,a[0]+2*np.pi]); footprints[b]=float(np.degrees(2*np.pi-gaps.max()))
    passed=(ar<MAX_ALIGN and all(bridge.values()) and all(v>=HARD for v in targetmin.values()) and pairmin>=HARD and all(v>=HARD for v in selfmin.values()))
    out={"design":str(cif),"PASS":passed,"target_alignment_CA_RMSD_A":round(ar,4),"alignment_CA_count":nm,
         "binder_chains":",".join(binds),"assignment":";".join(f"{b}->{assignment[b][0]}/{assignment[b][1]}" for b in binds),
         "both_bridge":all(bridge.values()),"binder_binder_BB_min_A":round(pairmin,3),"binder_binder_BB_clash":pairmin<HARD}
    for b in binds:
        out.update({f"{b}_length":lens[b],f"{b}_bridges":bridge[b],f"{b}_target_BB_min_A":round(targetmin[b],3),
                    f"{b}_target_BB_clash":targetmin[b]<HARD,f"{b}_self_repeat_BB_min_A":round(selfmin[b],3),
                    f"{b}_self_repeat_BB_clash":selfmin[b]<HARD,f"{b}_CA_footprint_deg":round(footprints[b],1)})
        for c in ["K","A","B"]: out[f"contacts_{b}_{c}"]=cc[(b,c)]
    return out

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--root",required=True); ap.add_argument("--json-out"); ap.add_argument("--csv-out"); args=ap.parse_args()
    ring=old.load_structure("inputs/target_ring_C11.pdb"); penta=old.load_structure("inputs/target_penta_soluble.pdb")
    A=ca_map(ring,"A"); tnext=fit_map(A,ca_map(ring,"B")); tprev=fit_map(A,ca_map(ring,"K"))
    rows=[]
    for cif in sorted(Path(args.root).glob("**/*.cif.gz")):
        try: rows.append(one(cif,ring,penta,tnext,tprev))
        except Exception as e: rows.append({"design":str(cif),"PASS":False,"error":repr(e)})
    print(json.dumps(rows,indent=2))
    if args.json_out: json.dump(rows,open(args.json_out,"w"),indent=2)
    if args.csv_out and rows:
        keys=sorted(set().union(*(r.keys() for r in rows))); w=csv.DictWriter(open(args.csv_out,"w",newline=""),fieldnames=keys); w.writeheader(); w.writerows(rows)
if __name__=="__main__": main()
