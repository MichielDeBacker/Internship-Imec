from pathlib import Path
import json,re,numpy as np
from scipy.spatial import cKDTree
import check_rfd3_symmetry as cr

cif=Path("out/rfd3_c11_fullsym_smoke_50056893/rfd3_c11_symmetry_smoke_c11_50_single_0_model_0.cif.gz")
js=Path(str(cif)[:-7]+".json")
st=cr.load_structure(cif); meta=json.load(open(js))
ring=cr.load_structure("inputs/target_ring_C11.pdb")
# align merged fixed target back to original ring through source mapping
T=st[st.chain_id=="a"]
out_ca={int(r):x for r,x in zip(T[T.atom_name=="CA"].res_id,T[T.atom_name=="CA"].coord)}
ring_ca={(str(c),int(r)):x for c,r,x in zip(ring[ring.atom_name=="CA"].chain_id,ring[ring.atom_name=="CA"].res_id,ring[ring.atom_name=="CA"].coord)}
def parse(s):
 m=re.fullmatch(r"(.+?)(-?\d+)",str(s)); return m.group(1),int(m.group(2))
mov=[];fix=[]
for src,dst in meta["diffused_index_map"].items():
 sc,sr=parse(src); dc,dr=parse(dst)
 if dc=="a" and dr in out_ca and (sc,sr) in ring_ca:
  mov.append(out_ca[dr]);fix.append(ring_ca[(sc,sr)])
R,mc,fc,rms=cr.kabsch(np.array(mov),np.array(fix))
print("target_align_CA",len(mov),"RMSD",rms)
# aligned target and binders
bbnames=np.array(["N","CA","C","O"])
def aln(x): return cr.apply_transform(x,R,mc,fc)
target_bb=aln(T[np.isin(T.atom_name,bbnames)].coord)
# contact/clearance by original ring chains
ring_parts={c:ring[ring.chain_id==c] for c in np.unique(ring.chain_id)}
rows=[]
for ch in list("ABCDEFGHIJK"):
 b=st[st.chain_id==ch]; bb=b[np.isin(b.atom_name,bbnames)]
 xyz=aln(b.coord); bbxyz=aln(bb.coord)
 tmin=float(cKDTree(target_bb).query(bbxyz)[0].min())
 contacts={c:int((cKDTree(ring_parts[c].coord).query(xyz)[0]<5.0).sum()) for c in ring_parts}
 best=sorted(contacts.items(), key=lambda kv:kv[1], reverse=True)[:3]
 rows.append((ch,tmin,best,xyz,bbxyz,b))
 print("binder",ch,"targetBB",round(tmin,3),"best_contacts",best)
# binder-binder backbone mins, all pairs and adjacent alphabetical cyclic
print("ADJACENT_BB_MINS")
for i,ch in enumerate(list("ABCDEFGHIJK")):
 ch2=list("ABCDEFGHIJK")[(i+1)%11]
 a=rows[i][4]; b=rows[(i+1)%11][4]
 d=float(cKDTree(b).query(a)[0].min())
 print(ch,ch2,round(d,3))
# exact symmetry consistency: compare each binder CA to chain A CA under original ring A->subunit transform
A_ring=ring[(ring.chain_id=="A")&(ring.atom_name=="CA")]
def kabsch_to(ch):
 X=A_ring.coord; Y=ring[(ring.chain_id==ch)&(ring.atom_name=="CA")].coord
 return cr.kabsch(X,Y)
Aout=aln(st[(st.chain_id=="A")&(st.atom_name=="CA")].coord)
print("COPY_RMSD_FROM_A")
for ch in list("ABCDEFGHIJK"):
 Rt,mct,fct,rr=kabsch_to(ch)
 pred=cr.apply_transform(Aout,Rt,mct,fct)
 obs=aln(st[(st.chain_id==ch)&(st.atom_name=="CA")].coord)
 val=float(np.sqrt(np.mean(np.sum((pred-obs)**2,axis=1))))
 print(ch,round(val,4))
# ASU angular footprint using aligned chain A CA
centre,axis=cr.ring_frame(ring); z=axis/np.linalg.norm(axis); trial=np.array([1.,0,0]);
if abs(np.dot(trial,z))>0.9: trial=np.array([0.,1,0])
x=trial-np.dot(trial,z)*z; x/=np.linalg.norm(x); y=np.cross(z,x)
q=Aout-centre; ang=np.mod(np.arctan2(q@y,q@x),2*np.pi); aa=np.sort(ang); gaps=np.diff(np.r_[aa,aa[0]+2*np.pi]); width=2*np.pi-gaps.max()
print("A_CA_angular_footprint_deg",round(float(np.degrees(width)),2),"wedge",round(360/11,2))
