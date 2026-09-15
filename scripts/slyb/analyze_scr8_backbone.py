from pathlib import Path
import json,numpy as np
from scipy.spatial import cKDTree
import check_rfd3_symmetry as cr
root=Path('out/rfd3_penta_single_36_52_hotspotori_scr8_50056908')
template=cr.load_structure('inputs/target_penta_soluble.pdb')
ring=cr.load_structure('inputs/target_ring_C11.pdb')
CEN=np.array([112.891,112.867,87.187]); AXIS=np.array([-0.000390,0.000302,1.0]); AXIS/=np.linalg.norm(AXIS)
e1=np.cross(AXIS,[1,0,0]); e1/=np.linalg.norm(e1); e2=np.cross(AXIS,e1)
ringbb=ring[np.isin(ring.atom_name,['N','CA','C','O'])]
def camap(st,ch):
 a=st[(st.chain_id==ch)&(st.atom_name=='CA')]; return {int(r):x for r,x in zip(a.res_id,a.coord)}
def kab(mov,fix):
 ks=sorted(set(mov)&set(fix)); P=np.array([mov[k] for k in ks]); Q=np.array([fix[k] for k in ks]); pm=P.mean(0); qm=Q.mean(0)
 U,S,Vt=np.linalg.svd((P-pm).T@(Q-qm)); R=U@Vt
 if np.linalg.det(R)<0: U[:,-1]*=-1; R=U@Vt
 return R,pm,qm
A=camap(ring,'A'); T={c:kab(A,camap(ring,c)) for c in ['B','K']}
rows=[]
for cif in sorted(root.glob('*.cif.gz')):
 meta=json.load(open(str(cif)[:-7]+'.json')); st=cr.load_structure(cif); tc,R,mc,fc,rmsd,n=cr.build_alignment(st,meta,template)
 bc=[c for c in np.unique(st.chain_id) if c!=tc][0]; b=st[st.chain_id==bc]; bb=b[np.isin(b.atom_name,['N','CA','C','O'])]; xyz=cr.apply_transform(bb.coord,R,mc,fc)
 d=xyz-CEN; th=np.degrees(np.arctan2(d@e2,d@e1)); th=(th-th.mean()+180)%360-180; fp=float(th.max()-th.min())
 r={'name':cif.name,'len':len(set(b.res_id)),'rmsd':float(rmsd),'fp':fp,'fp_pass':fp<=30,'PASS':False}
 if r['fp_pass']:
  r['target']=float(cKDTree(ringbb.coord).query(xyz)[0].min())
  for c,(Rt,pm,qm) in T.items(): r[c]=float(cKDTree((xyz-pm)@Rt+qm).query(xyz)[0].min())
  allxyz=cr.apply_transform(b.coord,R,mc,fc); r['ca']=cr.contact_count(allxyz,template[template.chain_id=='A']); r['cb']=cr.contact_count(allxyz,template[template.chain_id=='B']); r['zmax']=float(((xyz-CEN)@AXIS).max())
  r['PASS']=r['target']>=2.2 and r['B']>=2.2 and r['K']>=2.2 and r['ca']>5 and r['cb']>5 and r['zmax']<=22
 rows.append(r)
for i,r in enumerate(rows): print(i,json.dumps(r,sort_keys=True))
print('SUMMARY',json.dumps({'n':len(rows),'fp_pass':sum(x['fp_pass'] for x in rows),'full_pass':sum(x['PASS'] for x in rows),'fps':[round(x['fp'],3) for x in rows]},sort_keys=True))
