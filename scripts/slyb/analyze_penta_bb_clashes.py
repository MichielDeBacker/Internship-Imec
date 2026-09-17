from pathlib import Path
import json, numpy as np
from scipy.spatial import cKDTree
import check_rfd3_symmetry as cr

root=Path('.')
ring=cr.load_structure('inputs/target_ring_C11.pdb')
template=cr.load_structure('inputs/target_penta_soluble.pdb')

def ca_map(st,ch):
    a=st[(st.chain_id==ch)&(st.atom_name=='CA')]
    return {int(r):x for r,x in zip(a.res_id,a.coord)}
def kab(moving,fixed):
    keys=sorted(set(moving)&set(fixed))
    P=np.array([moving[k] for k in keys]); Q=np.array([fixed[k] for k in keys])
    pm=P.mean(0); qm=Q.mean(0)
    U,S,Vt=np.linalg.svd((P-pm).T@(Q-qm)); R=U@Vt
    if np.linalg.det(R)<0:
        U[:,-1]*=-1; R=U@Vt
    return R,pm,qm,keys
A=ca_map(ring,'A')
trans={}
for ch in ['B','K']:
    R,mc,fc,keys=kab(A,ca_map(ring,ch)); trans[ch]=(R,mc,fc)
    P=np.array([A[k] for k in keys]); Q=np.array([ca_map(ring,ch)[k] for k in keys])
    P2=(P-mc)@R+fc
    print('A->'+ch,'RMSD',np.sqrt(np.mean(np.sum((P2-Q)**2,axis=1))))

bbnames=np.array(['N','CA','C','O'])
for arm in ['penta_50_single','penta_60_single']:
    rows=[]
    for cif in sorted((root/'out/rfd3_T200_penta_single_screen'/arm).glob('*.cif.gz')):
        meta=json.load(open(str(cif)[:-7]+'.json'))
        st=cr.load_structure(cif)
        target_chain,R,mc,fc,rmsd,nmatch=cr.build_alignment(st,meta,template)
        bchains=[c for c in np.unique(st.chain_id) if c!=target_chain]
        if len(bchains)!=1: continue
        b=st[st.chain_id==bchains[0]]
        b=b[np.isin(b.atom_name,bbnames)]
        xyz=cr.apply_transform(b.coord,R,mc,fc)
        best=None
        for ch,(Rt,mct,fct) in trans.items():
            nbr=(xyz-mct)@Rt+fct
            tree=cKDTree(nbr); d,jj=tree.query(xyz)
            i=int(np.argmin(d)); val=float(d[i]); j=int(jj[i])
            rec=(val,ch,i,j)
            if best is None or val<best[0]: best=rec
        val,ch,i,j=best
        rows.append((val,cif.name,int(b.res_id[i]),str(b.atom_name[i]),int(b.res_id[j]),str(b.atom_name[j]),ch))
    rows.sort(reverse=True)
    vals=[x[0] for x in rows]
    print('\n'+arm,'N',len(rows))
    print('BB neighbor min: min/median/max',round(min(vals),3),round(float(np.median(vals)),3),round(max(vals),3))
    print('BEST 8 largest clearance')
    for r in rows[:8]:
        print(round(r[0],3),r[1],'self',f'{r[2]}:{r[3]}','vs',r[6],f'{r[4]}:{r[5]}')

print('\n=== FOOTPRINT ANALYSIS ===')
# PCA-derived ring frame from existing checker
centre,axis=cr.ring_frame(ring)
# construct orthonormal plane basis
z=axis/np.linalg.norm(axis)
trial=np.array([1.0,0,0])
if abs(np.dot(trial,z))>0.9: trial=np.array([0,1.0,0])
x=trial-np.dot(trial,z)*z; x=x/np.linalg.norm(x); y=np.cross(z,x)
for arm in ['penta_50_single','penta_60_single']:
    vals=[]
    for cif in sorted((root/'out/rfd3_T200_penta_single_screen'/arm).glob('*.cif.gz')):
        meta=json.load(open(str(cif)[:-7]+'.json')); st=cr.load_structure(cif)
        target_chain,R,mc,fc,rmsd,nmatch=cr.build_alignment(st,meta,template)
        bchains=[c for c in np.unique(st.chain_id) if c!=target_chain]
        b=st[st.chain_id==bchains[0]]; bb=b[np.isin(b.atom_name,bbnames)]
        xyz=cr.apply_transform(bb.coord,R,mc,fc)
        # backbone target clearance to complete ring
        tree=cKDTree(ring[np.isin(ring.atom_name,bbnames)].coord); dt=float(tree.query(xyz)[0].min())
        ca_xyz=cr.apply_transform(b[b.atom_name=='CA'].coord,R,mc,fc)
        q=ca_xyz-centre
        ang=np.mod(np.arctan2(q@y,q@x),2*np.pi)
        a=np.sort(ang); gaps=np.diff(np.r_[a,a[0]+2*np.pi]); width=2*np.pi-gaps.max()
        widthdeg=np.degrees(width)
        vals.append((widthdeg,dt,cif.name))
    widths=np.array([v[0] for v in vals]); dts=np.array([v[1] for v in vals])
    print(arm,'CA angular footprint deg min/median/max',round(widths.min(),1),round(np.median(widths),1),round(widths.max(),1),'C11 wedge',round(360/11,1))
    print(arm,'BB target min A min/median/max',round(dts.min(),3),round(np.median(dts),3),round(dts.max(),3),'BB-clear>=2.2',int((dts>=2.2).sum()),'/',len(dts))
    print('smallest footprints')
    for w,d,n in sorted(vals)[:5]: print(round(w,1),'deg targetBB',round(d,3),n)
