#!/usr/bin/env python3
"""Source-only ring-plane contour QA; atom projections do not replace 3D sampling."""
import pathlib,json,math
import gemmi,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=pathlib.Path('data/archive/heldout_reference');(OUT/'density_slices').mkdir(exist_ok=True)
def grid(p,k):
 b=gemmi.as_refln_blocks(gemmi.cif.read_file(str(OUT/f'{p.lower()}_validation_{k}_map_coef.cif.gz')))[0]
 f,ph=('pdbx_FWT','pdbx_PHWT') if k=='2fo-fc' else ('pdbx_DELFWT','pdbx_DELPHWT')
 g=b.transform_f_phi_to_map(f,ph,sample_rate=3);a=np.asarray(g);a[:]=(a-a.mean())/a.std();return g
for c in json.load(open(OUT/'source_diagnostics.json'))['cases']:
 p=c['case_id'];g=grid(p,'2fo-fc');dg=grid(p,'fo-fc');panels=[]
 for r in c['residues']:
  for alt in r['alternate_ids'] or ['']:panels.append((r,alt))
 n=len(panels);cols=min(3,n);rows=math.ceil(n/cols);fig,axs=plt.subplots(rows,cols,figsize=(5.0*cols,4.7*rows+.4),squeeze=False)
 for ax,(r,alt) in zip(axs.flat,panels):
  atoms=[a for a in r['atom_map'] if a['alt_id'] in ['',alt]];x=np.array([a['xyz_A'] for a in atoms]);names=[a['atom_name'] for a in atoms];by=dict(zip(names,x));comp=r['ccd_id'];ring=(['C1','C2','C3','C4','C5','C6'] if comp=='OPO' else ['C2','C3','C4','C5','O5'] if comp=='FRU' else ['C1','C2','C3','C4','C5','O5']);rx=np.array([by[a] for a in ring]);center=rx.mean(0);u,s,v=np.linalg.svd(rx-center);basis=v[:2];pr=(x-center)@basis.T;lim=max(3.5,np.max(np.abs(pr))+.8);a=np.arange(-lim,lim+.01,.08);xx,yy=np.meshgrid(a,a);z=np.empty_like(xx);dz=np.empty_like(xx)
  for i in range(len(a)):
   for j in range(len(a)):
    pos=gemmi.Position(*(center+xx[i,j]*basis[0]+yy[i,j]*basis[1]));z[i,j]=g.interpolate_value(pos);dz[i,j]=dg.interpolate_value(pos)
  ax.imshow(z,extent=(-lim,lim,-lim,lim),origin='lower',cmap='Blues',vmin=0,vmax=3,alpha=.42);ax.contour(xx,yy,z,levels=[1],colors=['#1762ab'],linewidths=.8)
  if dz.min()<=-3:ax.contour(xx,yy,dz,levels=[-3],colors=['#e04a32'],linewidths=1)
  if dz.max()>=3:ax.contour(xx,yy,dz,levels=[3],colors=['#279b4b'],linewidths=1)
  keys={a['native_atom_key']:i for i,a in enumerate(atoms)}
  for e in c['native_bond_graph']:
   if e['atom1'] in keys and e['atom2'] in keys:
    ij=[keys[e['atom1']],keys[e['atom2']]];ax.plot(pr[ij,0],pr[ij,1],color='#555',lw=1,zorder=2)
  ax.scatter(pr[:,0],pr[:,1],s=15,c=[{'C':'#555','O':'#d13e2e','N':'#2c55b8'}[a['element']] for a in atoms],zorder=3)
  for i,name in enumerate(names):ax.annotate(name,pr[i],fontsize=6,xytext=(2,2),textcoords='offset points')
  supported=sum(a['two_mFo_DFc_max_within_0p5A_sigma']>=1 for a in atoms);neg=min(a['mFo_DFc_atom_sigma'] for a in atoms);rscc=r['validation'].get('rscc');fit=f'{rscc:.3f}' if rscc is not None else 'unknown'
  ax.set_title(f"{r['deposited']['label_chain']}:{r['branch_num']} {comp} {alt} | RSCC {fit}\n3D supported atoms {supported}/{len(atoms)}; min diff {neg:.2f}σ",fontsize=8.5);ax.set_aspect('equal');ax.set_xlabel('ring-plane coordinate (Å)',fontsize=8);ax.set_ylabel('Å',fontsize=8);ax.tick_params(labelsize=7)
 for ax in list(axs.flat)[n:]:ax.axis('off')
 fig.suptitle(p+' fixed primary copy: wwPDB validation density\nBlue 2mFo−DFc 1σ; red/green mFo−DFc −/+3σ; projected atoms; alternate panels remain distinct',fontsize=11);fig.tight_layout(rect=(0,0,1,.86));fig.savefig(OUT/'density_slices'/f'{p.lower()}_ring_planes.png',dpi=140,bbox_inches='tight',pad_inches=.1);plt.close(fig);print(p,'rendered')
