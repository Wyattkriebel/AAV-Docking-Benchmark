#!/usr/bin/env python3
"""Render source-only ring-plane 2mFo-DFc 1sigma/difference +/-3sigma contour QA.
Atoms are projected onto each ring plane; contours are a slice, not a projection.
Full three-dimensional atom-centred sphere sampling is in source_diagnostics.json.
"""
import pathlib,json,math
import gemmi,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=pathlib.Path('data/archive/atomic_reference');(OUT/'density_slices').mkdir(exist_ok=True)
def grid(p,k):
 b=gemmi.as_refln_blocks(gemmi.cif.read_file(str(OUT/f'{p}_validation_{k}_map_coef.cif.gz')))[0]
 f,ph=('pdbx_FWT','pdbx_PHWT') if k=='2fo-fc' else ('pdbx_DELFWT','pdbx_DELPHWT')
 g=b.transform_f_phi_to_map(f,ph,sample_rate=3);a=np.asarray(g);a[:]=(a-a.mean())/a.std();return g
for c in json.load(open(OUT/'source_diagnostics.json'))['cases']:
 p=c['case_id'].lower();g=grid(p,'2fo-fc');dg=grid(p,'fo-fc');n=len(c['residues']);cols=min(3,n);rows=math.ceil(n/cols)
 fig,axs=plt.subplots(rows,cols,figsize=(4.4*cols,4.2*rows),squeeze=False)
 for ax,r in zip(axs.flat,c['residues']):
  atoms=r['atom_map'];x=np.array([a['xyz_A'] for a in atoms]);names=[a['atom_name'] for a in atoms];by=dict(zip(names,x));ring=['C2','C3','C4','C5','C6','O6'] if r['prepared_resname']=='SIA' else ['C1','C2','C3','C4','C5','O5'];rx=np.array([by[a] for a in ring]);center=rx.mean(0);u,s,v=np.linalg.svd(rx-center);basis=v[:2];pr=(x-center)@basis.T;lim=max(3.5,np.max(np.abs(pr))+.8);a=np.arange(-lim,lim+.01,.08);xx,yy=np.meshgrid(a,a);z=np.empty_like(xx);dz=np.empty_like(xx)
  for i in range(len(a)):
   for j in range(len(a)):
    pos=gemmi.Position(*(center+xx[i,j]*basis[0]+yy[i,j]*basis[1]));z[i,j]=g.interpolate_value(pos);dz[i,j]=dg.interpolate_value(pos)
  ax.imshow(z,extent=(-lim,lim,-lim,lim),origin='lower',cmap='Blues',vmin=0,vmax=3,alpha=.42)
  ax.contour(xx,yy,z,levels=[1],colors=['#1762ab'],linewidths=.8)
  if dz.min()<=-3:ax.contour(xx,yy,dz,levels=[-3],colors=['#e04a32'],linewidths=1)
  if dz.max()>=3:ax.contour(xx,yy,dz,levels=[3],colors=['#279b4b'],linewidths=1)
  for n1,n2,_ in r['CCD_internal_heavy_bonds']:
   ij=[names.index(n1),names.index(n2)];ax.plot(pr[ij,0],pr[ij,1],color='#555',lw=1,zorder=2)
  ax.scatter(pr[:,0],pr[:,1],s=15,c=[{'C':'#555','O':'#d13e2e','N':'#2c55b8'}[a['element']] for a in atoms],zorder=3)
  for i,name in enumerate(names):ax.annotate(name,pr[i],fontsize=6,xytext=(2,2),textcoords='offset points')
  ds=r['local_density'];ax.set_title(f"{r['deposited']['author_chain']}{r['deposited']['author_resnum']} {r['prepared_resname']} | RSCC {r['validation']['rscc']:.3f}\n3D supported atoms {ds['supported_atoms']}/{ds['total_atoms']}; min diff {ds['min_atom_difference_sigma']:.2f}σ",fontsize=9)
  ax.set_aspect('equal');ax.set_xlabel('ring-plane coordinate (Å)',fontsize=8);ax.set_ylabel('Å',fontsize=8);ax.tick_params(labelsize=7)
 for ax in list(axs.flat)[n:]:ax.axis('off')
 fig.suptitle(c['case_id']+' selected released copy: wwPDB validation density\nBlue 2mFo−DFc 1σ; red/green mFo−DFc −/+3σ; atoms projected onto ring-plane slice',fontsize=11)
 fig.tight_layout(rect=(0,0,1,.92));fig.savefig(OUT/'density_slices'/f'{p}_ring_planes.png',dpi=140);plt.close(fig)
 print(p,'rendered')
