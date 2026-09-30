"""Rice CH4 production by irrigated-rice share of the cell, O-0 (2014) against O-1 (2015)."""
import numpy as np, xarray as xr
P='/share/home/dq076/mode/Methane/cases/paper_v2'
cf=xr.open_dataset(f'{P}/v260925w/g2_v2q/landdata/diag/cropfrac_elm_2005.nc')
f=cf.cropfrac_elm.where(cf.cropfrac_elm>-1e30).fillna(0)
r61=f.sel(TypeIndex=47).values; r62=f.sel(TypeIndex=48).values
share=np.where(r61+r62>0, r62/np.maximum(r61+r62,1e-12), np.nan)
out={}
for c,y in (('v260927b/g2_o0',2014),('v260927g/g2_oo1',2015)):
    n=c.split('/')[1]
    t=xr.open_dataset(f'{P}/{c}/history/{n}_hist_tracer_{y}.nc')
    h=xr.open_dataset(f'{P}/{c}/history/{n}_hist_{y}.nc')
    A=h['landarea'].values*1e6 if 'landarea' in h else None
    pr=(t.f_methane_prod_tot_rice.fillna(0)).mean('time').values*16.043*86400*365   # g m-2 land yr-1
    out[n]=pr*A/1e12
bins=[0,0.1,0.3,0.5,0.7,0.9,1.01]
print('irrigated share  cells  O-0 P  O-1 P  ratio')
for a,b in zip(bins[:-1],bins[1:]):
    m=(share>=a)&(share<b)
    p0=np.nansum(out['g2_o0'][m]); p1=np.nansum(out['g2_oo1'][m])
    print(f'{a:.1f}-{b:.1f}  {m.sum():6d} {p0:6.2f} {p1:6.2f} {p1/max(p0,1e-9):5.2f}')
