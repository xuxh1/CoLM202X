"""Rice budget per year from the global tracer history: emission, production,
oxidation, E/P, the three transport channels, paddy area, and emission by
region (South, Southeast and East Asia). Q-15."""
import os, sys; os.environ['HDF5_USE_FILE_LOCKING']='FALSE'
sys.path.insert(0,'/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2/scripts')
import numpy as np, xarray as xr
import v2budget as V
B=V.B
cases=[('paper_v2/v260924c/g2_v2',[1998,1999,2000,2001]),
       ('paper_v2/v260925e/g2_v2e',[2002]),('paper_v2/v260925j/g2_v2f',[2003,2004]),
       ('paper_v2/v260925t/g2_v2g',[2005]),
       ('paper_v2/v260925v/g2_v2i',[2006]),
       ('paper_v2/v260925w/g2_v2j',[2006]),
       ('v260910c/g2_paper_fp2',[2010,2011,2019])]
# optional: rice_yearly.py <case under cases/> <y0> <y1> replaces the list above
if len(sys.argv) > 3:
    cases=[(sys.argv[1], list(range(int(sys.argv[2]), int(sys.argv[3])+1)))]
print('case year | E_rice P_rice O_rice E/P | aere ebul diff | paddyArea_Mkm2 | E by region: SAsia SEAsia EAsia other')
for spec,years in cases:
    ver,name=spec.rsplit('/',1); c=B.Case(ver,name)
    ay=years[0]
    A_land,A_act,A_lake=V.load_areas(c,ay)
    for y in years:
        if not os.path.exists(f'{c.dir}/history/{c.name}_hist_tracer_{y}.nc'): continue
        tr=xr.open_dataset(f'{c.dir}/history/{c.name}_hist_tracer_{y}.nc')
        w=xr.DataArray(B.SEC_PER_MONTH[:tr.sizes['time']],dims='time')
        tg=lambda v,a: float((v*a*w).sum())*B.M_CH4
        E=tg(tr['f_methane_surf_flux_rice'],A_land)
        P=tg(tr['f_methane_prod_tot_rice'],A_act); O=tg(tr['f_methane_oxid_tot_rice'],A_act)
        ae=tg(tr['f_methane_surf_aere_rice'],A_act); eb=tg(tr['f_methane_surf_ebul_rice'],A_act); di=tg(tr['f_methane_surf_diff_rice'],A_act)
        ar=float((tr['f_methane_area_rice'].mean('time')*A_land).sum())/1e12
        lat=tr['lat']; lon=tr['lon']
        fl=(tr['f_methane_surf_flux_rice']*A_land*w).sum('time')*B.M_CH4
        sa=float(fl.where((lat>5)&(lat<35)&(lon>60)&(lon<92)).sum())
        se=float(fl.where((lat>-11)&(lat<25)&(lon>=92)&(lon<150)).sum())
        ea=float(fl.where((lat>=25)&(lat<50)&(lon>=100)&(lon<150)).sum())
        print(f'{c.name} {y} | {E:6.1f} {P:6.1f} {O:6.1f} {E/P if P else 0:5.2f} | {ae:5.1f} {eb:5.1f} {di:5.1f} | {ar:5.2f} | {sa:5.1f} {se:5.1f} {ea:5.1f} {E-sa-se-ea:5.1f}')
        tr.close()
