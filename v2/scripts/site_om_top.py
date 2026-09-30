import glob, netCDF4 as nc, numpy as np
S = '/share/home/dq076/mode/Methane/data/FLUXNET-CH4/Sitedata'
sites = 'BR-Npw CA-SCB DE-Hte DE-SfN DE-Zrk FI-Lom FI-Si2 FI-Sii FR-LGt JP-BBY MY-MLM NZ-Kop SE-Deg US-Los US-Myb US-Sne US-Srr US-StJ US-Tw1 US-Tw4 US-Tw5 US-Uaf US-WPT BW-Gum BW-Nxr ID-Pag NL-Hor RU-Ch2 RU-Che RU-Cok US-A03 US-A10 US-Atq US-BZB US-BZF US-Beo US-Bes US-DPW US-ICs US-Ivo US-LA2 US-NC4 US-NGB US-NGC US-ORv US-OWC'.split()
for s in sites:
    p = glob.glob(f'{S}/{s}_*_site.nc')
    if not p: print(s, 'no file'); continue
    d = nc.Dataset(p[0])
    om = d['soil_OM_density'][:] if 'soil_OM_density' in d.variables else None
    lat = float(d['latitude'][...]) if 'latitude' in d.variables else np.nan
    print(f"{s:7s} lat {lat:6.1f} OM top {'none' if om is None else np.round(np.ravel(om)[:3],1)}")
