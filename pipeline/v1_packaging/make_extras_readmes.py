import json, os, re, shutil, datetime, sys
V1='/Volumes/home/wlsp/wetlsp-v1'; PUB='/Volumes/home/wlsp/wetlsp-ess-dive-parquet/AT-Nsd'
BK=os.path.expanduser('~/wetlsp_v1_test2/readme_bak_20261114'); os.makedirs(BK,exist_ok=True)
extras=[l.split()[0] for l in open(os.path.expanduser('~/wetlsp_v1_test2/extras_sites.txt')) if l.strip()]+['FI-Lom','FI-Si2','HK-MPM','US-A03']
extras=sorted(set(extras)); now=datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')
md0=open(f'{PUB}/README_parquet.md').read(); js0=json.load(open(f'{PUB}/README_parquet.json'))
V1NOTE=("\n## WetLSP v1.0 notes\n\nThis site was processed with WetLSP pipeline version 1.0 (`product_version` = \"v1.0\"). "
"Differences from the layout above:\n\n"
"- The time series has an extra column `qc_raw` (int8; null for `series = 'spline'`): per-observation QC code for raw EVI "
"(0 ok, 1 invalid (<=0, fill, EVI outside [-0.2, 1]), 2 despike, 3 second-spline QC, 4 in-season below dormant floor, 5 floored off-season). "
"Raw rows set to NA by v1.0 QC are kept with `evi = null` so their `qc_raw` is visible.\n"
"- The meta table has an extra key `product_version` (\"v1.0\").\n"
"- Fetch radius (`radius_m`): the catalog `final_fetch_r` where one exists; otherwise 300 m for extra-tropical sites (|lat| >= 23.436 deg) "
"and 600 m for tropical sites (|lat| < 23.436 deg). The value used is in the meta table key `radius_m`.\n"
"- This site is not in the earlier published WetLSP package, so no comparison figures against it are provided.\n"
"- See `README_v1.md` in this folder for processing provenance.\n")
def fix_v1readme(p):
    t=open(p).read(); o=t
    t=t.replace("For candidate extra sites (`candidate_extras/`) the radius","For extra sites not in the published package the radius")
    t=t.replace('"Fetch radius for candidate extra sites"','"Fetch radius for extra sites"').replace("## Fetch radius for candidate extra sites","## Fetch radius for extra sites")
    t=t.replace("Candidate extra sites that lack","Extra sites that lack")
    if t!=o:
        s=os.path.basename(os.path.dirname(p)); shutil.copy2(p,f'{BK}/{s}_README_v1.md'); open(p,'w').write(t); return 1
    return 0
n=sum(fix_v1readme(f'{V1}/{d}/README_v1.md') for d in sorted(os.listdir(V1)) if os.path.isfile(f'{V1}/{d}/README_v1.md'))
print('README_v1 fixed:',n)
left=[d for d in os.listdir(V1) if os.path.isfile(f'{V1}/{d}/README_v1.md') and 'candidate' in open(f'{V1}/{d}/README_v1.md').read().lower()]
print('still mention candidate:',left)
made=[]
for s in extras:
    d=f'{V1}/{s}'
    if not os.path.isdir(d): print('MISSING',s); continue
    stem=s.replace('-','_')
    need=[f'{stem}_pixels_{k}.parquet' for k in ('geom','meta','timeseries')]
    if not all(os.path.isfile(f'{d}/{f}') for f in need): print('NOT FULL',s); continue
    md=md0.replace('# WetLSP Pixel-Level EVI Time Series (Parquet)','# WetLSP v1.0 Pixel-Level EVI Time Series (Parquet)',1)
    md=md.replace('| evi      | float64 | Enhanced Vegetation Index                                    |',
      '| evi      | float64 | Enhanced Vegetation Index                                    |\n| qc_raw   | int8    | v1.0 QC code for raw EVI (null for spline; see v1.0 notes)   |')
    md=md.replace('| pixel_id_source | Source of `pixel_id` (e.g. `cell`)                           |',
      '| pixel_id_source | Source of `pixel_id` (e.g. `cell`)                           |\n| product_version | WetLSP product version (`v1.0`)                              |')
    md=md.replace('├── README_parquet.json                 # machine-readable version of this file',
      '├── README_parquet.json                 # machine-readable version of this file\n└── README_v1.md                        # v1.0 processing provenance')
    md=md.replace('## Citation & Attribution',V1NOTE+'\n## Citation & Attribution')
    js=json.loads(json.dumps(js0))
    js['title']='WetLSP v1.0 Pixel-Level EVI Time Series (Parquet)'; js['version']='1.0.0'; js['product_version']='v1.0'; js['site_id']=s
    js['generated']={'date_utc':now,'format':'application/json','status':'final','note':'WetLSP v1.0 site README, adapted from the published package README_parquet.json v2.1.0'}
    js['v1_changes']={'timeseries_extra_columns':[{'name':'qc_raw','type':'int8','description':'Per-observation QC code for raw EVI; null for series=spline. 0 ok, 1 invalid (<=0, fill, EVI outside [-0.2, 1]), 2 despike, 3 second-spline QC, 4 in-season below dormant floor, 5 floored off-season. Raw rows set to NA by v1.0 are kept with evi=null.'}],
      'meta_extra_keys':[{'key':'product_version','description':'WetLSP product version (v1.0)'}],
      'fetch_radius_rule':'Catalog final_fetch_r where available; otherwise 300 m extra-tropical (|lat| >= 23.436 deg), 600 m tropical (|lat| < 23.436 deg). Value used recorded in meta key radius_m.',
      'comparison_figures':'Not provided: site is not in the earlier published WetLSP package.',
      'provenance_file':'README_v1.md'}
    open(f'{d}/README_parquet.md','w').write(md); json.dump(js,open(f'{d}/README_parquet.json','w'),indent=2); made.append(s)
print('README_parquet written:',len(made),made)
