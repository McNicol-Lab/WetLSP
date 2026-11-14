import sys, os, json, duckdb, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":15,"axes.titlesize":16,"axes.labelsize":15,"legend.fontsize":12,"xtick.labelsize":12,"ytick.labelsize":13,"figure.titlesize":20})
import matplotlib.dates as mdates
from netCDF4 import Dataset
site=sys.argv[1]; R=os.path.expanduser("~/wetlsp_v1_test2"); stem=site.replace("-","_")
A=os.environ.get("CMP_A",f"{R}/ref/{site}"); B=f"{R}/final/{site}"; O=os.environ.get("CMP_O",f"{R}/compare/{site}"); os.makedirs(O,exist_ok=True)
c=duckdb.connect(); c.sql("set memory_limit='8GB'"); c.sql(f"set temp_directory='{R}/work/{site}/duck_tmp'")
c.sql(f"create view a as select * from '{A}/{stem}_pixels_timeseries.parquet'")
c.sql(f"create view b as select * from '{B}/{stem}_pixels_timeseries.parquet'")
S={"site":site}
q=lambda s: c.sql(s).df()
S["pixels_v01"]=int(c.sql("select count(distinct pixel_id) from a").fetchone()[0]); S["pixels_v10"]=int(c.sql("select count(distinct pixel_id) from b").fetchone()[0])
raw="count(*) AS n_rows, count(evi) AS nonNA, coalesce(sum((evi=0)::int),0) AS zeros, coalesce(sum((evi<=0)::int),0) AS le0, coalesce(sum((evi<-0.2 or evi>1)::int),0) AS out_of_range"
S["raw_v01"]=q(f"select {raw} from a where series='raw'").iloc[0].astype(int).to_dict()
S["raw_v10"]=q(f"select {raw} from b where series='raw'").iloc[0].astype(int).to_dict()
qc=q("select coalesce(qc_raw,-1) code, count(*) n from b where series='raw' group by 1 order by 1"); tot=qc.n.sum()
S["qc_raw_counts"]={int(k):int(v) for k,v in zip(qc.code,qc.n)}; S["qc_raw_pct"]={int(k):round(100*v/tot,3) for k,v in zip(qc.code,qc.n)}
S["pct_raw_flagged_nonzero_code"]=round(100*qc[qc.code>0].n.sum()/tot,3)
nv=qc[qc.code.isin([0,2,3,4,5])].n.sum(); S["n_raw_valid_input"]=int(nv); S["pct_removed_codes2to4_of_valid"]=round(100*qc[qc.code.between(2,4)].n.sum()/nv,3)
S["pct_despike2_of_valid"]=round(100*qc[qc.code==2].n.sum()/nv,3); S["pct_floored5_of_valid"]=round(100*qc[qc.code==5].n.sum()/nv,3)
S["code1_where_v01_had_value"]=int(c.sql("select count(*) from a join b using(pixel_id,series,date) where series='raw' and b.qc_raw=1 and a.evi is not null").fetchone()[0])
S["pct_removed_codes1to4"]=round(100*qc[qc.code.between(1,4)].n.sum()/tot,3); S["pct_despike_code2"]=round(100*qc[qc.code==2].n.sum()/tot,3); S["pct_floored_code5"]=round(100*qc[qc.code==5].n.sum()/tot,3)
S["v10_code_on_v01_zero_obs"]={int(k):int(v) for k,v in q("select coalesce(b.qc_raw,-1) code, count(*) n from a join b using(pixel_id,series,date) where series='raw' and a.evi=0 group by 1 order by 1").values}
# spline metrics per pixel-year
m=lambda t: f"select pixel_id, year, max(evi) peak, median(evi) filter (where month(date) between 6 and 8) mid, count(evi) n from {t} where series='spline' and year between 2021 and 2024 group by 1,2"
c.sql(f"create table ma as {m('a')}"); c.sql(f"create table mb as {m('b')}")
J=q("select year, a.peak pa, b.peak pb, a.mid ma_, b.mid mb_, a.n na, b.n nb from ma a join mb b using(pixel_id,year)")
S["spline"]={}
for y,g in J.groupby("year"):
    dp=(g.pb-g.pa).dropna(); dm=(g.mb_-g.ma_).dropna()
    S["spline"][int(y)]={"n_pix":len(g),"peak_v01_med":float(g.pa.median()),"peak_v10_med":float(g.pb.median()),"dpeak_med":float(dp.median()),"dpeak_p10":float(dp.quantile(.1)),"dpeak_p90":float(dp.quantile(.9)),
        "mid_v01_med":float(g.ma_.median()),"mid_v10_med":float(g.mb_.median()),"dmid_med":float(dm.median()),"dmid_p10":float(dm.quantile(.1)),"dmid_p90":float(dm.quantile(.9))}
ca=q("select \"year\" AS yr, count(*) AS pix_v01, avg(n) AS days_v01 from ma group by 1"); cb=q("select \"year\" AS yr, count(*) AS pix_v10, avg(n) AS days_v10 from mb group by 1")
cov=pd.DataFrame({"yr":[2021,2022,2023,2024]}).merge(ca,how="left").merge(cb,how="left").fillna(0).rename(columns={"yr":"year"})
S["coverage"]=cov.to_dict(orient="records")
# phenometrics from nc
V=["OGI","50PCGI","Peak","50PCGD","OGMn","EVImax","EVIamp","NumCycles"]; P={}
def rd(p):
    d=Dataset(p); d.set_auto_mask(True); out={v:np.ma.filled(d[v][:].astype(float),np.nan) for v in V if v in d.variables}
    out["_shape"]=d["OGI"].shape; out["_pv"]=getattr(d,"product_version",""); d.close(); return out
for y in range(2021,2025):
    fa=f"{A}/{site}-wetlsp-{y}.nc"; fb=f"{B}/{site}-wetlsp-{y}.nc"
    if not (os.path.exists(fa) and os.path.exists(fb)): P[y]="missing "+("v01 " if not os.path.exists(fa) else "")+("v10" if not os.path.exists(fb) else ""); continue
    a,b=rd(fa),rd(fb); r={"shape_v01":a["_shape"],"shape_v10":b["_shape"],"pv":[a["_pv"],b["_pv"]]}
    for v in V:
        r[v]={"v01_med":float(np.nanmedian(a[v])) if np.isfinite(a[v]).any() else None,"v10_med":float(np.nanmedian(b[v])) if np.isfinite(b[v]).any() else None,
              "valid_v01":round(float(np.isfinite(a[v]).mean()),4),"valid_v10":round(float(np.isfinite(b[v]).mean()),4)}
        if a["_shape"]==b["_shape"]:
            dd=(b[v]-a[v]); dd=dd[np.isfinite(dd)]; r[v]["pix_diff_med"]=float(np.median(dd)) if dd.size else None; r[v]["pix_diff_p90abs"]=float(np.quantile(abs(dd),.9)) if dd.size else None
    P[y]=r
S["phenometrics"]=P
warn=[]
if S["pct_removed_codes2to4_of_valid"]>5: warn.append(f"v1.0 QC removed {S['pct_removed_codes2to4_of_valid']}% of valid raw obs (codes 2-4, >5%)")
if S["pct_despike2_of_valid"]>5: warn.append(f"despike alone removed {S['pct_despike2_of_valid']}% of valid obs (>5%)")
if S["raw_v10"]["zeros"]>0: warn.append("zeros remain in v1.0 raw")
if S["raw_v10"]["out_of_range"]>0: warn.append("out-of-range EVI remains in v1.0 raw")
for rec in S["coverage"]:
    if rec["pix_v10"] < 0.9*max(rec["pix_v01"],1): warn.append(f"{rec['year']}: v1.0 spline pixels {rec['pix_v10']} < 90% of v0.1 {rec['pix_v01']}")
for y,r in P.items():
    if isinstance(r,str): warn.append(f"{y} nc {r}")
    elif r["OGI"]["valid_v10"] < 0.5*max(r["OGI"]["valid_v01"],1e-9): warn.append(f"{y}: nc OGI valid fraction {r['OGI']['valid_v10']} vs {r['OGI']['valid_v01']}")
S["warnings"]=warn
json.dump(S,open(f"{O}/summary.json","w"),indent=1,default=str)
# --- plot 1: example pixels (most changed midseason in 2022/23) raw+spline
yr=max(S["spline"], key=lambda y:S["spline"][y]["n_pix"]) if S["spline"] else 2022
nv=q(f"select pixel_id, count(*) filter (where qc_raw=0 and evi is not null) nv from b where series='raw' and year={yr} group by 1")
lo,hi=nv.nv.quantile([.4,.6]); c.register("nvok",nv[(nv.nv>=lo)&(nv.nv<=hi)])
cand=q(f"select pixel_id, abs(b.mid-a.mid) d from ma a join mb b using(pixel_id,year) join nvok using(pixel_id) where year={yr} and a.mid is not null and b.mid is not null")
dlo,dhi=cand.d.quantile([.4,.6]); pool=cand[(cand.d>=dlo)&(cand.d<=dhi)]
ex=pool.sample(n=min(6,len(pool)),random_state=42)[["pixel_id"]].iloc[:0]
med=pool.sample(n=min(6,len(pool)),random_state=42)[["pixel_id"]]
pids=list(ex.pixel_id)+list(med.pixel_id); c.register("pp",pd.DataFrame({"pixel_id":pids}))
da=q(f"select * from a join pp using(pixel_id) where year={yr}"); db=q(f"select * from b join pp using(pixel_id) where year={yr}")
fig,axs=plt.subplots(2,3,figsize=(16,9),sharey=True); axs=axs.ravel()
for ax,pid in zip(axs,pids):
    ra=da[(da.pixel_id==pid)&(da.series=="raw")]; rb=db[(db.pixel_id==pid)&(db.series=="raw")]
    ax.scatter(ra.date,ra.evi,s=14,c="0.6",label="v0.1 (published) raw")
    bad=rb[rb.qc_raw.fillna(0).between(1,4)]; ax.scatter(rb.date,rb.evi,s=6,c="k",label="v1.0 (patched) raw, kept")
    ax.scatter(bad.date,[-0.18]*len(bad),marker="x",c="r",s=14,label="v1.0 removed (qc 1-4)")
    for t,col,lab in [(da,"tab:blue","v0.1 (published) spline"),(db,"tab:orange","v1.0 (patched) spline")]:
        s=t[(t.pixel_id==pid)&(t.series=="spline")].sort_values("date"); ax.plot(s.date,s.evi,c=col,lw=1.5,label=lab)
    ax.set_title(f"pixel {pid} ({'typical'})"); ax.set_ylim(-0.22,max(1.0,float(np.nanmax(pd.concat([da[da.pixel_id==pid].evi,db[db.pixel_id==pid].evi])))+0.05)); ax.xaxis.set_major_formatter(mdates.DateFormatter('%b')); ax.grid(alpha=.3)
axs[0].set_ylabel("EVI"); axs[3].set_ylabel("EVI"); axs[0].legend(loc="upper left",fontsize=10); fig.suptitle(f"{site} {yr}: representative pixels (median obs count and change), WetLSP v0.1 (published) vs v1.0"); fig.tight_layout(); fig.savefig(f"{O}/{site}_pixels_{yr}.png",dpi=130); plt.close(fig)
# --- plot 2: distributions
fig,axs=plt.subplots(1,3,figsize=(16,9))
axs[0].hist((J.pb-J.pa).dropna().clip(-.5,.5),bins=100,color="tab:purple"); axs[0].set_title("Peak EVI change\n(v1.0 − v0.1, per pixel-year)"); axs[0].set_xlabel("Δ EVI"); axs[0].axvline(0,c="k")
axs[1].hist((J.mb_-J.ma_).dropna().clip(-.5,.5),bins=100,color="tab:green"); axs[1].set_title("Jun–Aug median EVI change\n(v1.0 − v0.1)"); axs[1].set_xlabel("Δ EVI"); axs[1].axvline(0,c="k")
lab={-1:"NA",0:"0 ok",1:"1 invalid",2:"2 despike",3:"3 2nd-spline",4:"4 in-season<floor",5:"5 floored"}
axs[2].bar([lab.get(k,str(k)) for k in S["qc_raw_pct"]],list(S["qc_raw_pct"].values()),color="tab:red"); axs[2].set_yscale("log"); axs[2].set_title("v1.0 QC code\n(% of raw observations)"); plt.setp(axs[2].get_xticklabels(),rotation=40,ha="right")
fig.suptitle(f"{site}: spline changes and QC flags, WetLSP v0.1 vs v1.0 (2021–2024)"); fig.tight_layout(); fig.savefig(f"{O}/{site}_distributions.png",dpi=130); plt.close(fig)
# --- plot 3: phenometric medians by year + valid fraction
ys0=[y for y in P if not isinstance(P[y],str)]
if ys0:
  fig,axs=plt.subplots(1,2,figsize=(16,9)); ys=[y for y in P if not isinstance(P[y],str)]
  for v,col in zip(["OGI","50PCGI","Peak","50PCGD","OGMn"],plt.cm.tab10.colors):
      axs[0].plot(ys,[P[y][v]["v01_med"] for y in ys],"o--",c=col,alpha=.6); axs[0].plot(ys,[P[y][v]["v10_med"] for y in ys],"s-",c=col,label=v)
  axs[0].set_title("Median phenology date\n(dashed = v0.1 published, solid = v1.0 patched)"); axs[0].set_ylabel("day of year"); axs[0].legend(); axs[0].set_xticks(ys)
  w=.35; x=np.arange(len(ys))
  axs[1].bar(x-w/2,[P[y]["OGI"]["valid_v01"] for y in ys],w,label="v0.1 (published)"); axs[1].bar(x+w/2,[P[y]["OGI"]["valid_v10"] for y in ys],w,label="v1.0 (patched)")
  axs[1].set_xticks(x,ys); axs[1].set_title("Share of map with a valid\ngreen-up date (OGI)"); axs[1].legend()
  fig.suptitle(f"{site}: phenometrics, WetLSP v0.1 vs v1.0"); fig.tight_layout(); fig.savefig(f"{O}/{site}_phenometrics.png",dpi=130); plt.close(fig)
# --- headline slide
r0,r1=S["raw_v01"],S["raw_v10"]; sp=S["spline"].get(yr,{}); ph=P.get(yr,{}) if not isinstance(P.get(yr),str) else {}
fig=plt.figure(figsize=(16,9)); gs=fig.add_gridspec(2,3,width_ratios=[1,1,0.9])
for k,pid in enumerate(pids[:4]):
    ax=fig.add_subplot(gs[k//2,k%2])
    ra=da[(da.pixel_id==pid)&(da.series=="raw")]; rb=db[(db.pixel_id==pid)&(db.series=="raw")]
    ax.scatter(ra.date,ra.evi,s=30,facecolors="none",edgecolors="tab:blue",label="v0.1 (published) raw")
    ax.scatter(rb.date,rb.evi,s=12,c="tab:orange",label="v1.0 (patched) raw, kept")
    rm_=rb[rb.qc_raw.fillna(0).between(2,4)].merge(ra[["date","evi"]],on="date",suffixes=("","_a"))
    ax.scatter(rm_.date,rm_.evi_a,marker="x",s=70,c="red",lw=2,label="removed by v1.0 QC")
    for t,col,lab in [(da,"tab:blue","v0.1 (published) spline"),(db,"tab:orange","v1.0 (patched) spline")]:
        s_=t[(t.pixel_id==pid)&(t.series=="spline")].sort_values("date"); ax.plot(s_.date,s_.evi,c=col,lw=2.5,label=lab)
    ax.set_ylim(-0.22,max(1.0,float(np.nanmax(pd.concat([da[da.pixel_id==pid].evi,db[db.pixel_id==pid].evi])))+0.05)); ax.grid(alpha=.3); ax.xaxis.set_major_formatter(mdates.DateFormatter('%b')); ax.set_title(f"pixel {pid}, {yr}")
    if k==0: ax.legend(loc="upper left",fontsize=11)
    if k%2==0: ax.set_ylabel("EVI")
tx=fig.add_subplot(gs[:,2]); tx.axis("off")
lines=[f"{site}, {yr}","",f"Raw zeros: {r0['zeros']:,} → {r1['zeros']:,}",f"Raw EVI ≤ 0: {r0['le0']:,} → {r1['le0']:,}",
 f"Out of range: {r0['out_of_range']:,} → {r1['out_of_range']:,}",f"Valid obs removed by QC: {S['pct_removed_codes2to4_of_valid']:.1f}%",f"  (despike {S['pct_despike2_of_valid']:.1f}%)",f"Floored off-season: {S['pct_floored5_of_valid']:.1f}%",f"Any QC flag (codes 1-5): {S['pct_raw_flagged_nonzero_code']:.1f}%",""]
if sp: lines+=[f"Peak EVI (median): {sp['peak_v01_med']:.3f} → {sp['peak_v10_med']:.3f}",f"Jun–Aug EVI: {sp['mid_v01_med']:.3f} → {sp['mid_v10_med']:.3f}",f"Per-pixel shift: peak {sp['dpeak_med']:+.3f}, mid {sp['dmid_med']:+.3f}",""]
if ph: lines+=["Median dates (DOY):"]+[f"  {lab}: {ph[v]['v01_med']:.0f} → {ph[v]['v10_med']:.0f}" for v,lab in [("OGI","green-up"),("Peak","peak"),("50PCGD","senescence")] if ph[v]['v01_med'] is not None and ph[v]['v10_med'] is not None]
tx.text(0.02,0.98,"\n".join(lines),va="top",fontsize=17,family="DejaVu Sans",transform=tx.transAxes,bbox=dict(boxstyle="round,pad=0.6",fc="#f4f4f4",ec="0.5"))
fig.suptitle(f"{site}: WetLSP v0.1 (published) vs v1.0 (patched)"); fig.tight_layout(); fig.savefig(f"{O}/{site}_headline.png",dpi=150); plt.close(fig)
print(json.dumps(S,indent=1,default=str))
