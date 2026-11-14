# Consolidate into single sorted parquet files with published write settings + geom CRS metadata.
import sys, re, os, json, hashlib, pyarrow as pa, pyarrow.parquet as pq, pyarrow.dataset as ds, pyarrow.compute as pc
site=sys.argv[1]; R=os.path.expanduser("~/wetlsp_v1_test2"); T=f"{R}/work/{site}/pq_ds"; O=f"{R}/final/{site}"; stem=site.replace("-","_")
KW=dict(compression="zstd", use_dictionary=True, write_statistics=True)
meta=pq.read_table(f"{T}/meta.parquet"); kv=dict(zip(meta["key"].to_pylist(), meta["value"].to_pylist())); wkt=kv["crs_wkt"]
epsg=re.findall(r'ID\["EPSG",(\d+)\]', wkt)[-1]
md={b"crs_wkt":wkt.encode(), b"crs_epsg":f"EPSG:{epsg}".encode(), b"crs_coordinate_columns":b"x,y",
    b"crs_note":f"Projected CRS of columns x,y (metres); copied verbatim from key crs_wkt in {stem}_pixels_meta.parquet.".encode()}
def w(t, name):
    p=f"{O}/{stem}_{name}.parquet"; tmp=p+".tmp"; pq.write_table(t, tmp, row_group_size=131072, **KW); os.replace(tmp,p); return p
g=ds.dataset(f"{T}/geom").to_table().sort_by("pixel_id"); g=g.replace_schema_metadata(md); w(g,"pixels_geom")
w(meta,"pixels_meta")
# Streaming consolidation (replaces the DuckDB global sort, which ran out of temp disk on 2000 m sites):
# chunk pixel_id ranges are disjoint, so sorting each chunk and writing chunks in pixel_id order gives the
# same global order (pixel_id, series desc, date) with bounded memory.
import glob
p=f"{O}/{stem}_pixels_timeseries.parquet"; tmp=p+".tmp"
SCH=pa.schema([("pixel_id",pa.int32()),("series",pa.string()),("date",pa.date32()),("year",pa.int32()),("evi",pa.float64()),("qc_raw",pa.int8())])
fs=[]
for f in glob.glob(f"{T}/ts/*.parquet"):
    mm=pc.min_max(pq.read_table(f,columns=["pixel_id"])["pixel_id"]); fs.append((mm["min"].as_py(),mm["max"].as_py(),f))
fs.sort(); assert all(fs[k][1]<fs[k+1][0] for k in range(len(fs)-1)), "chunk pixel ranges overlap"
wr=pq.ParquetWriter(tmp, SCH, **KW)
for _,_,f in fs:
    t=pq.read_table(f).select(SCH.names).cast(SCH).sort_by([("pixel_id","ascending"),("series","descending"),("date","ascending")])
    wr.write_table(t, row_group_size=131072); del t
wr.close()
os.replace(tmp,p); n=pq.ParquetFile(p).metadata.num_rows
print(json.dumps({"site":site,"geom_rows":g.num_rows,"ts_rows":n}))
