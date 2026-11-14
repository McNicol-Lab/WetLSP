import json,sys,os
site=sys.argv[1]; R=os.path.expanduser("~/wetlsp_v1_test2"); C=R+"/code"; W=f"{R}/work/{site}"
p=json.load(open(C+"/pipeline/wetlsp-parameters.json"))
s=p["setup"]
s.update(workDir=W, geojsonDir=W+"/geojson", dataDir=R+"/raw_flat", outDir=W+"/output",
         rScripts=C+"/pipeline/code", rFunctions=C+"/pipeline/wetlsp-functions.R",
         productTable=C+"/metadata/wetlsp-layers.csv", writeCSV=False)
_f=R+"/power/numcores02"
if os.path.exists(_f): s["numCores"]=int(open(_f).read().split()[0])  # harness: day/night cores for step 02 (compute only)
assert s["pipelineVersion"]=="1.0" and p["phenology_parameters"]["vegetation_index"]=="evi"
os.makedirs(W+"/pipeline",exist_ok=True); os.makedirs(W+"/geojson/AMFLX",exist_ok=True)
json.dump(p,open(W+"/pipeline/wetlsp-parameters.json","w"),indent=2)
g=W+f"/geojson/AMFLX/{site}.geojson"
if not os.path.exists(g): os.symlink(C+f"/geojson/AMFLX/{site}.geojson",g)
print("params ok",site)
