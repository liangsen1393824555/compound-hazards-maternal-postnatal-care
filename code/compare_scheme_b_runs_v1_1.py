"""Exact and numeric comparison of the two sealed Scheme-B runs."""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path

def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest().upper()
def leaves(x,p=''):
 out={}
 if isinstance(x,dict):
  for k,v in x.items():out.update(leaves(v,f'{p}.{k}' if p else k))
 elif isinstance(x,list):
  for i,v in enumerate(x):out.update(leaves(v,f'{p}[{i}]'))
 elif isinstance(x,(int,float)) and not isinstance(x,bool):out[p]=float(x)
 return out
def csvnums(p):
 out={}
 with p.open(encoding='utf-8-sig',newline='') as f:
  for i,row in enumerate(csv.DictReader(f)):
   for k,v in row.items():
    try:out[f'{i}.{k}']=float(v)
    except (ValueError,TypeError):pass
 return out
def diff(a,b):
 return None if set(a)!=set(b) else max((abs(a[k]-b[k]) for k in a),default=0.0)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--run-a',type=Path,required=True);ap.add_argument('--run-b',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
 files=['data/ucdp_scheme_b_nigeria.csv','data/ucdp_scheme_b_burkina_faso.csv','data/restricted_scheme_b_nigeria.csv','data/restricted_scheme_b_burkina_faso.csv','data/scheme_b_derivation_manifest.json','sealed/scheme_b_sealed_results.json','sealed/scheme_b_coefficients.csv','sealed/scheme_b_covariances.csv','sealed/scheme_b_estimates.csv','sealed/R_sessionInfo.txt']
 hashes={n:{'sha256':digest(a.run_a/n),'exact_match':digest(a.run_a/n)==digest(a.run_b/n)} for n in files}
 ja=json.loads((a.run_a/'sealed/scheme_b_sealed_results.json').read_text(encoding='utf-8'));jb=json.loads((a.run_b/'sealed/scheme_b_sealed_results.json').read_text(encoding='utf-8'))
 ds={n:diff(csvnums(a.run_a/n),csvnums(a.run_b/n)) for n in ['sealed/scheme_b_coefficients.csv','sealed/scheme_b_covariances.csv','sealed/scheme_b_estimates.csv']}
 num=diff(leaves(ja),leaves(jb));ok=all(x['exact_match'] for x in hashes.values()) and num==0 and all(v==0 for v in ds.values())
 payload={'scope':['derived Scheme-B country event files','restricted samples','support and overlap gates','coefficients','covariances','standardized risks','contrasts','confidence intervals','supportive unadjusted P values','tables','session information'],'criterion':'Exact SHA-256 equality and zero numeric difference','artifact_hashes':hashes,'numeric_max_absolute_difference':num,'table_max_absolute_differences':ds,'pass':ok}
 a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
