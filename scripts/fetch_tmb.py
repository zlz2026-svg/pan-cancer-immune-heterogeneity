# -*- coding: utf-8 -*-
# 从 cBioPortal 查询 TCGA PanCancer Atlas 各癌种 mutation count，计算 TMB
import urllib.request, json, csv

cancers = ["ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC",
           "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD",
           "PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]

WES_MB = 38.0  # WES 标准覆盖 38 Mb

results = {}
for ca in cancers:
    study = f"{ca.lower()}_tcga_pan_can_atlas_2018"
    url = (f"https://www.cbioportal.org/api/studies/{study}/clinical-data"
           f"?clinicalDataType=SAMPLE&attributeId=MUTATION_COUNT&projection=DETAILED")
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            data = json.loads(r.read().decode())
        counts = [int(d["value"]) for d in data if d.get("value") and d["value"].isdigit()]
        if counts:
            median_tmb = sorted(counts)[len(counts)//2] / WES_MB
            results[ca] = {"n": len(counts), "median_mut_count": sorted(counts)[len(counts)//2], "median_TMB": round(median_tmb, 2)}
            print(f"{ca}: n={len(counts)}, median TMB={median_tmb:.2f}")
        else:
            print(f"{ca}: 无数据")
    except Exception as e:
        print(f"{ca}: 失败 {e}")

with open("results/cancer_tmb.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["cancer","n","median_mut_count","median_TMB"])
    for ca, v in results.items():
        w.writerow([ca, v["n"], v["median_mut_count"], v["median_TMB"]])
print("DONE, 保存", len(results), "个癌种")
