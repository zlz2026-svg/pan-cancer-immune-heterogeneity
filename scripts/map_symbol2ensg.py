# -*- coding: utf-8 -*-
# 查询签名基因 symbol -> ENSG ID 映射
import mygene, json

mg = mygene.MyGeneInfo()

signatures = {
    "CD8_cyto": ["CD8A","CD8B","CD3D","CD3E","GZMA","GZMB","PRF1","IFNG","PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","TOX"],
    "TLS": ["CXCL13","CCL19","CCL21","CD79A","MS4A1","LTB","TNFRSF17"],
    "CAF": ["COL1A1","COL1A2","COL3A1","ACTA2","FAP","PDGFRB","POSTN","FN1","MMP2"],
    "TAN": ["CSF3R","S100A8","S100A9","FCGR3B","CXCR2","OLR1"],
    "cDC2": ["CD1C","FCER1A","CLEC10A","CD1E"],
    "IFNg_path": ["IFNG","CXCL9","CXCL10","CXCL11","GBP1","STAT1","IRF1","GBP2"],
    "Prolifer": ["MKI67","PCNA","TOP2A","BIRC5","CCNB1","MCM2","AURKA","CDK1"],
}

all_symbols = sorted(set(s for g in signatures.values() for s in g))
print("查询基因数:", len(all_symbols))

out = mg.querymany(all_symbols, scopes="symbol", fields="ensembl.gene,symbol", species="human", verbose=False)

# symbol -> ensg（取第一个）
mapping = {}
for r in out:
    sym = r.get("query")
    ens = r.get("ensembl")
    if not sym or not ens:
        print("  未命中:", sym, "->", r.get("notfound"))
        continue
    if isinstance(ens, list):
        ens = ens[0]
    gene_id = ens.get("gene") if isinstance(ens, dict) else ens
    mapping[sym] = gene_id

print("命中:", len(mapping), "/", len(all_symbols))

# 输出
with open("symbol2ensg.json", "w") as f:
    json.dump(mapping, f, indent=2)

# 打印各签名命中情况
for nm, genes in signatures.items():
    hit = [g for g in genes if g in mapping]
    miss = [g for g in genes if g not in mapping]
    print(f"{nm}: 命中 {len(hit)}/{len(genes)}", ("缺 " + ",".join(miss)) if miss else "")
print("DONE")
