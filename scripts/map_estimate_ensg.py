# -*- coding: utf-8 -*-
# ESTIMATE 签名基因 symbol -> ENSG 映射
import mygene, csv

mg = mygene.MyGeneInfo()

genes = []
with open("data/estimate_genes.tsv") as f:
    next(f)  # header
    for line in f:
        g, s = line.strip().split("\t")
        genes.append((g, s))

symbols = [g for g, _ in genes]
out = mg.querymany(symbols, scopes="symbol", fields="ensembl.gene,symbol", species="human", verbose=False)

mapping = {}
for r in out:
    sym = r.get("query")
    ens = r.get("ensembl")
    if not sym or not ens:
        print("  未命中:", sym)
        continue
    if isinstance(ens, list):
        ens = ens[0]
    gene_id = ens.get("gene") if isinstance(ens, dict) else ens
    mapping[sym] = gene_id

print("命中:", len(mapping), "/", len(symbols))

with open("data/estimate_genes_ensg.tsv", "w") as f:
    f.write("gene\tset\tensg\n")
    for g, s in genes:
        if g in mapping:
            f.write(f"{g}\t{s}\t{mapping[g]}\n")
print("DONE")
