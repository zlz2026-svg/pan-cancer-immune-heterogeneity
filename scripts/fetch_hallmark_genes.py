# M2 替代增殖度量：下载 Hallmark E2F_TARGETS + G2M_CHECKPOINT，合并去重，剔除 Prolifer 8 基因，映射 ENSG
# 用法：cd 到项目根目录后 python scripts/fetch_hallmark_genes.py
import requests, re, csv

H = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'}

# Prolifer 签名的 8 个基因（主分析调节变量 prolif_abs 用的）
prolifer8 = {"MKI67","PCNA","TOP2A","BIRC5","CCNB1","MCM2","AURKA","CDK1"}

def fetch(gs):
    r = requests.get(f'https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/{gs}.html', headers=H, timeout=40)
    r.raise_for_status()
    genes = re.findall(r'<td valign="top">([A-Z0-9]+)</td>', r.text)
    return set(genes)

e2f = fetch('HALLMARK_E2F_TARGETS')
g2m = fetch('HALLMARK_G2M_CHECKPOINT')
print(f'E2F_TARGETS: {len(e2f)} genes | G2M_CHECKPOINT: {len(g2m)} genes')

merged = e2f | g2m
print(f'合并去重后: {len(merged)} genes')

overlap = merged & prolifer8
alt = merged - prolifer8
print(f'与 Prolifer 8 基因重叠(将被剔除): {sorted(overlap)}')
print(f'剔除后 prolif_alt 基因集: {len(alt)} genes')

# mygene 映射 symbol -> ENSG
import mygene
mg = mygene.MyGeneInfo()
res = mg.querymany(sorted(alt), scopes='symbol', fields='ensembl.gene', species='human',
                    returnall=True, as_dataframe=False)
rows = []
for r in res['out']:
    sym = r.get('query')
    ensg = None
    if 'ensembl' in r and r['ensembl']:
        ens = r['ensembl']
        if isinstance(ens, list):
            ensg = ens[0].get('gene') if ens else None
        else:
            ensg = ens.get('gene')
    if ensg:
        rows.append((sym, ensg))

print(f'成功映射 ENSG: {len(rows)} / {len(alt)}')

with open('data/hallmark_prolif_alt_ensg.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['symbol', 'ensg'])
    w.writerows(rows)

# 记录剔除日志
with open('results/sens/m2_gene_selection_log.txt', 'w') as f:
    f.write(f'E2F_TARGETS n={len(e2f)}, G2M_CHECKPOINT n={len(g2m)}\n')
    f.write(f'merged unique n={len(merged)}\n')
    f.write(f'overlap with Prolifer8 removed: {sorted(overlap)}\n')
    f.write(f'prolif_alt final genes n={len(alt)}, mapped ENSG n={len(rows)}\n')

print('DONE. data/hallmark_prolif_alt_ensg.tsv written.')
