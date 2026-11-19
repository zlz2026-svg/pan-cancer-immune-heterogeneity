# GLASS 胶质瘤独立队列外部验证 — 下载表达 + 临床 + 签名评分
# 数据源：cBioPortal difg_glass_2019 (GLASS Consortium, Nature 2019)
import urllib.request, json, csv

BASE = "https://www.cbioportal.org/api"

def get(url):
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.loads(r.read().decode())

def post(url, body):
    req = urllib.request.Request(url, data=json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read().decode())

# 签名基因集（与主分析一致）
signatures = {
    "CD8_cyto":  ["CD8A","CD8B","CD3D","CD3E","GZMA","GZMB","PRF1","IFNG","PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","TOX"],
    "TLS":       ["CXCL13","CCL19","CCL21","CD79A","MS4A1","LTB","TNFRSF17"],
    "CAF":       ["COL1A1","COL1A2","COL3A1","ACTA2","FAP","PDGFRB","POSTN","FN1","MMP2"],
    "TAN":       ["CSF3R","S100A8","S100A9","FCGR3B","CXCR2","OLR1"],
    "cDC2":      ["CD1C","FCER1A","CLEC10A","CD1E"],
    "IFNg_path": ["IFNG","CXCL9","CXCL10","CXCL11","GBP1","STAT1","IRF1","GBP2"],
    "Prolifer":  ["MKI67","PCNA","TOP2A","BIRC5","CCNB1","MCM2","AURKA","CDK1"],
}
all_syms = sorted(set(g for v in signatures.values() for g in v))
print(f"签名基因数(去重): {len(all_syms)}")

# 1. 样本列表（有表达数据的样本）
samples = get(f"{BASE}/studies/difg_glass_2019/samples?projection=DETAILED")
sample_ids = [s["sampleId"] for s in samples]
patient_of = {s["sampleId"]: s["patientId"] for s in samples}
print(f"样本数: {len(sample_ids)}")

# 2. 患者 OS 临床
clin = get(f"{BASE}/studies/difg_glass_2019/clinical-data?clinicalDataType=PATIENT&projection=DETAILED")
os_status = {}
os_months = {}
for c in clin:
    pid = c["patientId"]
    if c["clinicalAttributeId"] == "OS_STATUS":
        os_status[pid] = c["value"]
    elif c["clinicalAttributeId"] == "OS_MONTHS":
        os_months[pid] = c["value"]
print(f"OS_STATUS 记录: {len(os_status)}, OS_MONTHS 记录: {len(os_months)}")

# 3. symbol -> entrez (mygene)
import mygene
mg = mygene.MyGeneInfo()
res = mg.querymany(all_syms, scopes="symbol", fields="entrezgene", species="human")
sym2entrez = {}
for r in res:
    if "entrezgene" in r:
        sym2entrez[r["query"]] = str(r["entrezgene"])
miss = [s for s in all_syms if s not in sym2entrez]
print(f"symbol->entrez 命中: {len(sym2entrez)}/{len(all_syms)}, 缺失: {miss}")

# 4. fetch 表达数据
entrez_ids = sorted(set(sym2entrez.values()), key=int)
profile = "difg_glass_2019_rna_seq_mrna"
body = {"sampleIds": sample_ids, "entrezGeneIds": [int(x) for x in entrez_ids]}
md = post(f"{BASE}/molecular-profiles/{profile}/molecular-data/fetch", body)
print(f"表达数据条数: {len(md)}")

# 构建 样本×基因 矩阵（symbol 维度）
entrez2sym = {v: k for k, v in sym2entrez.items()}
expr = {}
for row in md:
    sid = row["sampleId"]
    eg = str(row["entrezGeneId"])
    if eg not in entrez2sym:
        continue
    sym = entrez2sym[eg]
    expr.setdefault(sid, {})[sym] = float(row["value"])
print(f"有表达数据的样本数: {len(expr)}")

# 5. 签名评分（z-score 标准化 + 均值，与主分析一致）
import statistics
def sig_score(syms, expr_by_sample):
    out = {}
    for sid, gdict in expr_by_sample.items():
        vals = [gdict[g] for g in syms if g in gdict]
        if len(vals) < 2:
            out[sid] = None
            continue
        # 跨基因 z-score（单样本内不做，主分析是跨样本 z-score；这里简化：先算样本间 z-score 需全矩阵）
        out[sid] = sum(vals) / len(vals)  # 先算均值，z-score 在 R 里做
    return out

# 先输出原始 样本×基因 矩阵 + 临床，签名评分在 R 里做（与主分析完全一致）
rows = []
header = ["sample", "patient", "OS", "OS.time"] + all_syms
with open("results/glass_expr_clin.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(header)
    for sid in sample_ids:
        pid = patient_of[sid]
        if sid not in expr:
            continue
        os_v = os_status.get(pid, "")
        os_t = os_months.get(pid, "")
        if os_v == "1:DECEASED":
            os_bin = 1
        elif os_v == "0:LIVING":
            os_bin = 0
        else:
            os_bin = ""
        row = [sid, pid, os_bin, os_t] + [expr[sid].get(g, "") for g in all_syms]
        w.writerow(row)

n_rows = sum(1 for _ in open("results/glass_expr_clin.csv")) - 1
print(f"输出结果: results/glass_expr_clin.csv, 样本行数 {n_rows}")
