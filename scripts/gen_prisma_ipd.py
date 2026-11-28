#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""生成 PRISMA-IPD 清单（npj 要求 meta 提交）——单一数据源 TCGA 二次分析，不适用条目标注 N/A。"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
from docx import Document
from docx.shared import Pt
from docx.oxml.ns import qn

FP = r'E:\1A-AI work\泛癌免疫异质性归因\投稿材料包\npjPrecisionOncology\06-PRISMA清单\【1128】PRISMA-IPD_Checklist.docx'

d = Document()
style = d.styles['Normal']
style.font.name = 'Times New Roman'
style.font.size = Pt(10)
style._element.rPr.rFonts.set(qn('w:eastAsia'), '宋体')

# 标题
h = d.add_paragraph()
r = h.add_run('PRISMA-IPD Checklist')
r.bold = True
r.font.size = Pt(14)

sub = d.add_paragraph()
sub.add_run('Cross-cancer heterogeneity in the prognostic effects of immune microenvironment signatures: a pan-cancer random-effects meta-analysis with meta-regression attribution')

note = d.add_paragraph()
nr = note.add_run('Note: this study is an individual-data meta-analysis and meta-regression of a single public dataset (TCGA PanCancer Atlas), not a systematic review of multiple studies. Items concerning literature search, study selection, and risk-of-bias across studies are therefore marked N/A (not applicable).')
nr.italic = True
nr.font.size = Pt(9)

# 表格
rows = [
    ("Title", "1", "Identify the report as a meta-analysis", "Title page"),
    ("Abstract", "2", "Structured summary (background, methods, results, conclusions)", "Abstract"),
    ("Introduction: Rationale", "3", "Rationale in the context of what is known", "Introduction"),
    ("Introduction: Objectives", "4", "Explicit statement of the questions addressed", "Introduction (last paragraph)"),
    ("Methods: Protocol and registration", "5", "Protocol existence / registration number", "Methods; OSF preregistration https://osf.io/39htc"),
    ("Methods: Eligibility criteria", "6", "Inclusion and exclusion criteria", "Methods (Study design and data)"),
    ("Methods: Information sources", "7", "Description of all information sources and search date", "N/A \u2014 single public dataset (TCGA), no literature search"),
    ("Methods: Search", "8", "Full electronic search strategy for at least one database", "N/A"),
    ("Methods: Study selection", "9", "Process for selecting studies", "N/A \u2014 no multi-study selection"),
    ("Methods: Data collection", "10", "How IPD were requested, collected and managed", "Methods; data downloaded from UCSC Xena Toil hub"),
    ("Methods: Data items", "11", "List and define all variables sought", "Methods (expression, survival, age/sex/stage/grade, moderators)"),
    ("Methods: IPD integrity", "A1", "Aspects of IPD subject to data checking", "Methods (sample matching to \u221201 primary tumors; QC filters)"),
    ("Methods: Risk of bias in individual studies", "12", "Methods of risk-of-bias assessment", "N/A \u2014 single dataset, no per-study bias assessment"),
    ("Methods: Outcomes and effect measures", "13", "All outcomes and principal effect measures", "Methods (OS/DSS; log hazard ratio)"),
    ("Methods: Synthesis", "14", "Meta-analysis and meta-regression methods", "Methods (random-effects REML, Knapp\u2013Hartung, I\u00b2, permutation)"),
    ("Methods: Risk of bias across studies", "15", "Assessment of bias affecting cumulative evidence", "N/A"),
    ("Methods: Additional analyses", "16", "Sensitivity / subgroup / meta-regression analyses", "Methods (Sensitivity analyses)"),
    ("Results: Study selection", "17", "Numbers of studies screened/included, with flow diagram", "N/A \u2014 no literature screening; cohort flow in Fig. 1a"),
    ("Results: Study characteristics", "18", "Characteristics of each included study", "Table S1 (per-cancer cohort characteristics)"),
    ("Results: Risk of bias within studies", "19", "Risk of bias per study", "N/A"),
    ("Results: Results of individual studies", "20", "Effect estimates per study", "Fig. 2 (per-cancer hazard ratios)"),
    ("Results: Synthesis of results", "21", "Results of each meta-analysis with I\u00b2", "Results; Table 2; Fig. 2"),
    ("Results: Risk of bias across studies", "22", "Results of assessment across studies", "N/A"),
    ("Results: Additional analysis", "23", "Results of additional analyses", "Results; Table S3"),
    ("Discussion: Summary of evidence", "24", "Main findings and strength of evidence", "Discussion"),
    ("Discussion: Limitations", "25", "Limitations at study and review level", "Discussion (Limitations)"),
    ("Discussion: Conclusions", "26", "General interpretation and implications", "Discussion"),
    ("Funding", "27", "Sources of funding", "Acknowledgements [to be completed]"),
    ("IPD integrity reporting", "A3", "Important issues found in checking IPD, or none", "N/A \u2014 no issues found; see Methods"),
]

table = d.add_table(rows=1, cols=4)
table.style = 'Table Grid'
hdr = table.rows[0].cells
for i, txt in enumerate(["Section / topic", "Item #", "Checklist item", "Reported on page / section"]):
    hdr[i].paragraphs[0].add_run(txt).bold = True

for sec, num, item, loc in rows:
    cells = table.add_row().cells
    cells[0].text = sec
    cells[1].text = num
    cells[2].text = item
    cells[3].text = loc

d.save(FP)
print('已保存:', FP)

# 复检
d2 = Document(FP)
print('表格行数(含表头):', len(d2.tables[0].rows))
na = sum(1 for row in d2.tables[0].rows if 'N/A' in row.cells[3].text)
print('标注 N/A 的条目:', na)
