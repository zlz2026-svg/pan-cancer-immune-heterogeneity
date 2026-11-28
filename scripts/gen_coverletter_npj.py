#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""生成 npj Precision Oncology 版 Cover Letter（参考王丹 BBA 五段式，改 npj 语境）
作者团队占位符：on behalf of the authorship team 后的名单待老板填。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
from docx import Document
from docx.shared import Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn

FP = r'E:\1A-AI work\泛癌免疫异质性归因\投稿材料包\npjPrecisionOncology\01-CoverLetter\【1128】CoverLetter_npjPrecisionOncology_V1.docx'

d = Document()
style = d.styles['Normal']
style.font.name = 'Times New Roman'
style.font.size = Pt(11)
style._element.rPr.rFonts.set(qn('w:eastAsia'), '宋体')

def para(text, bold=False, align=None, space_after=6):
    p = d.add_paragraph()
    r = p.add_run(text)
    r.bold = bold
    if align is not None:
        p.alignment = align
    p.paragraph_format.space_after = Pt(space_after)
    return p

para('November 28, 2026', space_after=12)
para('The Editors')
para('npj Precision Oncology')
para('Nature Portfolio', space_after=12)
para('Dear Editor,', space_after=12)

para('We are pleased to submit our original research article, entitled \u201cCross-cancer heterogeneity in the prognostic effects of immune microenvironment signatures: a pan-cancer random-effects meta-analysis with meta-regression attribution\u201d, for consideration as an Article in npj Precision Oncology.')

para('Relevance to the journal.', bold=True, space_after=2)
para('npj Precision Oncology publishes precision-oncology research that dissects tumor heterogeneity, the tumor microenvironment and molecular biomarkers with translational implications. Our study addresses a long-standing contradiction in the field \u2014 why the same immune microenvironment signatures point in opposite prognostic directions across cancers \u2014 and resolves it quantitatively. Rather than introducing another biomarker, we attribute the between-cancer heterogeneity (\u03c4\u00b2) of seven immune signatures to tumor proliferation, non-tumoral content, germ-layer origin and mutational burden by meta-regression, an explanation-first framing that aligns with the journal\u2019s emphasis on mechanism-oriented precision oncology.')

para('What is new.', bold=True, space_after=2)
para('Prior pan-cancer studies have documented that immune signatures are heterogeneous across cancers but stopped short of quantifying why. We show that tumor proliferation is the strongest candidate driver of this heterogeneity (27.6\u201337.5% of between-cancer variance for the CD8 T-cell, CAF and IFN-\u03b3 signatures; permutation p < 0.05), consistent with a stoichiometric dilution of immune cells by proliferating malignant cells \u2014 a mechanism we further validate at single-cell resolution in a 481,738-cell pan-cancer atlas (malignant-cell fraction inversely correlated with immune-cell fraction, Spearman r = \u22120.885). Effect directions were reproduced in two independent cohorts (GLASS glioma; ICGC OV-AU ovarian cancer). The practical implication is that immune signatures should be corrected for proliferation and non-tumoral content before cross-cancer comparison or clinical deployment.')

para('Declarations.', bold=True, space_after=2)
para('This manuscript is original, has not been published previously, and is not under consideration elsewhere. All authors have read and approved the submission and declare no competing interests. The meta-analysis here is a statistical method (individual-data meta-analysis and meta-regression of a single public dataset, TCGA), not a systematic review of the literature; accordingly, a completed PRISMA-IPD checklist accompanies this submission, with items inapplicable to a single-dataset secondary analysis marked N/A. The manuscript contains 4 main figures, 3 supplementary figures, 2 main tables, 3 supplementary tables and 14 references. All analyses used only publicly available data (TCGA, GLASS, ICGC OV-AU and a pan-cancer single-cell atlas); code is archived on Zenodo (https://doi.org/10.5281/zenodo.24549465).')

para('Institutional standing.', bold=True, space_after=2)
para('The corresponding author\u2019s institution, Zhejiang Provincial People\u2019s Hospital, is a large tertiary hospital in China with a strong program in head-and-neck and thyroid cancer, and operates the Zhejiang Provincial Clinical Research Center for Malignant Tumors.')

para('We believe this work will be of interest to the journal\u2019s readership, and we would welcome the Editor\u2019s guidance on scope or format.')

para('Thank you for your consideration.', space_after=18)

para('Sincerely,', space_after=2)
para('Li-zhuo Zhang, MD, PhD (Corresponding author)', space_after=2)
para('Department of Head and Neck Surgery, Zhejiang Provincial People\u2019s Hospital, Hangzhou, China', space_after=2)
para('E-mail: zlz2020@foxmail.com; ORCID: 0000-0002-4460-1743', space_after=12)

para('on behalf of the authorship team:', space_after=2)
para('[authorship team \u2014 to be completed]')

d.save(FP)
print('已保存:', FP)

# 复检
d2 = Document(FP)
full = '\n'.join(p.text for p in d2.paragraphs)
print('段落数:', len(d2.paragraphs))
print('含 PRISMA-IPD:', 'PRISMA-IPD' in full)
print('含 4 main figures:', '4 main figures' in full)
print('含 Zenodo DOI:', 'zenodo.24549465' in full)
