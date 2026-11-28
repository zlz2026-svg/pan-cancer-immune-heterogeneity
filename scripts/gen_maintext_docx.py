#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""正文 V16 md → Word 投稿稿（忠实转换，学术论文格式，npj 初投版）
- 去版本注释块、分隔线
- # → Title; ## → Heading1; ### → Heading2
- markdown 表格 → Word 表格
- [^N] 正文引用 → 上标 [N]；[^N]: → 参考文献编号条目
- **bold** / *italic* / `code` inline 处理
"""
import sys, re
sys.stdout.reconfigure(encoding='utf-8')
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn

SRC = r'E:\1A-AI work\泛癌免疫异质性归因\02-最终交付\【1120】正文V17.md'
FP = r'E:\1A-AI work\泛癌免疫异质性归因\投稿材料包\npjPrecisionOncology\02-正文MainText\【1128】正文V17投稿稿.docx'

lines = open(SRC, encoding='utf-8').read().split('\n')

# 删除版本注释块（> 开头的行），保留标题和正文
lines = [l for l in lines if not l.strip().startswith('>')]

d = Document()
# 页面边距（npj 接收后要求，初投宽松）
for sec in d.sections:
    sec.left_margin = Inches(1)
    sec.right_margin = Inches(1)
    sec.top_margin = Inches(1)
    sec.bottom_margin = Inches(1)

# 默认字体
normal = d.styles['Normal']
normal.font.name = 'Times New Roman'
normal.font.size = Pt(11)
normal._element.rPr.rFonts.set(qn('w:eastAsia'), '宋体')

# 标题样式
for sname, sz in [('Heading 1', 12), ('Heading 2', 11)]:
    st = d.styles[sname]
    st.font.name = 'Times New Roman'
    st.font.size = Pt(sz)
    st.font.bold = True
    st.font.color.rgb = RGBColor(0, 0, 0)

INLINE_CODE = re.compile(r'`([^`]+)`')

def add_runs(p, text):
    """解析 **bold** / *italic* / `code` / [^N] 引用，逐段加 run"""
    # 先处理 [^N] 脚注引用 -> 上标 [N]
    parts = re.split(r'(\[\^\d+\])', text)
    for part in parts:
        if not part:
            continue
        m = re.match(r'\[\^(\d+)\]$', part)
        if m:
            r = p.add_run(f'[{m.group(1)}]')
            r.font.superscript = True
            continue
        # 处理 **bold**
        sub_parts = re.split(r'(\*\*[^*]+\*\*)', part)
        for sp in sub_parts:
            if not sp:
                continue
            if sp.startswith('**') and sp.endswith('**'):
                r = p.add_run(sp[2:-2]); r.bold = True
                continue
            # 处理 `code`
            code_parts = INLINE_CODE.split(sp)
            for i, cp in enumerate(code_parts):
                if not cp:
                    continue
                if i % 2 == 1:
                    r = p.add_run(cp); r.font.name = 'Courier New'
                    r.font.size = Pt(10)
                else:
                    # 处理 *italic*（单星号，排除表格里的，此处正文）
                    it_parts = re.split(r'(\*[^*\n]+\*)', cp)
                    for ip in it_parts:
                        if not ip:
                            continue
                        if ip.startswith('*') and ip.endswith('*') and len(ip) > 2:
                            r = p.add_run(ip[1:-1]); r.italic = True
                        else:
                            p.add_run(ip)

def add_md_table(rows):
    """rows: list of list of str"""
    n_cols = max(len(r) for r in rows)
    tbl = d.add_table(rows=len(rows), cols=n_cols)
    tbl.style = 'Table Grid'
    for i, row in enumerate(rows):
        for j in range(n_cols):
            cell = tbl.rows[i].cells[j]
            txt = row[j].strip() if j < len(row) else ''
            cell.text = txt
            for para in cell.paragraphs:
                for run in para.runs:
                    run.font.size = Pt(9)
                    run.font.name = 'Times New Roman'
            if i == 0:
                for para in cell.paragraphs:
                    for run in para.runs:
                        run.bold = True
    return tbl

# 逐行解析
i = 0
n = len(lines)
in_table = []
table_rows = []

def flush_table():
    global table_rows
    if table_rows:
        # 去分隔行 |---|
        rows = [r for r in table_rows if not re.match(r'^\|[\s:\-|]+\|$', r)]
        if rows:
            add_md_table([re.split(r'(?<!\\)\|', r.strip().strip('|')) for r in rows])
        table_rows = []
        d.add_paragraph()  # 表格后空行

while i < n:
    line = lines[i].rstrip()
    s = line.strip()

    # 表格行
    if s.startswith('|'):
        if table_rows and not s.startswith('|'):
            flush_table()
        table_rows.append(s)
        i += 1
        continue
    else:
        flush_table()

    # 空行 / 分隔线
    if s == '':
        i += 1
        continue
    if s == '---':
        i += 1
        continue

    # 标题
    if s.startswith('### '):
        d.add_heading(s[4:], level=2)
        i += 1
        continue
    if s.startswith('## '):
        d.add_heading(s[3:], level=1)
        i += 1
        continue
    if s.startswith('# '):
        # 主标题
        p = d.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r = p.add_run(s[2:])
        r.bold = True
        r.font.size = Pt(14)
        i += 1
        continue

    # 参考文献 [^N]:
    m = re.match(r'^\[\^(\d+)\]:\s*(.*)$', s)
    if m:
        p = d.add_paragraph()
        p.paragraph_format.left_indent = Inches(0.3)
        p.paragraph_format.first_line_indent = Inches(-0.3)
        r = p.add_run(f'{m.group(1)}. ')
        add_runs(p, m.group(2))
        p.paragraph_format.space_after = Pt(4)
        i += 1
        continue

    # 表格标题 **Table ...**
    if s.startswith('**') and ('Table' in s or 'Supplementary Figure' in s or s.startswith('**Fig.')):
        p = d.add_paragraph()
        add_runs(p, s)
        i += 1
        continue

    # 普通段落
    p = d.add_paragraph()
    p.paragraph_format.space_after = Pt(6)
    add_runs(p, s)
    i += 1

flush_table()
d.save(FP)
print('已保存:', FP)

# 复检
d2 = Document(FP)
print('段落数:', len(d2.paragraphs))
print('表格数:', len(d2.tables))
full = '\n'.join(p.text for p in d2.paragraphs)
print('含 ABSTRACT:', 'ABSTRACT' in full)
print('含 REFERENCES:', 'REFERENCES' in full)
print('含 AI USE:', 'AI USE DECLARATION' in full)
print('参考文献条数(编号段):', len(re.findall(r'^\d+\. ', full, re.M)))
