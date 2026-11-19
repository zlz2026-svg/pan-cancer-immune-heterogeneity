# Fig1 研究设计流程图（matplotlib 纵向流程）
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(8, 10))
ax.set_xlim(0, 10); ax.set_ylim(0, 13)
ax.axis("off")

def box(x, y, w, h, text, fc, tc="black", fs=9.5, bold=False):
    b = FancyBboxPatch((x - w/2, y - h/2), w, h,
                       boxstyle="round,pad=0.08", fc=fc, ec="black", lw=1.2)
    ax.add_patch(b)
    ax.text(x, y, text, ha="center", va="center", fontsize=fs,
            color=tc, fontweight="bold" if bold else "normal", wrap=True)

def arrow(x1, y1, x2, y2):
    a = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                        mutation_scale=14, lw=1.4, color="black")
    ax.add_patch(a)

# 数据源
box(5, 12.2, 6.6, 1.1, "TCGA PanCancer Atlas\n33 cancer types\n(RNA-seq + survival)", "#e8f0fe", bold=True)
# 排除
box(5, 10.4, 6.6, 1.0, "Exclusion\nLAML (bone-marrow, no primary);\nTGCT (<5 death events)", "#fdecea")
# 31 癌种
box(5, 8.7, 6.6, 0.95, "31 cancer types\n8,924 primary tumours", "#e8f0fe", bold=True)
# 签名
box(5, 7.15, 6.6, 1.0, "7 immune/proliferation signatures\nwithin-cancer z-scored scores", "#e8f0fe")
# Cox
box(5, 5.6, 6.6, 1.0, "Per-cancer Cox regression\n(217 models; age, sex, stage where available)", "#e8f0fe")
# meta
box(5, 4.05, 6.6, 1.0, "Random-effects meta-analysis\npooled HR + I\u00b2 heterogeneity", "#e8f0fe")
# 归因
box(5, 2.5, 6.6, 1.15, "Univariate meta-regression attribution\nproliferation / germ layer /\nstromal-immune / TMB", "#fff3e0", bold=True)
# 验证
box(5, 0.85, 6.6, 0.95, "External validation\nGLASS glioma cohort (n = 79)", "#e8f0fe")

arrow(5, 11.65, 5, 10.95)
arrow(5, 9.9, 5, 9.2)
arrow(5, 8.2, 5, 7.7)
arrow(5, 6.65, 5, 6.15)
arrow(5, 5.1, 5, 4.6)
arrow(5, 3.55, 5, 3.1)
arrow(5, 1.95, 5, 1.35)

plt.tight_layout()
plt.savefig("figures/Fig1_flow.png", dpi=300, bbox_inches="tight", facecolor="white")
print("Fig1_flow.png saved")
