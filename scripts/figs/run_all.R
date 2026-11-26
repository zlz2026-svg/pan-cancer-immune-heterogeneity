# run_all.R -- rebuild every manuscript figure from the raw analysis outputs.
#
#   cd <project root> ; "E:/R/R-4.3.3/bin/Rscript.exe" scripts/figs/run_all.R
#
# Order matters only through 01_data.R, which every figure script sources itself;
# each step is run in its own process so a failure cannot leave a half-written
# device open. Any non-zero exit stops the chain.

root <- getwd()
if (!file.exists(file.path(root, "data", "cancer_meta.csv")))
  stop("Run this from the project root (data/cancer_meta.csv not found).")

rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) rscript <- file.path(R.home("bin"), "Rscript")

steps <- c(
  "scripts/figs/01_data.R",   # data layer + 12 assertions
  "scripts/figs/fig1.R",      # design: flow, cohort landscape, hypotheses
  "scripts/figs/fig2.R",      # 31 x 7 HR matrix + pooled strip
  "scripts/figs/fig3.R",      # attribution heatmap + 4 meta-regression panels
  "scripts/figs/fig4.R",      # direction reversals + GLASS validation
  "scripts/figs/figS2.R",     # robustness: alternative metric, stage, leave-one-out
  "scripts/figs/figS1.R",     # supplementary: dilution schematic
  "scripts/figs/figGA.R"      # graphical abstract: three-column narrative
)

run <- function(cmd, args, label) {
  cat("\n=== ", label, " ===\n", sep = "")
  t0 <- Sys.time()
  status <- system2(cmd, args, stdout = TRUE, stderr = TRUE)
  took <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  cat(paste(status, collapse = "\n"), "\n")
  bad <- attr(status, "status")
  if (!is.null(bad) && bad != 0) {
    cat("FAILED after", took, "s\n")
    quit(status = bad)
  }
  cat("ok (", took, "s)\n", sep = "")
}

for (s in steps) run(rscript, s, s)

# image-level QA: size, embedded font sizes, background, colour-vision safety
qa <- Sys.which("python")
if (nzchar(qa)) {
  run(qa, "scripts/figs/99_qa.py", "99_qa.py")
  # rendered geometry: glyph boxes must not overlap, leave the page, or repeat
  run(qa, "scripts/figs/98_audit.py", "98_audit.py")
} else {
  cat("\npython not on PATH -- skipping 99_qa.py and 98_audit.py\n")
}

cat("\nAll figures rebuilt into figures/.\n")
