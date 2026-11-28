# Verify sensitivity results: (1) B2_repro reproduces main B2; (2) sample changes from grade
main <- read.csv("results/specB/per_cancer_three.csv")
main_b2 <- main[main$version == "B2", c("cancer","signature","beta","se")]
sens <- read.csv("results/sens_grade_dss/per_cancer_sens.csv")
repro <- sens[sens$version == "B2_repro", c("cancer","signature","beta","se","n")]
grade <- sens[sens$version == "B2_grade", c("cancer","signature","beta","se","n")]

m <- merge(main_b2, repro, by = c("cancer","signature"), suffixes = c("_main","_repro"))
cat("beta max |diff| (main B2 vs repro):", max(abs(m$beta_main - m$beta_repro), na.rm = TRUE), "\n")
cat("se   max |diff| (main B2 vs repro):", max(abs(m$se_main - m$se_repro), na.rm = TRUE), "\n")
cat("rows main/repro:", nrow(main_b2), nrow(repro), "\n")

mg <- merge(main_b2, grade, by = c("cancer","signature"), suffixes = c("_main","_grade"))
mg <- mg[!duplicated(mg$cancer), c("cancer","n_main","n_grade")]
cat("\nPer-cancer n: main vs +grade (rows where n differs = grade was added):\n")
print(mg[mg$n_main != mg$n_grade, ])
cat("\nCancer types with grade added:", nrow(mg[mg$n_main != mg$n_grade, ]), "\n")
