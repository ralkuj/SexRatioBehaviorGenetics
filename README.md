# SexRatioBehaviorGenetics
Code related to the Brief Report

Harper K, Kuja-Halkola R, Verweij KJH, Magnusson PKE, Zietsch BP. No Evidence of Genetic Basis to Variation in Human Offspring Sex Ratio. Beh Gen 2026.

Rebutting some of the conclusions in the article

Song S, Zhang J. In search of the genetic variants of human sex ratio at birth: was Fisher wrong about sex ratio evolution? Proc Biol Sci. Oct 2024;291(2033):20241876. doi:10.1098/rspb.2024.1876

Which includes a critique of the conclusions from paper 

Zietsch BP, Walum H, Lichtenstein P, Verweij KJH, Kuja-Halkola R. No genetic contribution to variation in human offspring sex ratio: a total population study of 4.7 million births. Proc Biol Sci. 2020;287(1921):20192849. DOI: 10.1098/rspb.2019.2849

Software version used:
R version 4.3.2 (2023-10-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux 9.4 (Plow)

Authors of code shared:
Ralf Kuja-Halkola, 01, 02, 03, 04.
Kaitlyn Harper, 05.

The simulations in first files (01) is free of external data. The simulation in second and fourth file (02 and 04) uses data on family structure from Zietsch et al, this data cannot be shared openly.

Files and short description:

01.SNP_simulations_Song&Zhang_20250225.R
- Implements the simulation from Song&Zhang to estimate the statistical power to detect a non-zero heritability olf offspring sex-ratio. Implements the original wrong p-value threshold by Song&Zhang  of 5e-8 and a correct threshold of 0.05.

02.SNP_simulations_FamilyStructure_20251209.R
- Implements the one SNP approach of Song&Zhang, but using real family structure from Zietsch et al. This yield power calculations that are trustworthy. First all male-male siblings' offspring simulated, then firstborn to all male-male siblings.

03.SNP_simulations_Plotting_20251210.R
- Plots the above simulations.

04.SexRatio_simulations_Normal_liability_20250225.R
- Performs a simulation using the relative structure from Zietsch et al to find statistical power for detecting non-null heritability under assumed model with offspring sex allocation drawn from a bernoulli ditribution with probability distributed according to a normal distribution centered aroun 0.5 and varying variance. The variances include the variance used in Song&Zhang, but also smaller and larger variances.

05.SZ_SLiM_Plotting.R
- Plots Song and Zhang's SLiM Fisher's Principle simulation results.

Simulation results:
SimulationSNP_1000repeats_SongOriginal_wrong_pval_20241105.Rds
- Simulation  results based on Song&Zhangs original simulation with incorrect p-value threshold. Produced in 01.SNP_simulations_Song&Zhang_20250225.R.

SimulationSNP_1000repeats_SongOriginal_correct_pval_20241105.Rds
- Simulation based on Song&Zhangs original simulation with correct p-value threshold.  Produced in 01.SNP_simulations_Song&Zhang_20250225.R.

SimulationSNP_FamilyStructure_male_All_1000repeats_20251209.Rds
SimulationSNP_FamilyStructure_male_All_Xtra0075_1000repeats_20251209.Rds
- Simulation based on Song&Zhangs one SNP approach, using one known family structure and all offspring. Produced in 02.SNP_simulations_FamilyStructure_20251209.R.

SimulationSNP_FamilyStructure_male_First_1000repeats_20251209.Rds
SimulationSNP_FamilyStructure_male_First_Xtra0075_1000repeats_20251209.Rds
- Simulation based on Song&Zhangs one SNP approach, using one known family structure and only first-born offspring. Produced in 02.SNP_simulations_FamilyStructure_20251209.R.

Simulation_Realdata_Male_1000repeats_20241107.Rds
- Simulation mirroring the Zietsch et al paper, using known family structure, for different assumed variances of offspring sex probability. Produced in 04.SexRatio_simulations_Normal_liability_20250225.R.

SongZhang-original_(height_sims)summary.tsv
- Song & Zhang's original height simulations results.

SongZhang-original_Simulation_SLiM_human_directional_combined.csv
- Song & Zhang's original SLiM Fisher's Principle simulation results (combined per the script 05.SZ_SLiM_Plotting).
