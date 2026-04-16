# DEG bootstrapping analysis
Code to run bootstrapping based DEG analysis as well as visualize the results.

## Run Trials
1. Run the DEG bootstrapping trials:
    1. Run `seq 20 | xargs -n 1 -P 10 Rscript A_DEG_Bootstrapping_ttmt.R`
    2. Run `seq 20 | xargs -n 1 -P 10 Rscript A_DEG_Bootstrapping_geno.R`
3. Run the DEG data collater
    1. Collate treatment DEGs using: `Rscript B_DEG_Collater_ttmt.R`
    2. Collate genotype DEGs using: `Rscript B_DEG_Collater_geno.R`
4. Plot and save the results
    1. Plot treatment DEGs by following the instructions in `C_DEG_Plotting_ttmt.ipynb`
    2. Plot genotype DEGs by following the instructions in `C_DEG_Plotting_geno.ipynb`
