# Aggregate CellRanger outputs across all samples
cellranger aggr --id=aggr-all \
                --csv=aggr-all-new.csv \
                --normalize=none
