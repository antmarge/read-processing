# Read processing and read-based analyses for low-coverage data

## Margaret Antonio 2019.07.17


### Description
Basic pipeline for processing read data to get genotype
calls and/or pseudohaploid calls. 

#### Running the pipeline
1. Go to snakemake_variables.py. Check all filepaths (reference files, software, etc.) to make sure they are valid.
2. Create the sample_list file with your sample ids and filepaths
3. Do a dry run of the snakemake pipeline with `snakemake -np`
4. Check cluster configurations in sm_script.sbatch and sm_slurm_config.json. The files included were written for Stanford's Sherlock slurm cluster. 
5. Suggested: run pipeline in a interactive session in a screen using the line included in sm_script.sbatch


#### Things to add / do
1. Example files
2. Complete final LASER rule
3. Other read-based analyses
4. Pseudohaploid genotypes
 
