# pysb_on_hpc
Information of using pysb on greene hpc

## File you need to run pysb on hpc
1. sif file: https://drive.google.com/file/d/17PgFDEaqffqBuJ7qydHd8OQPM4DKg4h5/view?usp=sharing
2. py file: example.py
3. sbatch file: example.sbatch

## Steps for submitting jobs
### 1. The files needed for simulation
The provided sif file contains the packages required for running pysb. 

The py file is an example code to simulate the KRAS G12C nucleotide exchange assay from the publication below. The example code is testing how the GTP and GDP concentration affect the drug IC50( the example showing here is AMG-510). The "load_model_G12C_NEA()" is a function that define the network and the kinetic parameters in rule-based manner. The "NEA_AMG510" is the function mimics the protocol of experiments in the publication( All the steps done in experiments). The "search_amg_ic50" is a function that use bisection method to find IC50 of the drug with given conditions. More inofrmation about using PySB can be found in: https://pysb.org/ 
>Kopra, Kari, et al. "Thermal shift assay for small GTPase stability screening: evaluation and suitability." International journal of molecular sci

The sbatch file is the script submitting job to hpc. The sbatch file will first require resources from hpc, load sif file and then run py file. More detail about slurm job submitting could be found in tutorial: https://sites.google.com/nyu.edu/nyu-hpc/training-support/tutorials/slurm-tutorial

ps. don't use or load plotting package like Matplotlib in py file.

### 2. Run the simulation
a. Check py file is runable.

b. Chekc sbatch file has the reasonable settings of request for resources, using the right path to load the sif file, and the right path of py file.

c. Go the repository where I save the sbatch file. 

```
cd /path_to_sif_file
```
d. Submitting batch job
```
sbatch example.sbatch
```
e. If you want to check the status of your job. (the command will showing all the job status of the user with a window that will renew the status every 2 seconds)
```
watch -n 2 squeue -u <your_user_id>
```

## Citation

If you use PySB or the provided sif file in your work, please cite the following publication:

>Lopez, Carlos F., et al. "Programming biological models in Python using PySB." Molecular systems biology 9.1 (2013): 646.
