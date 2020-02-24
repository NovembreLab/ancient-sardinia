# ancient-sardinia

**@authors: Joe Marcus and Harald Ringbauer**

This repo holds the snakemake pipeline and scripts used to analyze ancient DNA data from Sardinia in Marcus et al. 2020.

## Install python packages:

Either use conda, or install requeire packages with **pip**:

```
pip3 install --user scipy numpy pandas snakemake cython pysam click xlrd xlwt openpyxl XlsxWriter jupyter pyfaidx scikit-learn matplotlib seaborn scikit-allel shapely
```

submit-snakemake.sh             # Loads conda environment with packages (if available)
submit-snakemake_noconda.sh     # Loads python packages internally
