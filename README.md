# The Domain-Linker-Domain (DLD) Dataset

The **DLD dataset** identifies Independent Domain Linkers (IDLs) by extracting high-quality domain annotations from SCOP2, filtering for multi-domain proteins, and refining linker regions using DSSP-based secondary structure analysis.

A two-step smoothing process distinguishes **1,640 Independent Domain Linkers (IDLs)** and **647 Dependent Domain Linkers (DDLs)**, alongside other protein regions. The dataset is mapped to UniProt, providing a comprehensive resource for studying **disordered flexible linkers** and supporting **computational prediction methods**.

---

## 📂 Collecting the Dataset

The code files are numbered sequentially (`0_` to `14_`).  
To reproduce the DLD dataset, please run the code files in order:
E.g. 

```bash
python 1_*.py
...
python 14_*.py
```
