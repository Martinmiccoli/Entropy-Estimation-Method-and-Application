# 🔢 Entropy Estimators and Gene Networks

This repository contains a comprehensive implementation of various entropy estimation techniques and their applications in gene network analysis, including the **James-Stein shrinkage entropy estimator** and the ARACNE algorithm for network inference.

---

## ✨ Features

- 🚀 **Entropy Estimation**: Implementation of:
  - Maximum Likelihood (ML) Estimator
  - Miller-Madow Estimator
  - Bayesian Estimators with multiple priors
  - Chao-Shen Estimator
  - James-Stein Shrinkage Estimator
- 🔍 **Mutual Information**: Efficient calculation of mutual information for discrete data.
- 🧬 **Gene Network Analysis**:
  - ARACNE algorithm for pruning indirect interactions.
  - Visualization of gene interaction networks using `ggraph` and `igraph`.
- 📊 **Simulation Framework**: Generate synthetic datasets under various scenarios and compare estimator performance using metrics like MSE and bias.

---

## 📦 Installation

To use the tools provided in this repository, ensure you have R installed with the following libraries:
- `gtools`
- `ggplot2`
- `entropy`
- `gridExtra`
- `ggraph`
- `igraph`
- `infotheo`
- `parmigene`
- `GeneNet`

---

## 📊 Results

The project includes a detailed simulation framework to evaluate the performance of different entropy estimators under various scenarios, such as:

1) Dirichlet distribution with small alpha values.
2) Uniform Dirichlet distribution.
3) Sparse Dirichlet distribution with structural zeros.
4) Zipf power-law distribution.

---

## 📖 More Info
- [📑 Paper](https://github.com/Martinmiccoli/Entropy-Estimation-Method-and-Application/blob/a874716333f15567a82e7e92b8d6ea594e7e952d/12.EntropyEstimators.pdf)
- [🕸️ Ecoli Gene Network](https://github.com/Martinmiccoli/Entropy-Estimation-Method-and-Application/blob/a874716333f15567a82e7e92b8d6ea594e7e952d/ecoli_gene_network.pdf)
- [🕷️ Aracne Comparison](https://github.com/Martinmiccoli/Entropy-Estimation-Method-and-Application/blob/main/mi_aracne_comparison.png)
- [📈 PPT Presentation](https://github.com/Martinmiccoli/Entropy-Estimation-Method-and-Application/blob/a874716333f15567a82e7e92b8d6ea594e7e952d/PPT%20Entropy%20Inference%20and%20the%20James-Stein%20Estimator.pptx)
