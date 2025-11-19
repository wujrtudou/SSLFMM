# SSLFMM

Semisupervised Learning under a Mixed Missingness Mechanism in Finite Mixture Models.

---

## Overview

`SSLFMM` implements a semi-supervised learning framework for finite mixture models under a **mixed label-missingness mechanism**. The package is designed for settings where class labels are partially observed and the missingness mechanism is **informative**, combining:

- **Missing Completely at Random (MCAR)**, and  
- **Entropy-based Missing At Random (MAR)**.

Estimation is performed using an **Expectation–Conditional Maximisation (ECM)** algorithm with robust initialisation routines for stable convergence. The package provides tools for:

- simulating partially labelled data under a mixed MCAR / entropy-MAR mechanism,  
- fitting semi-supervised finite mixture models under this mechanism,  
- computing theoretical Bayes classification error rates, and  
- extracting posterior cluster probabilities and entropy measures.

The methodology is related to the statistical perspective and informative missingness behaviour discussed in:

- Ahfock and McLachlan (2020) \<doi:10.1007/s11222-020-09971-5\>  
- Ahfock and McLachlan (2023) \<doi:10.1016/j.ecosta.2022.03.007\>

---

## Installation

### From CRAN (once available)

```r
install.packages("SSLFMM")
