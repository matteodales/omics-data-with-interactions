# Modeling interaction effects with multi-omics data

This repository hosts the code for the novel statistical method described in my master thesis _"Modeling interaction effects with multi-omics data"_, for multi-source omics data involving interactions.

# Abstract

Modeling with multi-omics data, consisting for example of genomics,
proteomics or metabolomics, presents multiple challenges such as the high-dimensionality of the problem ($p \gg n$) and the need for integration between multiple data sources. In this thesis, a new interaction model that includes multiple sources of data is developed from the integration of two existing methods, Pliable Lasso and Cooperative Learning. The new model is tested both on simulation studies and on real multi-omics datasets for labor onset and cancer treatment response prediction. In both settings, it outperforms pliable lasso paired with other methods of data fusion in terms of prediction performance and is able to perform variable selection of relevant features and interactions.

**Keywords:** Cooperative learning, pliable lasso, interaction models, multi-omics, personalized medicine.

# About the repository

* The "implementation" folder contains the implementation of the model.
* The "simulation_studies" folder contains the code to reproduce the results of the simulation studies section
* The "real_multiomics_studies" folder contains the code to reproduce the results of the real multiomics studies section
* The "figures" folder contains the figures present in the thesis and the code to generate them

# References

Daisy Yi Ding, Shuangning Li, Balasubramanian Narasimhan, Robert Tibshirani. Cooperative Learning for Multi-view Analysis. 2022. [[link]](https://arxiv.org/abs/2112.12337)

Robert Tibshirani, Jerome Friedman. A Pliable Lasso. 2019. [[link]](https://arxiv.org/abs/1712.00484)
