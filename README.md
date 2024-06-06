# Penalty_Based_Lagrangian_Bilevel

- `hyperparam_opt` folder contains experiments on SVM hyper parameter optimization.
  - Three algorithms are in py files in `algorithms` folder.
  - Datasets include diabetes dataset and fourclass dataset.
  - Please run the lasted jupter notebook.

- `Transportation` folder contains experiments on some Network planning problems.
- `toy_example.ipynb` presents 2 toy examples to show the penalty-based Lagrangian Algorithm can effectively work.
/%
- `SVM_diabetes.ipynb`
  is an experiment that Liuyuan is working on. Liuyuan is trying to reperform the experiments on a SVM hyperparameter optimization problem as discussed both by
  LV-HBA [Yao et al., 2023] and GA [Xu and Zhu, 2023]. The details of their experiments can be seen also at
  https://github.com/SUSTech-Optimization/LV-HBA
  https://github.com/xsy786912649/Efficient-gradient-approximation-method-for-constrained-bilevel-optimization-problem
  `diabete.txt` is the data file that the SVM experiment is performed on.
- Folder `hyperparam_opt` contains all the code associated with the hyperparameter optimization experiments. More precisely, the algorithms are coded inside the `algorithms` folder, and called from the `SVM_Tests.ipynb` Jupyter Notebook to create the figures from the paper.

%\
