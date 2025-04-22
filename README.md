## Selective Multiple Testing

### Context
We provide a one-stop collection of resources on covariate selective inference with Family-Wise Error Rate (FWER) control on large panel asset pricing data and model, as described in [Selective Multiple Testing: Inference for Large Panels with Many Covariates](https://papers.ssrn.com/sol3/Papers.cfm?abstract_id=4315891). Specifically, we enable users to perform the rolling-window estimations described in \S 7.3 of the paper. Any constant window model or alternative window size can be produced by changing the window parameter or estimation in the code.

* data: 
 1. Response variables: $Y=R^{660 \times 243}$ test portfolios excess returns downloaded and processed from [Kenneth French's website](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html). The cross-section size is 243 as test portfolios and there are 660 monthly observations. We regress out the market factor from each of the invidiual factors.
 2. Covariates: $X=R^{660 \times 114}$ asset pricing high-minus-low factors that are downloaded and processed from [Hou, K., Xue, C. and Zhang, L., 2020. Replicating anomalies. The Review of financial studies, 33(5), pp.2019-2133](https://global-q.org/index.html). The covariates' dimension is 114 and there are 660 monthly observations. We regress out the market factor from each of the individual factors.

* R code:
  * step1: `multiple_fixed_window.R` performs a rolling window regression to generate a $P=R^{114 \times 243}$ matrix for each rolling window. Specifically, `posi_lognorm_pval_enforce_dimension` admits a $(X,Y)$ tuple and optional priors, and returns the post-LASSO selection inference valid p-values from [Jason D. Lee. Dennis L. Sun. Yuekai Sun. Jonathan E. Taylor. "Exact post-selection inference, with application to the lasso." Ann. Statist. 44 (3) 907 - 927, June 2016](https://projecteuclid.org/journals/annals-of-statistics/volume-44/issue-3/Exact-post-selection-inference-with-application-to-the-lasso/10.1214/15-AOS1371.full);
  * step2: `select.R` performs the _Selective Multiple Testing_ selection described in the paper;
  * step3: `eval_performance_multiple_fixed_window.R` performs evaluations described in the paper.
* python code on selection:
  * python: `funs.py` provides minimal stand-alone function that only requires `pandas` and `numpy` to perform our _Selective Multiple Testing_ selection method given a matrix of p-values, controlling for Family-Wise Error Rate (FWER).

### Additional resources

For a method-focused code base, we provide a python version of the  _Selective Multiple Testing_ in the Github repository (https://github.com/yfan7/panel_CPD/tree/master) for our accompanying paper _Large Dimensional Change Point Detection with
FWER Control as Automatic Stopping_(https://drive.google.com/file/d/15SotyMqpWBUTrwaCpzNGron2F4uz1wdL/view).

### Usage

To cite this code, in addition to the data sources, please use the following citation:
```
@article{pelger2022inference,
  title={Inference for Large Panel Data with Many Covariates},
  author={Pelger, Markus and Zou, Jiacheng},
  journal={arXiv preprint arXiv:2301.00292},
  year={2022}
}
```
