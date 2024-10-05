<h5>Development Stage: <span style="color:red">Research</span></h5>

<h2>Estimating Weighted-Adjusted RMSTs in Causal Survival Analysis</h2>

1) Calculate the calibration weights `EW.p.M(df, varX, gf, M)`
2) Estimate the proposed weighted-adjusted RMSTs:
   - weighted Kaplan-Meier RMST estimator `W.KM.Est(df, p, tau)`
   - weighted G-Formula RMST estimator `W.GF.Est(df, p, varX, gf, tau)`
   - weighted Hajek RMST estimator `W.HJ.Est(df, p, tau)`
   - weighted Augmented RMST estimator `W.AG.Est(df, p, varX, gf, tau)`

<h3>Install</h3>

In your R console
<pre>
library(devtools)
install_github("https://github.com/IntegrativeStats/CausalRMST")
</pre>

<h3>Examples</h3>

see `Example Code.R`

<h3>Reference</h3>

Hua, K., Hong, H., & Wang, X. (2024), ``Inference of treatment effect and its regional modifiers using restricted mean survival time in multi-regional clinical trials" arXiv:2404.08128.
