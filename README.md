<h2>Elastic Integrative Analysis for Heterogeneous Treatment Effect</h2>


A test-based dynamic borrowing framework combining
a randomized clinical trial (RCT) and a real-world evidence (RWE) study,
in which a preliminary test statistic is used to gauge the comparability
and reliability of the RWE and to decide whether or not to use the RWE in
an integrative analysis. The parameter of interest is $\psi$,
which quantifies how the treatment effect varies over the treatment
modifiers.


<h4>Usage</h4>
<pre>
elasticHTE(data.rct, data.rwe, ...,
           outcome.type = c("cont", "bin"),
           ps.rct = NULL,
           sieve.degree = 2L,
           outcome.method = c("glm", "SL"),
           outcome.controls = list(family = "gaussian"),
           ps.method = c("glm", "SL"),
           ps.controls = list(family = "quasibinomial"),
           n.pert = 100L,
           fixed = FALSE,
           n.gamma = 1000L,
           n.boot = 100L,
           thres.psi = NULL)
</pre>

<h5>Abbreviated Formal Argument Descriptions</h5>

- data.rct: The value object returned by *dataInput()* for the
data from a randomized clinical trial (RCT). 
- data.rwe: The value object returned by *dataInput()* for the
data from a real-world evidence (RWE) study. 
- ...: Ignored
- outcome.type: The type of outcome. 
- ps.rct: Optional input providing a vector of known propensity
  scores P(A=1) for the RCT dataset.
- sieve.degree: The order of the polynomial defining the sieve model. 
- outcome.method: The regression method for outcomes.
- outcome.controls: Additional inputs provided to the specified outcome regression method.
- ps.method: The regression method for propensity score analysis.
- ps.controls: Additional inputs provided to the specified propensity score regression method.
- n.pert: The number of perturbations to generate when
estimating the variance.
- fixed: How to select the tuning parameter $c_{\gamma}$. FALSE, the default, selects an adaptive
selection strategy; TRUE selects a fixed threshold strategy.
- n.gamma: The number of samples to generate to estimate
$c_{\gamma}$, the threshold.
- n.boot: The number of bootstrap samples to generate
when estimating the confidence intervals.
- thres.psi: The threshold for constructing
adaptive confidence intervals.

<h5>Abbreviate Returned Object Description</h5>

A list with components:

- psi: The estimated \eqn{\psi: `psi`} associated
with the treatment modifiers under various models.
- ve: The estimated standard error for
\eqn{\psi: `psi`}.
- CIs.inf, CIs.sup: The estimated confidence intervals for `psi`.
- CI.settings: A list of the settings used in the confidence interval
procedure.
- Tstat: The estimated test statistic.
- conservative: $I(Tstat < thres.psi)$
- nuispar: A list providing the selected $\gamma$ and
its corresponding threshold value $c_{\gamma}$; $I(c_{\gamma}  > Tstat)$ and its
p-value; $\eta$, where Tstat = $\eta^T \eta$; and a list of the settings used in the selection procedure.

<h3>Examples</h3>

<pre>
  
# load provided illustrative toy dataset with continuous outcome
data("elasticToy.cont")

# conduct the elastic integrative analysis with defaults
result.cont <- elasticHTE(data.rct = dataInput(elasticToy.cont.rct,
                                               outcome.model = Y ~ (X1+X2)*A,
                                               ps.model = A ~ X1 + X2),
                          data.rwe = dataInput(elasticToy.cont.rwe,
                                               outcome.model = Y ~ (X1+X2)*A,
                                               ps.model = A ~ X1 + X2))

# load provided illustrative toy dataset with binary outcome
data("elasticToy.bin")

# conduct the elastic integrative analysis with defaults
result.bin <- elasticHTE(data.rct = dataInput(elasticToy.bin.rct,
                                              outcome.model = Y ~ (X1+X2)*A,
                                              ps.model = A ~ X1 + X2),
                         data.rwe = dataInput(elasticToy.bin.rwe,
                                              outcome.model = Y ~ (X1+X2)*A,
                                              ps.model = A ~ X1 + X2),
                         outcome.type = "bin")
</pre>

<h3>Install</h3>

In your R console

<pre>
library(devtools)
install_github("IntegrativeStats/ElasticIntegrative")
</pre>
