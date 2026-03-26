# MOSAIC eDNA Occupancy Model — Explainer and Methods

---

## Part 1: Plain-Language Explainer

### The core question

When we collect a water sample and run eDNA analysis, a negative result doesn't necessarily mean an animal wasn't there. There are at least two separate ways we can "miss" a species: the animal may not have shed detectable eDNA into that particular water sample at all, or the eDNA may have been present in the water but simply not amplified in one or more of the PCR replicates we ran in the lab. These are distinct biological and technical processes, and a good model needs to account for both.

### The two-level model

We use a **two-level occupancy model** that separates these two sources of uncertainty.

**Level 1 — Was the eDNA captured in the filter?**
The first level estimates the probability that a given water filter actually captured eDNA from a given marine mammal family. We call this the *presence/capture probability* (φ, "phi"). It answers the question: given all the environmental and methodological conditions of this particular sample, how likely is it that eDNA from this family ended up in the filter?

**Level 2 — Did PCR detect it?**
Given that eDNA *was* captured in the filter, the second level estimates the probability that any individual PCR replicate successfully amplified and detected it. We call this the *detection probability per replicate* (p). This reflects lab-level variability: even from the same filter extract, some replicates may succeed and others fail.

The observed data for each filter is simply the number of PCR replicates (out of typically 4) that returned a positive result.

### What drives each level?

**Presence/capture (φ)** is modelled as a function of several covariates:

- **Species identity**: each of the six focal families (dolphins, porpoises, baleen whales, sea lions, seals, and mackerel sharks) has its own baseline capture probability, reflecting real differences in their abundance and eDNA shedding rates.
- **Cruise**: a fixed adjustment accounts for the fact that the two surveys were conducted in different years, when overall marine mammal abundance may have differed.
- **Coastal vs. offshore**: each family gets its own inshore/offshore offset, because the spatial distribution of nearshore versus pelagic species differs substantially.
- **Collection method**: MOSAIC-D (CTD cast) is the reference method. MOSAIC-A and MOSAIC-B (direct in-line filtration) are treated as equivalent to each other and together receive a single contrast relative to MOSAIC-D. MOSAIC-C (carboy collection) receives a separate contrast relative to MOSAIC-D.
- **Filter volume**: more water filtered generally means more eDNA captured.
- **Distance from shore**: a continuous measure of how far offshore the sample was collected.

**Detection per replicate (p)** is modelled as a function of:

- **Filtration time**: longer filtration times may affect eDNA integrity and concentration, influencing PCR success. This effect is shared across all families.

### Sharing information across species

A key feature of the model is that it is *hierarchical*: the six families are not modelled completely independently. Instead, their baseline capture probabilities (and inshore/offshore offsets) are assumed to come from a common distribution. This means that information about how common or rare marine predators generally are in this system is shared across species. Rare families with few detections borrow strength from more commonly detected ones, producing more stable estimates than if each family were fitted alone. The spread of that shared distribution is itself estimated from the data.

### How uncertainty is handled

The model is fully Bayesian: rather than producing single point estimates, it produces a full posterior distribution over every unknown quantity. This naturally propagates all sources of uncertainty — from the hyperparameters down through the family-level effects and all the way to the filter-level probabilities — into the final estimates. Model fitting used Markov Chain Monte Carlo (MCMC) sampling via NIMBLE, run as three independent chains to verify convergence.

---

## Part 2: Methods Section

### Statistical model

We modelled eDNA detections of six focal marine families -5 marine mammals, and one family of sharks (*Delphinidae*, *Phocoenidae*, *Balaenopteridae*, *Otariidae*, *Phocidae*, and *Alopiidae) using a two-level hierarchical Bayesian occupancy model implemented in NIMBLE (de Valpine et al. 2017). The model simultaneously estimates the probability that eDNA from a given family was present and captured on a given filter (the *presence/capture* submodel) and the probability that a single PCR technical replicate detected that eDNA given it was captured (the *detection* submodel). Separating these two processes is necessary because a zero count of positive replicates can arise either from the absence of eDNA in the sample or from PCR failure despite its presence.

#### Observation model

Let $i = 1, \ldots, N$ index water filters, $f = 1, \ldots, F$ index families, and $r_i$ denote the number of PCR technical replicates processed for filter $i$. The observed count of positive replicates $n_{\text{det},i,f}$ is modelled as:

$$z_{i,f} \sim \text{Bernoulli}(\phi_{i,f})$$

$$p^{\text{eff}}_{i,f} = p_i \cdot z_{i,f}$$

$$n_{\text{det},i,f} \sim \text{Binomial}\!\left(p^{\text{eff}}_{i,f},\ r_i\right)$$

where $z_{i,f} \in \{0, 1\}$ is a latent binary variable indicating whether family $f$ eDNA was present and captured on filter $i$, $\phi_{i,f}$ is the corresponding presence/capture probability, and $p_i$ is the per-replicate detection probability. The deterministic node $p^{\text{eff}}_{i,f}$ equals $p_i$ when eDNA is present ($z_{i,f} = 1$) and collapses to zero when it is absent ($z_{i,f} = 0$), ensuring that positive replicates are impossible in the absence of captured eDNA.

#### Presence/capture submodel

The logit-scale linear predictor for $\phi_{i,f}$ includes taxon-level random effects, a cruise fixed effect, and covariates shared across taxa:

$$\text{logit}(\phi_{i,f}) = \alpha_f + d_{\text{cruise}} \cdot \text{cruise2}_i + \beta_f \cdot \text{inshore}_i + d_{AB} \cdot \mathbb{1}[\text{MOSAIC-A/B}]_i + d_C \cdot \mathbb{1}[\text{MOSAIC-C}]_i + d_{\text{vol}} \cdot \widetilde{\text{vol}}_i + d_{\text{dist}} \cdot \widetilde{\text{dist}}_i$$

Here $\alpha_f$ is a taxon-specific intercept; $\text{cruise2}_i \in \{0,1\}$ is a binary indicator for the second survey cruise (with the first cruise as reference), and $d_{\text{cruise}}$ is the corresponding log-odds contrast; $\beta_f$ is a taxon-specific inshore offset, with $\text{inshore}_i = 1$ for samples collected in coastal waters and 0 for offshore; $\mathbb{1}[\text{MOSAIC-A/B}]_i$ is an indicator equal to 1 for both MOSAIC-A and MOSAIC-B (direct in-line filtration methods treated as equivalent), with shared contrast $d_{AB}$ relative to the reference method MOSAIC-D (CTD cast); $d_C$ is the contrast for MOSAIC-C (carboy filtration) relative to MOSAIC-D; and $\widetilde{\text{vol}}_i$ and $\widetilde{\text{dist}}_i$ are standardised (zero-mean, unit-variance) filter volume and distance from shore, with shared log-odds slopes $d_{\text{vol}}$ and $d_{\text{dist}}$.

#### Detection submodel

Per-replicate detection probability is modelled as a function of filtration time, shared across all taxa:

$$\text{logit}(p_i) = \delta + d_{\text{time}} \cdot \widetilde{\text{time}}_i$$

where $\widetilde{\text{time}}_i$ is standardised filtration time and $\delta$ is the baseline log-odds of detection.

#### Hierarchical priors

Taxon-specific intercepts and inshore offsets are drawn from shared distributions, enabling partial pooling across families:

$$\alpha_f = \mu_\alpha + \sigma_\alpha \cdot \tilde{z}_{\alpha,f}, \quad \tilde{z}_{\alpha,f} \sim \mathcal{N}(0,\, 1)$$

$$\beta_f = \mu_\beta + \sigma_\beta \cdot \tilde{z}_{\beta,f}, \quad \tilde{z}_{\beta,f} \sim \mathcal{N}(0,\, 1)$$

This non-centred parameterisation (Papaspiliopoulos et al. 2007) decouples the hyperparameters from the family-level deviates, improving MCMC mixing in the presence of the characteristic funnel-shaped geometry of hierarchical posteriors.

Hyperpriors were specified as:

$$\mu_\alpha \sim \mathcal{N}(0,\, 1.5), \quad \sigma_\alpha \sim \text{Exponential}(1)$$

$$\mu_\beta \sim \mathcal{N}(0.0,\, 1), \quad \sigma_\beta \sim \text{Exponential}(1)$$

Weakly informative priors were placed on all shared coefficients:

$$d_{AB},\ d_C,\ d_{\text{vol}},\ d_{\text{dist}},\ d_{\text{cruise}} \sim \mathcal{N}(0,\, 1)$$

$$\delta \sim \mathcal{N}(-1,\, 1), \quad d_{\text{time}} \sim \mathcal{N}(0,\, 1)$$

The prior on $\delta$ is centred at $-1$ on the logit scale (corresponding to a baseline per-replicate detection probability of approximately 27%), reflecting the expectation that individual PCR replicates have moderate rather than near-certain detection rates.

#### Model fitting

Posterior distributions were estimated by MCMC using the NIMBLE software package (v1.4.1; de Valpine et al. 2017) in R. Three independent chains were run for 120,000 iterations each, with the first 20,000 discarded as burn-in and the remainder thinned by a factor of 20, yielding 15,000 posterior samples per chain (45,000 total). Convergence was assessed using the Gelman–Rubin potential scale reduction factor ($\hat{R}$; Gelman & Rubin 1992), with $\hat{R} < 1.05$ for all monitored parameters considered indicative of adequate convergence. Effective sample sizes were also inspected to confirm sufficient posterior exploration.

#### References

de Valpine, P., Turek, D., Paciorek, C. J., Anderson-Bergman, C., Temple Lang, D., & Bodik, R. (2017). Programming with models: writing statistical algorithms for general model structures with NIMBLE. *Journal of Computational and Graphical Statistics*, 26(2), 403–413.

Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation using multiple sequences. *Statistical Science*, 7(4), 457–472.

Papaspiliopoulos, O., Roberts, G. O., & Sköld, M. (2007). A general framework for the parametrization of hierarchical models. *Statistical Science*, 22(1), 59–73.
