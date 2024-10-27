## Scripts and Workflow
- Modeling 
  - `src/bayes_mcmc_nonlatent_3w.R` model of outcome regression with proxy. 
  - `src/bayes_mcmc_latent_u1w3.R` model of latent adjustment and no adjustment. 
  - `src/bayes_mcmc_latent_u1w3_random_effect.R` model with random effects. 
  - `src/helper.R`
- Redlining data clean 
  - Key output file is `/data/cleaned/final_census1940.rds` generated from `src/redlining_data_clean.R`. 
  - Spatial spline files are `/data/cleaned/RedliningSplineMat[].rds` generated from `src/spline.R`. 
- Simulation study
  - Simulation with geometry in grids, follow files named `simu[]_data.R`, `simu[]_latent_u1w3.R`, `simu[]_nonlatent_w3.R`, `simu[]_post.R` etc. 
    - use `simu[]_data.R` to generate simulation data, use `simu[]_latent_u1w3.R`, `simu[]_nonlatent_w3.R` to run Bayesian MCMC models, and use `simu[]_post.R` to collect all simulation results into single file. 
  - Simulation with geometry of redlining maps, follow files named `simured1_data.R`, `simured1_latent_u1w3.R`, `simured1_nonlatent_w.R`, `simured1_post.R` etc.
- Redlining data analysis
  - Scripts 
    - `src/redlining_model_latent_u1w3.R` execute model of latent adjustment and no adjustment.
    - `src/redlining_model_nonlatent_w3.R` execute model of outcome regression with proxy. 
    - `src/redlining_post.R` post-process after modeling, create summary/figures/tables
    - `src/redlining.csh`
- Paper preparation
  - `src/redlining_data_stats.R` collect data statistics in the paper. 
  - `src/simu_post.R` collect simulation results.
  

## Data directory
file `/data/cleaned/final_census1940.rds`

| Variable              | Description                                                        |
|-- | -- |
| `state`                 | State name                                                        |
| `city`                  | City name,                                                        |
| `SPL`                   | Unique numerical ID for each city                                 |
| `neighborho`            | Redlining region unit                                             |
| `holc_grade`            | Redlining Grades from A to D                                      |
| `holc`                  | Numerical redlining Grades, 1 for A, 2 for B, 3 for C, 4 for D    |   
| `pm25`                  | PM$_{2.5}$ concentration in 2010                                  |
| `no2`                   | NO$_2$ concentration in 2010                                      |
| `bc.rate.employed`      | Box-Cox transformed rate of employment. The untransformed one is `rate.employed` |
| `pct.black.pop40.rank`  | Rank normalized percentage of Black Population. The untransformed one is `pct.black.pop40` |
| `wt.mean.rent`          | Count-weighted monthly home rent ($) in the redlining region                             |

 