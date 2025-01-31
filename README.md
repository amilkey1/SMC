# sequential multispecies coalescent
This program uses a sequential Monte Carlo approach to sample the posterior distribution of trees under the multispecies coalescent model.

The file proj.conf is used to specify program settings:

## The following settings apply to any start mode.
`seed`: specify random number seed. Default is 1. \
`theta`: specify the population-scaled mutation rate. If estimating theta, this will be the mean of the thetas drawn for each population in the species tree.\
`lambda`: specify the speciation rate\
`subset`: specify the gene partitions in the data file. For example:
* `subset = gene1:1-1000`
* `subset = gene2:1001-2000`\
  specifies that bases 1-1000 correspond to gene 1 and bases 1001-2000 correspond to gene 2. If no subset is specified, the entire file will be read as one gene

`model`: the model to use in first-level SMC. Current options are JC or HKY\
`kappa`: if HKY is chosen, set kappa\
`base_frequencies`: if HKY is chosen, set base frequencies\
`verbose`: 0, 1, or 2 determines level of output\
`outgroup`: set to species name. Default is no outgroup. \
`startmode`: set to smc or sim to perform SMC or simulate data

## The following settings apply to SMC option.
`datafile`: specify the name of the nexus file containing the raw sequence data\
`nparticles`: number of particles to use in the first-level SMC\
`particle_increase`: number of species particles per gene particle in second-level SMC\
`thin`: fraction of particles to preserve from first-level into second-level\
`save_every`: sets number of samples to save. ex. save_every = 100 saves 1 out of every 100 trees sampled\
`nthreads`: set to >1 for multithreading\
`theta_proposal_mean`: mean of inverse gamma proposal distribution when estimating `theta` mean to draw population thetas; if nothing provided, mean will be set to `theta`\
`theta_prior_mean`: mean of prior distribution when estimating `theta`. Can be the same as proposal or different.\
`lambda_prior_mean`: mean of exponential prior and proposal distribution when estimating `lambda`. If nothing is specified, lambda will be fixed (recommended for now).\
`fix_theta`: set to true to fix one theta for all populations. Default is false so a separate theta will be drawn for each population.\
`relative_rates`: relative substitution rate for each locus. Relative rates must average to 1.0. Ex. `relative_rates = 0.5, 1.5` for 2 loci\
`run_on_empty`: set to true to sample from the prior\
`save_memory`: set to true to throw away partials in likelihood calculations after they are used. This option may take longer but may use less memory.

### The following settings apply only if reading in gene newicks and only performing second-level SMC. Default is to perform first- and second-level SMC.
`gene_newicks`: set to true if specifying gene files to be read in\
`newick_path`: specify the path to a file containing gene newicks obtained from another program. Only second-level SMC will be performed. Gene newicks must be ultrametric and have no zero length branches.
`ngenes`: if setting `gene_newicks` to true, set the number of genes provided\

## The following settings apply to simulation option (`startmode = sim`)
`filename`: name of file to write simulated data to\
`nspecies`: number of species to simulate data for\
`ntaxaperspecies`: number of samples per species; ex. `2, 3, 5` specifies 2 samples from species A, 3 samples from species B, 5 samples from species C\
`fix_theta_for_simulations`: set to true to fix one theta for all populations. Default is true.\
`save_gene_trees_separately`: set to true to save gene newicks in separate files rather than all in one file. This makes it easier to read back in gene trees if using true gene trees to run a second-level only SMC.\
`simoccupancy`: probability that any given taxon will have data for any given locus; 1-simoccupancy is probability of all missing data for a taxon\
`simedgeratevar`: variance of lognormal relative rate distribution across edges in gene trees\
`simasrvshape`: shape of gamma among-site rate heterogeneity within a locus\
`simcomphet`: Dirichlet parameter governing compositional heterogeneity (default value results in compositional homogeneity)
