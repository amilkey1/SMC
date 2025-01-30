# sequential multispecies coalescent
This program uses a sequential Monte Carlo approach to sample the posterior distribution of trees under the multispecies coalescent model.

# conf tutorial
The file proj.conf is used to specify program settings:

`seed`: specify random number seed. Default is 1. \
`theta`: specify the population-scaled mutation rate\
`lambda`: specify the speciation rate\
`subset`: specify the gene partitions in the data file. For example:
  `subset = gene1:1-1000`
  `subset = gene2:1001-2000`
specifies that bases 1-1000 correspond to gene 1 and bases 1001-2000 correspond to gene 2. If no subset is specified, the entire file will be read as one gene\

`nparticles`: number of particles to use in the first-level SMC\
`model`: the model to use in first-level SMC. Current options are JC or HKY\
`kappa`: if HKY is chosen, set kappa\
`base_frequencies`: if HKY is chosen, set base frequencies\
`run_on_empty`: set to true to sample from the prior\
`nthreads`: set to >1 for multithreading\
`verbose`: 0, 1, or 2 determines level of output\
`save_memory`: set to true to throw away partials in likelihood calculations after they are used\
`outgroup`: set to species name\
`startmode`: set to smc or sim to perform SMC or simulate data\

The following settings apply to SMC option:
`datafile`: specify the name of the nexus file containing the raw sequence data
`gene_newicks`: set to true if specifying gene files to be read in
`ngenes`: if setting `gene_newicks`, set the number of genes provided
`particle_increase`: number of species particles per gene particle in second-level SMC
`thin`: fraction of particles to preserve from first-level into second-level
`save_every`: sets number of samples to save. ex. save_every = 100 saves 1 out of every 100 trees sampled
`theta_proposal_mean`: mean of proposal distribution when estimating `theta`; if nothing provided, mean will be set to `theta`
`theta_prior_mean`: mean of prior distribution when estimating `theta`
`lambda_prior_mean`: mean of prior distribution when estimating `lambda`

specify the name of a .txt file containing gene newicks obtained from another program. Gene newicks must be ultrametric and have no zero length branches.

The following settings apply to simulation option:
`filename`: name of file to write simulated data to
`nspecies`: number of species to simulate data for
`ntaxaperspecies`: number of samples per species
`fix_theta_for_simulations`: set to true to fix one theta for all populations
