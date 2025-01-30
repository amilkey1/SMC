# sequential multispecies coalescent
This program uses a sequential Monte Carlo approach to sample the posterior distribution of trees under the multispecies coalescent model.

# conf tutorial
The file proj.conf is used to specify program settings:

`datafile`: specify the name of the nexus file containing the raw sequence data

`nparticles`: number of particles to use in the SMC

proposal: specify the prior-prior or prior-post proposal (default is prior-post if nothing is specified)

seed: specify random number seed
  - default is 1

theta: specify the population-scaled mutation rate
  - if estimate_theta is set to true, this value will be the starting value for theta

speciation rate: specify the speciation rate

subset: specify the gene partitions in the data file

  for example:
  
      subset = gene1:1-1000
      
      subset = gene2:1001-2000

      specifies that bases 1-1000 correspond to gene 1 and bases 1001-2000 correspond to gene 2
      
  if no subset is specified, the entire file will be read as one gene

model: specify the GTR-family substitution model (current valid options are JC and HKY)

kappa: if a model that includes the transtion / transversion rate ratio is specified under "model", specify a value for kappa
  - default is 1.0

base_frequencies: if a non-JC model is specified under "model", specify values for the base frequencies in alphabetical order
  
  for example:
  
    base_frequencies = 0.2, 0.2, 0.3, 0.3

    specifies that the frequency of base A = 0.2, C=0.2, G=0.3, 0.3

estimate_theta: set to true if the program should use the specified theta as a starting value and estimate theta as a parameter (default is false)

gene_newicks: specify the name of a .txt file containing gene newicks obtained from another program
 
  - gene newicks must be ultrametric and have no zero length branches

species_newicks: specify the name of a .txt file containing a species newick obtained from another program
 
  - the program will read in the species newick topology and draw branch lengths from the species tree prior
  
  - if both gene and species newicks are specified, the program will calculate the coalescent likelihood for the species tree given those gene trees and exit

start_from_species_tree_prior: specify the program is to start by sampling from the species tree prior
  - do not specify any newicks if using this option

start_from_gene_tree_prior: specify the program is to start by sampling from the gene tree prior
  - do not specify any newicks if using this option

niterations: specify the number of times the program should alternate between filtering the gene trees and species trees
  - 1 iteration = estimate species trees and then estimate gene trees (or vice versa, depending on whether the program is started with species trees or gene trees)
  - if starting from species trees, must specify >1 iteration

run_on_empty: set to true for the program to run without data, sampling only from the prior
  - default is false

species_particles_per_gene_particle: specify number of species particles to be created per 1 gene particle
  - for example:
      nparticles = 10
      species_particle_per_gene_particle = 20
    The gene trees will be filtered using 10 particles, and the species trees will be filtered using 10*20 =       200 particles.

outgroup: specify the name of the outgroup
  - if nothing is specified, no outgroup will be used

estimate_theta: set to true to estimate theta
  - default value is false

estimate_lambda: set to true to estimate lambda
  - default value is false

ntries_theta: specify the number of values of theta to try during theta estimation
  - default value is 50

ntries_lambda: specify the number of values of lambda to try during lambda estimation
  - default value is 50
