# sequential multispecies coalescent
This program uses a sequential Monte Carlo approach to sample the posterior distribution of trees under the multispecies coalescent model.

# conf tutorial
The file proj.conf is used to specify program settings:

datafile: specify the name of the nexus file containing the raw sequence data

nparticles: number of particles to use in the SMC

seed: random seed

theta: specify the population-scaled mutation rate

speciation rate: specify the speciation rate

subset: specify the gene partitions in the data file

  for example:
  
      subset = gene1:1-1000
      
      subset = gene2:1001-2000

      specifies that bases 1-1000 correspond to gene 1 and bases 1001-2000 correspond to gene 2
      
  if no subset is specified, the entire file will be read as one gene

model: specify the GTR-family substitution model (current valid options are JC and HKY)

kappa: if a model that includes the transtion / transversion rate ratio is specified under "model", specify a value for kappa

base_frequencies: if a non-JC model is specified under "model", specify values for the base frequencies in alphabetical order
  
  for example:
  
    base_frequencies = 0.2, 0.2, 0.3, 0.3

    specifies that the frequency of base A = 0.2, C=0.2, G=0.3, 0.3

estimate_theta: set to true if the program should use the specified theta as a starting value and estimate theta as a parameter (default is false)

estimate_speciation_rate: set to true if the program should use the specified speciation rate as a starting value and estimate the speciation rate as a parameter (default is false)

gene_newicks: specify the name of a .txt file containing gene newicks obtained from another program
 
  gene newicks must be ultrametric and have no zero length branches

species_newicks: specify the name of a .txt file containing a species newick obtained from another program
 
  - the program will read in the species newick topology and draw branch lengths from the species tree prior
  
  - a species newick OR a gene newick may be specified but not both
  
  - if no newick is specified, the program will begin by sampling from the species tree prior

niterations: specify the number of times the program should alternate between filtering the gene trees and species trees
 
  - 1 iteration = estimate species trees and then estimate gene trees (or vice versa, depending on whether the program is started with species trees or gene trees)
