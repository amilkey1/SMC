import sys, os, re, subprocess as sub, shutil, numpy as np
from scipy.stats import randint,uniform,lognorm,describe
from math import log,exp,sqrt,pow
import math

# Settings you can change
method           = 'grid' # should be either 'uniform' or 'lognorm' or 'grid'
ntax           = [2,2,2,2,2] # number of taxa in each species

# These used only if method == 'uniform' or 'grid'
T_low            = 0.0       # smallest tree height (T) value 
T_high           = 1.0       # largest tree height (T) value
half_theta_low   = 0.0       # smallest theta/2 value
half_theta_high  = 0.5       # largest theta/2 value

# These used only if method == 'lornorm' 
Tmean            = 1.0       # mean tree height (T)
Tsd              = 0.7       # standard deviation of T
Rmean            = 0.2       # mean ratio of theta to T
Rsd              = 0.2       # standard deviation of theta/T ratios

nloci          = 10          # number of loci (conditionally independent given species tree)
seqlen         = 100       # number of sites in each gene
nreps          = 4          # number of simulation replicates
nparticles     = 100       # number of particles to use for SMC
simprogname    = 'smc'    # name of program used to simulate data (expected to be in $HOME/bin on cluster)
smcprogname    = 'smc'    # name of program used to perform SMC (expected to be in $HOME/bin on cluster)
beastprogname  = 'beast'     # name of program used to perform SMC (expected to be in $HOME/bin on cluster)
paupprogname = 'paup4a168_osx'    #name of PAUP*
smctreefname   = 'species_trees.trees' # name of species tree file for SMC
beasttreefname = 'species.trees'           # name of species tree file for BEAST
svdqtreefname = 'svd.tre'           # name of species tree file for BEAST
username       = 'aam21005'  # name of user on UConn HPC cluster
nodechoices    = [('general', 'epyc128'), ('priority','skylake')]
nodechoice     = 0          # 0-offset index into nodechoices
#partition      = 'general'   # specifies partition to use for HPC: either 'general' or 'priority'
#constraint     = 'epyc128'   # specifies constraint to use for HPC: e.g. 'skylake', 'epyc128', etc.
dirname        = 'g'         # name of directory created (script aborts if it already exists)
rnseed         = 123579      # overall pseudorandom number seed for everything except setting sim conf file
rnsimseed      = 123       # overall pseudorandom number seed for setting sim conf file
mcmciter       = 5000      # chain length for Beast MCMC
saveevery      = 1000         # MCMC storeevery modulus
preburnin      = 0        # MCMC burn in
storeevery     = 1000        # state storeevery modulus
screenevery    = 1000        # screen print modulus
genetreeevery  = 1000         # gene tree save modulus
spptreeevery   = 1000          # species tree save modulus (mcmciter/spptreeevery should equal nparticles

# Settings you can change but probably shouldn't
maxsimult   = None        # maximum number of jobs to run simultaneously (set to None if there is no maximum)

# Values obtained from settings
np.random.seed(seed=rnsimseed)

# Set up random variable for choosing seeds for each simulation
rnsimseeds = randint.rvs(1, 1000000, size=nreps)

# Values obtained from settings
np.random.seed(seed=rnseed)
nspecies = len(ntax)

# Set up random variable for choosing seeds for each replicate
rnseeds = randint.rvs(1, 1000000, size=nreps)

if method == 'lognorm':
    # Set up "standardized" scipy lognormal distribution for T
    Tsigsq = log(1. + pow(Tsd/Tmean,2.))
    Tsigma = sqrt(Tsigsq)
    Tmu = log(Tmean) - 0.5*Tsigsq
    Tshape = Tsigma
    Tscale = exp(Tmu)
    Trv = lognorm(Tshape, scale=Tscale)
    Tvect = Trv.rvs(size=nreps)
    Td = describe(Tvect, 0, 1)  # 0 = axis, 1 = df used in correcting variance)

    # Set up "standardized" scipy lognormal distribution for theta/T ratio
    Rsigsq = log(1. + pow(Rsd/Rmean,2.))
    Rsigma = sqrt(Rsigsq)
    Rmu = log(Rmean) - 0.5*Rsigsq
    Rshape = Rsigma
    Rscale = exp(Rmu)
    Rrv = lognorm(Rshape, scale=Rscale)
    Rvect = Rrv.rvs(size=nreps)
    Rd = describe(Rvect, 0, 1)  # 0 = axis, 1 = df used in correcting variance)

    thetas = [r*t for r,t in zip(Rvect, Tvect)]
    thetad = describe(thetas, 0, 1)

    nspp = len(ntax)
    phi = sum([1./k for k in range(2, nspp + 1)])
    lambdas = [phi/t for t in Tvect]
    lambdad = describe(lambdas, 0, 1)
elif method == 'uniform':
    thetas = [2.*q for q in uniform.rvs(loc=half_theta_low, scale=half_theta_high-half_theta_low, size=nreps)]
    thetad = describe(thetas, 0, 1)

    Tvect = uniform.rvs(loc=T_low, scale=T_high-T_low, size=nreps)
    Td = describe(Tvect, 0, 1)  # 0 = axis, 1 = df used in correcting variance)
    nspp = len(ntax)
    phi = sum([1./k for k in range(2, nspp + 1)])
    lambdas = [phi/t for t in Tvect]
    lambdad = describe(lambdas, 0, 1)
elif method == 'grid':
    nreps_to_use = math.isqrt(nreps)
#    if (nreps_to_use.is_integer() == False):
#        assert False, 'if choosing grid method, nreps must be a square of an int'
    thetas_placeholder = [half_theta_low + (x+1)*(half_theta_high-half_theta_low)/nreps_to_use for x in range(nreps_to_use)]
    thetas_placeholder = [i * 2.0 for i in thetas_placeholder]
    lambdas_placeholder = [T_low + (x+1)*(T_high-T_low)/nreps_to_use for x in range(nreps_to_use)]
    
    thetas = []
    Tvect = []
    
    # create grid
    for i in range(len(thetas_placeholder)):
        for x in range(len(thetas_placeholder)):
            thetas.append(thetas_placeholder[i])

    count = 0
    while (count < nreps_to_use):
        for i in range(len(lambdas_placeholder)):
            Tvect.append(lambdas_placeholder[i])
        count += 1
    
    thetad = describe(thetas, 0, 1)
    Td = describe(Tvect, 0, 1)
    nspp = len(ntax)
    phi = sum([1./k for k in range(2, nspp + 1)])
    lambdas = [phi/t for t in Tvect]
    lambdad = describe(lambdas, 0, 1)

else:
    assert False, 'method should be either "lognorm" or "uniform" or "grid" but you specified "%s"' % method

def inventName(k, lower_case):
    # If   0 <= k < 26, returns A, B, ..., Z,
    # If  26 <= k < 702, returns AA, AB, ..., ZZ,
    # If 702 <= k < 18278, returns AAA, AAB, ..., ZZZ, and so on.
    #
    # For example, k = 19009 yields ABCD:
    # ABCD 19009 = 26 + 26*26 + 26*26*26 + 0*26*26*26 + 1*26*26 + 2*26 + 3
    #              <------- base ------>   ^first       ^second   ^third ^fourth
    # base = (26^4 - 1)/25 - 1 = 18278
    #   26^1 + 26^2 + 26^3 = 26^0 + 26^1 + 26^2 + 26^3 - 1 = (q^n - 1)/(q - 1) - 1, where q = 26, n = 4
    #   n = 1 + floor(log(19009)/log(26))
    # fourth = ((19009 - 18278                           )/26^0) % 26 = 3
    # third  = ((19009 - 18278 - 3*26^0                  )/26^1) % 26 = 2
    # second = ((19009 - 18278 - 3*26^0 - 2*26^1         )/26^2) % 26 = 1
    # first  = ((19009 - 18278 - 3*26^0 - 2*26^1 - 1*26^2)/26^3) % 26 = 0

    # Find how long a name string must be
    logibase26 = 0
    if k > 0:
        logibase26 = log(k)/log(26)
    n = 1 + int(logibase26)
    letters = []
    base = (pow(26,n) - 1)/25.0 - 1
    cum = 0
    ordA = lower_case and ord('a') or ord('A')
    for i in range(n):
        ordi = int((k - base - cum)/pow(26,i) % 26)
        letters.append(chr(ordA + ordi))
        cum += ordi*pow(26,i)
    letters.reverse()
    invented_name = ''.join(letters)
    return invented_name

def writeNexusFile(fn, ntax, nchar, mask, taxa, sequences):
    if os.path.exists(fn):
        os.rename(fn, '%s.bak' % fn)
    longest = max([len(t) for t in taxa])
    taxonfmt = '  %%%ds' % longest
    f = open(fn, 'w')
    f.write('#nexus\n\n')
    f.write('begin data;\n')
    f.write('  dimensions ntax=%d nchar=%d;\n' % (ntax, nchar))
    f.write('  format datatype=dna gap=-;\n')
    f.write('  matrix\n')
    if mask is not None:
        f.write(taxonfmt % ' ')
        f.write('[%s]\n' % mask)
    for t in taxa:
        taxon_name = re.sub('\s+', '_', t)
        f.write(taxonfmt % taxon_name)
        f.write(' %s\n' % sequences[t])
    f.write('  ;\n')
    f.write('end;\n')
    f.close()

def createMainDir():
    if os.path.exists(dirname):
        sys.exit('dirname "%s" already exists; please move or rename it and try again' % dirname)
    else:
        os.mkdir(dirname)

def createRepDirName(rep_index):
    j = rep_index + 1
    repdirpath = 'rep%d' % j
    return repdirpath

def createRepDirPath(rep_index):
    j = rep_index + 1
    repdirname = 'rep%d' % j
    repdirpath = os.path.join(dirname, repdirname)
    return repdirpath

def createSimDirPath(rep_index):
    repdirpath = createRepDirPath(rep_index)
    simdirpath = os.path.join(repdirpath, 'sim')
    return simdirpath

def createSMCDirPath(rep_index):
    repdirpath = createRepDirPath(rep_index)
    smcdirpath = os.path.join(repdirpath, 'smc')
    return smcdirpath

def createBeastDirPath(rep_index):
    repdirpath = createRepDirPath(rep_index)
    beastdirpath = os.path.join(repdirpath, 'beast')
    return beastdirpath

def createRepDir(rep_index):
    repdirpath = createRepDirPath(rep_index)
    os.mkdir(repdirpath)

def createSimDir(rep_index):
    simdirpath = createSimDirPath(rep_index)
    os.mkdir(simdirpath)

def createSMCDir(rep_index):
    smcdirpath = createSMCDirPath(rep_index)
    os.mkdir(smcdirpath)

def createBeastDir(rep_index):
    beastdirpath = createBeastDirPath(rep_index)
    os.mkdir(beastdirpath)
    
def createGetSVDQHeight():
	svdqfn = os.path.join(dirname, 'get-svdq-lamdas.py')
		
	s = "i = 1\n"
	s += 'new_f = open("svdq_lambdas.txt", "x")\n'
	s += "for rep in range(%12d):\n" % nreps
	s += "	name = 'rep' + str(i)\n"
	s += '	count = 0\n'
	s += '	found = False\n'
	s += '	numerator = 0.0\n'
	s += "	name = 'rep' + str(i) + '/sim' + '/svdout.txt'\n"
	s += "	for a in range(%d+1):\n" % nspecies
	s += "		if a > 1:\n"
	s += "			numerator += 1. / a\n"
	s += '	for line in reversed(list(open(name))):\n'
	s += '		if line.startswith("Sum"):\n'
	s += '			count += 1\n'
	s += '			found = True\n'
	s += '		elif found:\n'
	s += '			count += 1\n'
	s += '		if count == 3:\n'
	s += '			height = float(line[39:49])\n'
	s += '			lamda = numerator / height\n'
	s += '			new_f.write(str(lamda)+"\\n")\n'
	s += '			break\n'
	s += "	i += 1\n"
	svdqf = open(svdqfn, 'w')
	svdqf.write(s)
	svdqf.close()


def createSimConf(rep_index):
    theta = thetas[rep_index]
    lamda = lambdas[rep_index]
    simdirpath = createSimDirPath(rep_index)
    fn = os.path.join(simdirpath, 'proj.conf')

    s  = ''
    s += 'filename  = sim.nex\n'
    s += 'startmode = sim\n'
    s += 'seed    = %d\n' % rnsimseeds[rep_index]
    s += '\n'
    cum = 0
    for g in range(nloci):
        locus = g + 1
        s += 'subset = locus%d[nucleotide]:%d-%d\n' % (locus, cum + 1, cum + seqlen)
        cum += seqlen
    s += '\n'
    s += 'theta  = %.2f\n' % theta
    s += 'lambda = %.2f\n' % lamda
    s += '\n'
    s += 'nspecies = %d\n' % nspecies
    s += 'ntaxaperspecies ='
    for spp in range(nspecies):
        s += str(ntax[spp])
        if spp != nspecies-1:
            s += ','
    s += '\n'
    s += 'verbose = 0\n'
    outf = open(fn, 'w')
    outf.write(s)
    outf.close()
    
def createSMCConf(rep_index):
    theta = thetas[rep_index]
    lamda = lambdas[rep_index]
    smcdirpath = createSMCDirPath(rep_index)
    smcconffn = os.path.join(smcdirpath, 'proj.conf')

    s  = ''
    s += 'startmode = smc\n'
    s += 'datafile  = ../sim/sim.nex\n'
    s += 'seed    = %d\n' % rnseeds[rep_index]
    s += '\n'
    s += '\n'
    cum = 0
    for g in range(nloci):
        locus = g + 1
        s += 'subset = locus%d[nucleotide]:%d-%d\n' % (locus, cum + 1, cum + seqlen)
        cum += seqlen
    s += 'theta  = %.2f\n' % theta
    s += 'lambda = %.2f\n' % lamda
    s += '\n'
    s += '\n'
    s += 'theta_proposal_mean = 1.0\n'
    s += 'theta_prior_mean = 1.0\n'
    s += '\n'
    s += '\n'
    s += 'nparticles = %d\n' % nparticles
    s += 'nthreads = 36\n'
    s += '\n'
    s += 'verbose = 1\n'
    s += 'run_on_empty = false\n'
    s += 'particle_increase = 5\n'
    s += 'thin=0.01\n'
    s += 'save_every = 1\n'
    s += 'save_gene_trees = false\n'

    smcconff = open(smcconffn, 'w')
    smcconff.write(s)
    smcconff.close()

def createBeastXML(rep_index):
    beastdirpath = createBeastDirPath(rep_index)
    beastxmlfn = os.path.join(beastdirpath, 'starbeast.xml')

    s  = ''
    s += '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
    s += '<beast beautitemplate=\'StarBeast3\' beautistatus=\'noAutoSetClockRate|noAutoUpdateFixMeanSubstRate\' namespace="beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.base.evolution.alignment:beast.pkgmgmt:beast.base.core:beast.base.inference:beast.base.evolution.tree.coalescent:beast.pkgmgmt:beast.base.core:beast.base.inference.util:beast.evolution.nuc:beast.base.evolution.operator:beast.base.inference.operator:beast.base.evolution.sitemodel:beast.base.evolution.substitutionmodel:beast.base.evolution.likelihood" required="BEAST.base v2.7.6:starbeast3 v1.1.8:BEASTLabs v2.0.2:ORC v1.1.2" version="2.7">\n'
    s += '\n'
    for gene in range(1,nloci+1):
        s += '    <data\n'
        s += 'id="gene%d"\n' % gene
        s += 'spec="Alignment"\n'
        s += 'name="alignment">\n'
        species = 0
        taxon = 0
        for n in ntax:
            species += 1
            species_name = inventName(species - 1, False)
            for t in range(n):
                taxon_name = inventName(taxon, True)
                taxon += 1
                s += '        <sequence id="seq_%s^%s%d" spec="Sequence" taxon="%s^%s" totalcount="4" value="__SEQ_%d_%s^%s__"/>\n' % (taxon_name, species_name, gene, taxon_name, species_name, gene, taxon_name, species_name)
        s += '    </data>\n\n'
    s += '    <map name="Uniform" >beast.base.inference.distribution.Uniform</map>\n\n'
    s += '    <map name="Exponential" >beast.base.inference.distribution.Exponential</map>\n\n'
    s += '    <map name="LogNormal" >beast.base.inference.distribution.LogNormalDistributionModel</map>\n\n'
    s += '    <map name="Normal" >beast.base.inference.distribution.Normal</map>\n\n'
    s += '    <map name="Beta" >beast.base.inference.distribution.Beta</map>\n\n'
    s += '    <map name="Gamma" >beast.base.inference.distribution.Gamma</map>\n\n'
    s += '    <map name="LaplaceDistribution" >beast.base.inference.distribution.LaplaceDistribution</map>\n\n'
    s += '    <map name="prior" >beast.base.inference.distribution.Prior</map>\n\n'
    s += '    <map name="InverseGamma" >beast.base.inference.distribution.InverseGamma</map>\n\n'
    s += '    <map name="OneOnX" >beast.base.inference.distribution.OneOnX</map>\n\n'
    s += '    <run id="mcmc" spec="MCMC" chainLength="%d" preBurnin="%d" storeEvery="%d">\n' % (mcmciter, preburnin, saveevery)
    s += '        <state id="state" spec="State" storeEvery="%d">\n' % storeevery
    s += '            <stateNode id="Tree.t:Species" spec="starbeast3.tree.SpeciesTree">\n'
    s += '                <taxonset id="taxonsuperset" spec="starbeast3.tree.StarBeast3TaxonSet">\n'
    taxon = 1
    for spp in range(len(ntax)):
        species = spp + 1
        species_name = inventName(species - 1, False)
        s += '                    <taxon id="%s" spec="TaxonSet">\n' % species_name
        for t in range(ntax[spp]):
            taxon_name = inventName(taxon - 1, True)
            s += '                        <taxon id="%s^%s" spec="Taxon"/>\n' % (taxon_name, species_name)
            taxon += 1
        s += '                    </taxon>\n'
    s += '                </taxonset>\n'
    s += '            </stateNode>\n'
    s += '            <parameter id="popSize" spec="parameter.RealParameter" lower="0.0" name="stateNode">1.0</parameter>\n'
    s += '            <parameter id="popMean" spec="parameter.RealParameter" lower="0.0" name="stateNode">1.0</parameter>\n'
    for gene in range(1,nloci+1):
        s += '            <tree id="Tree.t:gene%d" spec="beast.base.evolution.tree.Tree" name="stateNode">\n' % gene
        s += '                <taxonset id="TaxonSet.gene%d" spec="TaxonSet">\n' % gene
        s += '                    <alignment idref="gene%d"/>\n' % gene
        s += '                </taxonset>\n'
        s += '            </tree>\n'
    s += '        </state>\n'

    s += '        <init id="SBI" spec="starbeast3.core.StarBeastStartState" estimate="false" popMean="@popMean" speciesTree="@Tree.t:Species">\n'
    for gene in range(1,nloci+1):
        s += '            <gene idref="Tree.t:gene%d"/>\n' % gene
    s += '                        <parameter id="speciationRate.t:Species" spec="parameter.RealParameter" estimate="false" lower="SVDQLAMBDA" name="birthRate" upper="SVDQLAMBDA">SVDQLAMBDA</parameter>\n'
    s += '            <speciesTreePrior id="SpeciesTreePopSize.Species" spec="starbeast3.evolution.speciation.SpeciesTreePrior" bottomPopSize="@popSize" gammaParameter="@popMean" taxonset="@taxonsuperset" tree="@Tree.t:Species">\n'
    s += '                <populationModel id="speciesTreePopulationModel" spec="starbeast3.evolution.speciation.ConstantPopulations" populationSizes="@popSize" speciesTree="@Tree.t:Species"/>\n'
    s += '                <treePrior id="YuleModel.t:Species" spec="beast.base.evolution.speciation.YuleModel" birthDiffRate="@speciationRate.t:Species" tree="@Tree.t:Species"/>\n'
    s += '            </speciesTreePrior>\n'
    s += '            <sharedRateModel id="branchRatesModel.Species" spec="starbeast3.evolution.branchratemodel.SharedSpeciesClockModel">\n'
    s += '                <branchRateModel id="strictClockModel.Species" spec="starbeast3.evolution.branchratemodel.StrictClockModelSB3" tree="@Tree.t:Species">\n'
    s += '                    <parameter id="SpeciesTreeStrictClockRate" spec="parameter.RealParameter" estimate="false" lower="0.0" name="clock.rate">1.0</parameter>\n'
    s += '                </branchRateModel>\n'
    s += '            </sharedRateModel>\n'
    s += '        </init>\n'

    s += '        <distribution id="posterior" spec="CompoundDistribution">\n'
    s += '            <distribution id="speciescoalescent" spec="CompoundDistribution">\n'
    for gene in range(1,nloci+1):
        s += '                <distribution id="treePrior.t:gene%d" spec="starbeast3.evolution.speciation.GeneTreeForSpeciesTreeDistribution" populationModel="@speciesTreePopulationModel" speciesTree="@Tree.t:Species" speciesTreePrior="@SpeciesTreePopSize.Species" tree="@Tree.t:gene%d"/>\n' % (gene,gene)
    s += '            </distribution>\n'

    s += '            <distribution idref="SpeciesTreePopSize.Species"/>\n'

    s += '            <distribution id="prior" spec="CompoundDistribution">\n'
    s += '                <distribution idref="YuleModel.t:Species"/>\n'
    s += '                <prior id="popMean.prior" name="distribution" x="@popMean">\n'
    s += '                    <Exponential id="Exponential.11" name="distr">\n'
    s += '                        <parameter id="RealParameter.0" spec="parameter.RealParameter" estimate="false" name="mean">4.0</parameter>\n'
    s += '                    </Exponential>\n'
    s += '                </prior>\n'
    s += '            </distribution>\n'

    s += '            <distribution id="vectorPrior" spec="CompoundDistribution">\n'
    s += '                <prior id="constPopSizesPrior.Species" name="distribution" x="@popSize">\n'
    s += '                    <InverseGamma id="popPriorDistr.InverseGamma" beta="@popMean" name="distr">\n'
    s += '                        <alpha id="Function$Constant.0" spec="Function$Constant" value="2.0"/>\n'
    s += '                    </InverseGamma>\n'
    s += '                </prior>\n'
    s += '            </distribution>\n'

    s += '            <distribution id="likelihood" spec="CompoundDistribution" useThreads="true">\n'
    for gene in range(1,nloci+1):
        s += '                <distribution id="treeLikelihood.gene%d" spec="TreeLikelihood" data="@gene%d" tree="@Tree.t:gene%d">\n' % (gene, gene, gene)
        s += '                    <siteModel id="SiteModel.s:gene%d" spec="SiteModel">\n' % gene
        s += '                        <parameter id="mutationRate.s:gene%d" spec="parameter.RealParameter" estimate="false" name="mutationRate">1.0</parameter>\n' % gene
        s += '                        <parameter id="gammaShape.s:gene%d" spec="parameter.RealParameter" estimate="false" lower="0.0" name="shape">1.0</parameter>\n' % gene
        s += '                        <parameter id="proportionInvariant.s:gene%d" spec="parameter.RealParameter" estimate="false" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>\n' % gene
        s += '                        <substModel id="JC69.s:gene%d" spec="JukesCantor"/>\n' % gene
        s += '                    </siteModel>\n'
        s += '                    <branchRateModel id="GeneTreeClock.c:gene%d" spec="starbeast3.evolution.branchratemodel.StarBeast3Clock" geneTree="@treePrior.t:gene%d" sharedRateModel="@branchRatesModel.Species" tree="@Tree.t:gene%d">\n' % (gene, gene, gene)
        s += '                        <parameter id="clockRate.c:gene%d" spec="parameter.RealParameter" estimate="false" lower="0.0" name="clock.rate">1.0</parameter>\n' % gene
        s += '                    </branchRateModel>\n'
        s += '                </distribution>\n' # gene
    s += '            </distribution>\n' # likelihood
    s += '        </distribution>\n' # posterior

    s += '        <operator id="Reheight.t:Species" spec="starbeast3.operators.NodeReheight2" taxonset="@taxonsuperset" tree="@Tree.t:Species" weight="30.0">\n'
    for gene in range(1,nloci+1):
        s += '            <gene idref="treePrior.t:gene%d"/>\n' % gene
    s += '        </operator>\n'

    s += '        <operator id="CoordinatedExponential.t:Species" spec="starbeast3.operators.CoordinatedExponential" speciesTree="@Tree.t:Species" weight="15.0">\n'
    for gene in range(1,nloci+1):
        s += '            <gene idref="treePrior.t:gene%d"/>\n' % gene
    s += '        </operator>\n'

    s += '        <operator id="CoordinatedUniform.t:Species" spec="starbeast3.operators.CoordinatedUniform" speciesTree="@Tree.t:Species" weight="30.0">\n'
    for gene in range(1,nloci+1):
        s += '            <gene idref="treePrior.t:gene%d"/>\n' % gene
    s += '        </operator>\n'

    s += '        <operator id="TreeRootScaler.t:Species" spec="kernel.BactrianScaleOperator" rootOnly="true" scaleFactor="0.7" tree="@Tree.t:Species" upper="10.0" weight="3.0"/>\n'
    s += '        <operator id="BactrianNodeOperator.t:Species" spec="kernel.BactrianNodeOperator" tree="@Tree.t:Species" weight="3.0"/>\n'
    s += '        <operator id="AdaptableTopologyOperator.lengths.Species" spec="AdaptableOperatorSampler" weight="100.0">\n'
    s += '            <tree idref="Tree.t:Species"/>\n'
    s += '            <operator idref="BactrianNodeOperator.t:Species"/>\n'
    s += '            <operator id="TreeScaler.t:Species" spec="kernel.BactrianScaleOperator" scaleFactor="0.01" tree="@Tree.t:Species" upper="10.0" weight="1.0"/>\n'
    s += '            <operator idref="CoordinatedUniform.t:Species"/>\n'
    s += '            <operator idref="CoordinatedExponential.t:Species"/>\n'
    s += '            <operator id="updown.all" spec="operator.kernel.BactrianUpDownOperator" scaleFactor="0.75" weight="1.0">\n'
    s += '                <down idref="Tree.t:Species"/>\n'
    s += '                <down idref="popSize"/>\n'
    s += '                <down idref="popMean"/>\n'
    for gene in range(1,nloci+1):
        s += '                <down idref="Tree.t:gene%d"/>\n' % gene
    s += '            </operator>\n'
    s += '        </operator>\n'

    s += '        <operator id="PopSizeGibbsSampler.Species" spec="starbeast3.operators.PopSizeGibbsSampler" gammaprior="@popPriorDistr.InverseGamma" popSizes="@popSize" weight="50.0">\n'
    for gene in range(1,nloci+1):
        s += '            <gene idref="treePrior.t:gene%d"/>\n' % gene
    s += '        </operator>\n'

    s += '        <operator id="AdaptableOperatorSampler.popmean:Species" spec="AdaptableOperatorSampler" weight="5.0">\n'
    s += '            <parameter idref="popMean"/>\n'
    s += '            <operator id="Scale.popmean" spec="kernel.BactrianScaleOperator" parameter="@popMean" upper="10.0" weight="1.0"/>\n'
    s += '            <operator idref="updown.all"/>\n'
    s += '            <operator id="SampleFromPriorOperator.popmean" spec="orc.operators.SampleFromPriorOperator" parameter="@popMean" prior2="@popMean.prior" weight="1.0"/>\n'
    s += '        </operator>\n'

    s += '        <operator id="ParallelMCMCTreeOperator" spec="starbeast3.operators.ParallelMCMCTreeOperator" chainCoverage="1.0" learning="false" nregression="50" otherState="@state" runtime="1000.0" speciesTree="@Tree.t:Species" targetCPU="0.0" weight="1.0">\n'
    for gene in range(1,nloci+1):
        s += '            <distribution id="ParallelMCMCTreeOperatorLikelihood.gene%d" spec="starbeast3.operators.ParallelMCMCTreeOperatorTreeDistribution" geneprior="@treePrior.t:gene%d" tree="@Tree.t:gene%d" treelikelihood="@treeLikelihood.gene%d"/>\n' % (gene, gene, gene, gene)
    s += '            <schedule id="operatorSchedule" spec="starbeast3.core.OperatorScheduleRecalculator">\n'
    s += '                <subschedule id="operatorSubschedule" spec="OperatorSchedule" operatorPattern="^ParallelMCMCTreeOperator$" weight="1.0" weightIsPercentage="true"/>\n'
    s += '            </schedule>\n'
    s += '        </operator>\n'

    s += '        <logger id="tracelog" spec="Logger" fileName="starbeast3.log" logEvery="10000" model="@posterior" sort="smart">\n'
    s += '            <log idref="posterior"/>\n'
    s += '            <log idref="likelihood"/>\n'
    s += '            <log idref="prior"/>\n'
    s += '            <log idref="vectorPrior"/>\n'
    s += '            <log idref="speciescoalescent"/>\n'
    s += '            <log idref="speciationRate.t:Species"/>\n'
    s += '            <log id="TreeStat.Species" spec="beast.base.evolution.tree.TreeStatLogger" tree="@Tree.t:Species"/>\n'
    s += '            <log idref="YuleModel.t:Species"/>\n'
    s += '            <log idref="popMean"/>\n'
    s += '            <log idref="popSize"/>\n'

    clustertree_id = 2
    for gene in range(1,nloci+1):
        s += '            <log idref="treeLikelihood.gene%d"/>\n' % gene
        s += '            <log idref="treePrior.t:gene%d"/>\n' % gene
        s += '            <log id="TreeStat.t:gene%d" spec="beast.base.evolution.tree.TreeStatLogger" tree="@Tree.t:gene%d"/>\n' % (gene, gene)
        s += '            <log id="TreeDistanceNJ.t:gene%d" spec="beastlabs.evolution.tree.TreeDistanceLogger" tree="@Tree.t:gene%d">\n' % (gene, gene)
        s += '                <ref id="ClusterTree.%d" spec="beast.base.evolution.tree.ClusterTree" clusterType="neighborjoining" taxa="@gene%d"/>\n' % (clustertree_id, gene)
        s += '            </log>\n'
        s += '            <log id="TreeDistanceUPGMA.t:gene%d" spec="beastlabs.evolution.tree.TreeDistanceLogger" tree="@Tree.t:gene%d">\n' % (gene, gene)
        s += '                <ref id="ClusterTree.%d" spec="beast.base.evolution.tree.ClusterTree" clusterType="upgma" taxa="@gene%d"/>\n' % (clustertree_id + 1, gene)
        s += '            </log>\n'
        clustertree_id += 2
    s += '        </logger>\n'

    s += '        <logger id="speciesTreeLogger" spec="Logger" fileName="species.trees" logEvery="%d" mode="tree">\n' % spptreeevery
    s += '            <log id="SpeciesTreeLoggerX" spec="starbeast3.core.SpeciesTreeLogger" popSize="@popSize" speciesTreePrior="@SpeciesTreePopSize.Species" tree="@Tree.t:Species">\n'
    s += '                <treetop id="treeTopFinder" spec="beast.base.evolution.speciation.TreeTopFinder">\n'
    for gene in range(1,nloci+1):
        s += '                    <tree idref="Tree.t:gene%d"/>\n' % gene
    s += '                </treetop>\n'
    s += '            </log>\n'
    s += '        </logger>\n'
    s += '        <logger id="screenlog" spec="Logger" logEvery="%d">\n' % screenevery
    s += '            <log idref="posterior"/>\n'
    s += '            <log id="ESS.0" spec="util.ESS" arg="@posterior"/>\n'
    s += '            <log idref="likelihood"/>\n'
    s += '            <log idref="prior"/>\n'
    s += '        </logger>\n'
    for gene in range(1,nloci+1):
        s += '        <logger id="treelog.t:gene%d" spec="Logger" fileName="$(tree).trees" logEvery="%d" mode="tree">\n' % (gene, genetreeevery)
        s += '            <log id="TreeWithMetaDataLogger.t:gene%d" spec="beast.base.evolution.TreeWithMetaDataLogger" tree="@Tree.t:gene%d"/>\n' % (gene, gene)
        s += '        </logger>\n'
    s += '        <operatorschedule idref="operatorSchedule"/>\n'
    s += '    </run>\n'
    s += '\n'
    s += '</beast>\n'

    beastxmlf = open(beastxmlfn, 'w')
    beastxmlf.write(s)
    beastxmlf.close()

def createREADME():
    readmefn = os.path.join(dirname, 'README')
    
    readme  = 'Summary of parameters used in simulations\n'
    readme += '-----------------------------------------\n'
    if method == 'lognorm':
        readme += 'T (tree height):\n'
        readme += '  nobs = %d\n' % Td.nobs
        readme += '  mean = %.5f (expected value = %.5f)\n' % (Td.mean, Tmean)
        readme += '  s.d. = %.5f (expected value = %.5f)\n' % (sqrt(Td.variance),Tsd)
        readme += '  min  = %.5f\n' % (min(Tvect),)
        readme += '  max  = %.5f\n' % (max(Tvect),)
        readme += '\n'
        readme += 'R (theta/T ratio):\n'
        readme += '  nobs = %d\n' % Rd.nobs
        readme += '  mean = %.5f (expected value = %.5f)\n' % (Rd.mean, Rmean)
        readme += '  s.d. = %.5f (expected value = %.5f)\n' % (sqrt(Rd.variance),Rsd)
        readme += '  min  = %.5f\n' % (min(Rvect),)
        readme += '  max  = %.5f\n' % (max(Rvect),)
        readme += '\n'
        readme += 'theta:\n'
        readme += '  nobs = %d\n' % thetad.nobs
        readme += '  mean = %.5f\n' % thetad.mean
        readme += '  s.d. = %.5f\n' % (sqrt(thetad.variance),)
        readme += '  min  = %.5f\n' % (min(thetas),)
        readme += '  max  = %.5f\n' % (max(thetas),)
        readme += '\n'
        readme += 'lambda:\n'
        readme += '  nobs = %d\n' % lambdad.nobs
        readme += '  mean = %.5f\n' % lambdad.mean
        readme += '  s.d. = %.5f\n' % (sqrt(lambdad.variance),)
        readme += '  min  = %.5f\n' % (min(lambdas),)
        readme += '  max  = %.5f\n' % (max(lambdas),)
        readme += '\n'
    elif method == 'uniform':
        readme += 'theta:\n'
        readme += '  nobs = %d\n' % thetad.nobs
        readme += '  mean = %.5f\n' % thetad.mean
        readme += '  s.d. = %.5f\n' % (sqrt(thetad.variance),)
        readme += '  min  = %.5f\n' % (min(thetas),)
        readme += '  max  = %.5f\n' % (max(thetas),)
        readme += '\n'
        readme += 'lambda:\n'
        readme += '  nobs = %d\n' % lambdad.nobs
        readme += '  mean = %.5f\n' % lambdad.mean
        readme += '  s.d. = %.5f\n' % (sqrt(lambdad.variance),)
        readme += '  min  = %.5f\n' % (min(lambdas),)
        readme += '  max  = %.5f\n' % (max(lambdas),)
        readme += '\n'
    elif method == 'grid':
        readme += 'theta:\n'
        readme += '  nobs = %d\n' % thetad.nobs
        readme += '  mean = %.5f\n' % thetad.mean
        readme += '  s.d. = %.5f\n' % (sqrt(thetad.variance),)
        readme += '  min  = %.5f\n' % (min(thetas),)
        readme += '  max  = %.5f\n' % (max(thetas),)
        readme += '\n'
        readme += 'lambda:\n'
        readme += '  nobs = %d\n' % lambdad.nobs
        readme += '  mean = %.5f\n' % lambdad.mean
        readme += '  s.d. = %.5f\n' % (sqrt(lambdad.variance),)
        readme += '  min  = %.5f\n' % (min(lambdas),)
        readme += '  max  = %.5f\n' % (max(lambdas),)
        readme += '\n'
    else:
        assert False, 'method should be either "lognorm" or "uniform" but you specified "%s"' % method

    readme += 'Installing SMC on the remote cluster\n'
    readme += '------------------------------------\n'
    readme += '\n'
    readme += '1. Place the smc executable in the $HOME/bin directory and ensure that this\n'
    readme += '   directory is on your PATH by editing your ".bashrc" file and adding this\n'
    readme += '   line at the bottom:\n'
    readme += '\n'
    readme += 'export PATH="$HOME/bin:$PATH"\n'
    readme += '\n'
    readme += '2. Ensure that the name of your executable is stored in the variable\n'
    readme += '   "smcprogname" in the "deploy.py" script.\n'
    readme += '\n'
    
    readme += 'Installing BEAST on the remote cluster\n'
    readme += '--------------------------------------\n'
    readme += '\n'
    readme += '1. Download the BEAST2 executable from https://www.beast2.org (e.g. Linux x86 version)\n'
    readme += '\n'
    readme += '2. Unpack and ensure that the "beast" run script is available at $HOME/beast/bin\n'
    readme += '\n'
    readme += '3. Modify the "beast" run script by changing the last line from (note ellipses:\n'
    readme == '   the beginning of the line is shown):\n'
    readme += '\n'
    readme += '"$JAVA" -Dlauncher.wait.for.exit=true...\n'
    readme += '\n'
    readme += '   to\n'
    readme += '\n'
    readme += '"$JAVA" -Dbeast.user.package.dir=/home/pol02003/beast-packages -Dlauncher.wait.for.exit=true...\n'
    readme += '\n'
    readme += '   where "/home/pol02003" should be replaced with the path to *your* home directory\n'
    readme += '\n'
    readme += '4. Now create the "beast-packages" directory:\n'
    readme += '\n'
    readme += 'mkdir $HOME/beast-packages\n'
    readme += '\n'
    readme += '   and install the starbeast3 package:\n'
    readme += '\n'
    readme += 'cd $HOME/beast/bin\n'
    readme += './packagemanager -dir $HOME/beast-packages -add starbeast3\n'
    readme += '\n'

    readme += 'Creating a simulation experiment on the remote cluster\n'
    readme += '-----------------------------------------------------\n'
    readme += '1. Remove the directory specified in the variable "dirname" (e.g. "g")\n'
    readme += '\n'
    readme += 'rm -rf g\n'
    readme += '\n'
    readme += '2. Run the deploy script:\n'
    readme += '\n'
    readme += 'python3 deploy.py\n'
    readme += '\n'
    readme += '4. Navigate into the "g" directory:\n'
    readme += '\n'
    readme += 'cd g\n'
    readme += '\n'
    
    readme += 'Running the simulation experiment on the remote cluster\n'
    readme += '-------------------------------------------------------\n'
    readme += '1. Simulate data and run SVD Quartets:\n'
    readme += '\n'
    readme += '. simulate.sh\n'
    readme += '\n'
    readme += '2. Copy data from sim.nex to starbeast.xml for each replicate:\n'
    readme += '\n'
    readme += 'python3 copydata.py\n'
    readme += '\n'
    readme += '3. Get SVD Quartets heights of each simulation and replace the true lambda in the smc conf file and beast xml file:\n'
    readme += '\n'
    readme += 'python3 get-svdq-lamdas.py\n'
    readme += 'python3 rreplace.py\n'
    readme += '\n'
    readme += '4. Start SMC runs:\n'
    readme += '\n'
    readme += 'sbatch smcslurm.sh\n'
    readme += '\n'
    readme += '5. Start BEAST runs:\n'
    readme += '\n'
    readme += 'sbatch beastslurm.sh\n'
    readme += '\n'
    readme += '6. Summarize results:\n'
    readme += '\n'
    
    readme += 'Summarizing results on the remote cluster\n'
    readme += '-----------------------------------------\n'
    readme += '1. Use td ("treedist") to compute KF  distances to true species tree for SMC and BEAST runs\n'
    readme += '\n'
    readme += '. smctd.sh\n'
    readme += '. beasttd.sh\n'
    readme += 'python3 crunch.py\n'
    readme += '\n'
    readme += '2. Carry out repeated-measures ANOVA on remote cluster\n'
    readme += '\n'
    readme += 'python3 anova-means.py  # ignores variation within runs, using the run mean as the data point\n'
    readme += 'python3 anova-full.py   # takes account of variation within runs, every posterior sample is used\n'
    readme += '\n'
    readme += '3. Get number of deep coalescences for each simulation\n'
    readme += '\n'
    readme += 'python3 get-deep-coal.py\n'
    readme += '4. Get time for each program\n'
    readme += '\n'
    readme += 'python3 get-times.py\n'
    readme += '\n'
    
    readme += 'Plotting results on local computer\n'
    readme += '-----------------------------------------\n'
    readme += 'Transfer the following files to local computer: \n'
    readme += '	rf-summary.txt\n'
    readme += '	kf-summary.txt\n'
    readme += '	simcond.R\n'
    readme += '	deep_coal.txt\n'
    readme += 'Run the following scripts to create plots:\n'
    readme += '	Rscript simcond.R\n'
    readme += '	Rscript rgl.R\n'
    readme += '	python3 plotdistances.py\n'
    readme += '\n'
    
    readmef = open(readmefn, 'w')
    readmef.write(readme)
    readmef.close()

def createRplot():
    plotfn = os.path.join(dirname, 'simcond.R')

    plotstuff  = 'cwd = system(\'cd "$( dirname "$0" )" && pwd\', intern = TRUE)\n'
    plotstuff += 'setwd(cwd)\n'
    plotstuff += 'pdf("simcond.pdf")\n'

    if method == 'lognorm':
        Tstr = ['%g' % t for t in Tvect]
        plotstuff += 'T = c(%s)\n' % ','.join(Tstr)

        Rstr = ['%g' % r for r in Rvect]
        plotstuff += 'R = c(%s)\n' % ','.join(Rstr)

        thetastr = ['%g' % q for q in thetas]
        plotstuff += 'theta = c(%s)\n' % ','.join(thetastr)

        lambdastr = ['%g' % l for l in lambdas]
        plotstuff += 'lambda = c(%s)\n' % ','.join(lambdastr)

        plotstuff += 'plot(theta, lambda, type="p", pch=19, main="Simulation conditions", xlab="theta", ylab="lambda")\n'
    elif method == 'uniform':
        Tstr = ['%g' % t for t in Tvect]
        plotstuff += 'T = c(%s)\n' % ','.join(Tstr)

        thetastr = ['%g' % q for q in thetas]
        plotstuff += 'theta = c(%s)\n' % ','.join(thetastr)

        plotstuff += 'plot(theta/2, T, type="p", pch=19, main="Simulation conditions", xlab="theta/2", ylab="T")\n'
    elif method == 'grid':
        Tstr = ['%g' % t for t in Tvect]
        plotstuff += 'T = c(%s)\n' % ','.join(Tstr)

        thetastr = ['%g' % q for q in thetas]
        plotstuff += 'theta = c(%s)\n' % ','.join(thetastr)

        plotstuff += 'plot(theta/2, T, type="p", pch=19, main="Simulation conditions", xlab="theta/2", ylab="T")\n'
    else:
        assert False, 'method should be either "lognorm" or "uniform" but you specified "%s"' % method

    plotstuff += 'dev.off()\n'

    plotstuff += 'rf <- read.table(file="rf-summary.txt")\n'
    plotstuff += 'rf_smc <- rf$V3\n'
    plotstuff += 'rf_beast <- rf$V6\n'
    plotstuff += '\n'
    plotstuff += 'kf <- read.table(file="kf-summary.txt")\n'
    plotstuff += 'kf_smc <- kf$V3\n'
    plotstuff += 'kf_beast <- kf$V6\n'
    plotstuff += '\n'
    plotstuff += '# make more plots\n'
    plotstuff += 'library(ggplot2)\n'
    plotstuff += 'df <- data.frame(T, theta)\n'
    plotstuff += '\n'

    plotstuff += 'max_kf_beast <- max(kf_beast)\n'
    plotstuff += 'max_kf_smc <- max(kf_smc)\n'
    plotstuff += 'max_kf <- max(c(kf_beast, kf_smc))\n'
    plotstuff += '\n'
    plotstuff += 'max_rf_beast <- max(rf_beast)\n'
    plotstuff += 'max_rf_smc <- max(rf_smc)\n'
    plotstuff += 'max_rf <- max(c(rf_beast, rf_smc))\n'
    plotstuff += '\n'
    plotstuff += 'myPalette<-colorRampPalette(c("purple","yellow"))\n'
    plotstuff += 'sc_kf <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, max_kf+1))\n'
    plotstuff += 'sc_rf <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, max_rf+1))\n'
    plotstuff += '\n'

    plotstuff += '# color smc by kf distances\n'
    plotstuff += 'pdf("smc_kf_distances.pdf")\n'
    plotstuff += 'p_kf_smc <- ggplot(df, aes(theta/2, T, color=kf_smc)) + sc_kf\n'
    plotstuff += 'p_kf_smc + geom_point()\n'
    plotstuff += 'dev.off()\n'
    plotstuff += '\n'
    plotstuff += '# color beast by kf distances\n'
    plotstuff += 'pdf("beast_kf_distances.pdf")\n'
    plotstuff += 'p_kf_beast <- ggplot(df, aes(theta/2, T, color=kf_beast)) + sc_kf\n'
    plotstuff += 'p_kf_beast + geom_point()\n'
    plotstuff += 'dev.off()\n'
    plotstuff += '\n'
    plotstuff += '# color smc by rf distances\n'
    plotstuff += 'pdf("smc_rf_distances.pdf")\n'
    plotstuff += 'p_rf_smc <- ggplot(df, aes(theta/2, T, color=rf_smc)) + sc_rf\n'
    plotstuff += 'p_rf_smc + geom_point()\n'
    plotstuff += 'dev.off()\n'
    plotstuff += '\n'
    plotstuff += '# color beast by rf distances\n'
    plotstuff += 'pdf("beast_rf_distances.pdf")\n'
    plotstuff += 'p_rf_beast <- ggplot(df, aes(theta/2, T, color=rf_beast)) + sc_rf\n'
    plotstuff += 'p_rf_beast + geom_point()\n'
    plotstuff += 'dev.off()\n'
    plotstuff += '\n'

    plotstuff += 'deep_coal <- read.csv("deep_coal.txt", header=FALSE)\n'
    plotstuff += 'deep_coal <- as.numeric(deep_coal)\n'
    plotstuff += '\n'
    plotstuff += '# plot number of deep coalescences vs KF distance -smc \n'
    plotstuff += 'pdf("deep_coal_kf_smc.pdf")\n'
    plotstuff += 'df_smc_kf <- data.frame(kf_smc, deep_coal)\n'
    plotstuff += 'p_smc_deep_kf <- ggplot(df_smc_kf, aes(deep_coal, kf_smc))\n'
    plotstuff += 'p_smc_deep_kf+geom_point() + theme_classic() + ylim(0, max_kf+1)\n'
    plotstuff += 'dev.off()\n'
    plotstuff += '\n'
    plotstuff += '# plot number of deep coalescences vs KF distance -beast\n'
    plotstuff += 'pdf("deep_coal_kf_beast.pdf")\n'
    plotstuff += 'df_beast_kf <- data.frame(kf_beast, deep_coal)\n'
    plotstuff += 'p_beast_deep_kf <- ggplot(df_beast_kf, aes(deep_coal, kf_beast))\n'
    plotstuff += 'p_beast_deep_kf+geom_point() + theme_classic() + ylim(0, max_kf+1)\n'
    plotstuff += 'dev.off()\n'
    plotstuff += '\n'
    plotstuff += '# plot number of deep coalescences vs RF distance - smc\n'
    plotstuff += 'pdf("deep_coal_rf_smc.pdf")\n'
    plotstuff += 'df_smc_rf <- data.frame(rf_smc, deep_coal)\n'
    plotstuff += 'p_smc_deep_rf <- ggplot(df_smc_rf, aes(deep_coal, rf_smc))\n'
    plotstuff += 'p_smc_deep_rf+geom_point() + theme_classic() + ylim(0, max_rf+1)\n'
    plotstuff += 'dev.off()\n'
    plotstuff += '\n'
    plotstuff += '# plot number of deep coalescences vs RF distance -beast\n'
    plotstuff += 'pdf("deep_coal_rf_beast.pdf")\n'
    plotstuff += 'df_beast_rf <- data.frame(rf_beast, deep_coal)\n'
    plotstuff += 'p_beast_deep_rf <- ggplot(df_beast_rf, aes(deep_coal, rf_beast))\n'
    plotstuff += 'p_beast_deep_rf+geom_point() + theme_classic() + ylim(0, max_rf+1)\n'
    plotstuff += 'dev.off()\n'
    plotstuff += '\n'
    plotstuff += '# make bar plot of differences between smc / true theta and beast / true theta\n'
    plotstuff += 'thetadf <- read.csv("theta-comparisons.txt", header=TRUE)\n'
    plotstuff += 'tibble::rowid_to_column(thetadf, "ID") # add indices\n'
    plotstuff += 'smc_dif <- thetadf$smc - thetadf$true\n'
    plotstuff += 'beast_dif <- thetadf$beast - thetadf$true\n'
    plotstuff += 'thetadf$smc_dif <- smc_dif\n'
    plotstuff += 'thetadf$beast_dif <- beast_dif\n'
    plotstuff += '\n'
    plotstuff += 'barplot(height = t(thetadf[c("smc_dif", "beast_dif")]), names.arg = thetadf$ID, beside = TRUE)\n'
    plotstuff += 'smc_mean <- mean(smc_dif)\n'
    plotstuff += 'beast_mean <- mean(beast_dif)\n'
    plotstuff += 'smc_mean\n'
    plotstuff += 'beast_mean\n'

    plotstuff == '\n'
    plotf = open(plotfn, 'w')
    plotf.write(plotstuff)
    plotf.close()

def createReplacePy():
	replacefn = os.path.join(dirname, 'rreplace.py')
	s = "import os\n"
	s += "from tempfile import NamedTemporaryFile\n"
	s += "svdq_lamdas = []\n"
	s += "lines = open('svdq_lambdas.txt', 'r').readlines()\n"
	s += "for line in lines:\n"
	s += "	svdq_lamdas.append(line)\n"
	s += "for i in range(%12d):\n" % nreps
	s += "	fname = 'rep' + str(i+1) + '/smc/proj.conf'\n"
	s += "	with open(fname) as fin, NamedTemporaryFile(dir='.', delete=False) as fout:\n"
	s += "		for line in fin:\n"
	s += "			if line.startswith('lambda ='):\n"
	s += "				line = 'lambda = ' + str(svdq_lamdas[i])\n"
	s += "			fout.write(line.encode('utf8'))\n"
	s += "		os.rename(fout.name, fname)\n"
	s += "		i += 1\n"
	s += "\n"
	s += "for i in range(%12d):\n" % nreps
	s += "	fname = 'rep' + str(i+1) + '/beast/starbeast.xml'\n"
	s += "	with open(fname) as fin, NamedTemporaryFile(dir='.', delete=False) as fout:\n"
	s += "		for line in fin:\n"
	s += "			if 'SVDQLAMBDA' in line: #if line contains placeholder\n"
	s += "				to_replace = str(svdq_lamdas[i])\n"
	s += "				to_replace = to_replace.strip()\n"
	s += '				line = "			'
	s += "<parameter id='speciationRate.t:Species' spec='parameter.RealParameter' estimate='false' lower='%s' name='birthRate' upper='%s'>%s</parameter>\\n	"
	s += '"'
	s += ' % (to_replace, to_replace, to_replace)\n'
#	s += "				line = '<parameter id='speciationRate.t:Species' spec='parameter.RealParameter' estimate='false' lower='%s' name='birthRate' upper='%s'>'%s'</parameter>\\n' % (to_replace, to_replace, to_replace)\n"
	s += "			fout.write(line.encode('utf8'))\n"
	s += "		os.rename(fout.name, fname)\n"
	s += "		i += 1\n"
	replacef = open(replacefn, 'w')
	replacef.write(s)
	replacef.close()


def createSimBash():
    simfn = os.path.join(dirname, 'simulate.sh')
    simbash = '#!/bin/bash\n'

    for rep in range(nreps):
        repdirpath = createRepDirName(rep)
        simbash += '\n'
        simbash += 'cd %s\n' % os.path.join(repdirpath, 'sim')
        simbash += '%s\n' % simprogname
        simbash += '%s svd-qage.nex\n' % paupprogname
        simbash += 'cd ~-\n'

    simf = open(simfn, 'w')
    simf.write(simbash)
    simf.close()

def createSMCSlurm():
    # see https://blog.ronin.cloud/slurm-job-arrays/
    smcslurmfn = os.path.join(dirname, 'smcslurm.sh')
    s  = ''
    s += '#!/bin/bash\n'
    s += '\n'
    partition  = nodechoices[nodechoice][0]
    constraint = nodechoices[nodechoice][1]
    if partition == 'general':
        s += '#SBATCH -p general\n'
    else:
        s += '#SBATCH -p priority\n'
        s += '#SBATCH -q pol02003sky\n'
    s += '#SBATCH -C \'%s\'\n' % constraint
    s += '#SBATCH -A pol02003\n'
    s += '#SBATCH --nodes=1\n'
    s += '#SBATCH --ntasks=34\n'
    if maxsimult is None:
        s += '#SBATCH --array=1-%d\n' % (nreps,)
    else:
        s += '#SBATCH --array=1-%d%%%d\n' % (nreps, maxsimult)
    s += '#SBATCH --job-name=smc\n'
    s += '#SBATCH -o smc-%a.out\n'
    s += '#SBATCH -e smc-%a.err\n'
    s += '#SBATCH --mem=50G\n'
    s += '\n'
    s += 'LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/lib"\n'
    s += 'export TIMEFORMAT="user-seconds %3U"\n'
    s += 'cd /home/%s/%s/rep${SLURM_ARRAY_TASK_ID}/smc\n' % (username, dirname)
    s += 'time $HOME/bin/%s\n' % smcprogname

    smcslurmf = open(smcslurmfn, 'w')
    smcslurmf.write(s)
    smcslurmf.close()

def createBeastSlurm():
    # see https://blog.ronin.cloud/slurm-job-arrays/
    beastslurmfn = os.path.join(dirname, 'beastslurm.sh')
    s  = ''
    s += '#!/bin/bash\n'
    s += '\n'
    partition = nodechoices[nodechoice][0]
    constraint = nodechoices[nodechoice][1]
    if partition == 'general':
        s += '#SBATCH -p general\n'
    else:
        s += '#SBATCH -p priority\n'
        s += '#SBATCH -q pol02003sky\n'
    s += '#SBATCH -C \'%s\'\n' % constraint
    s += '#SBATCH -A pol02003\n'
    s += '#SBATCH --nodes=1\n'
    if maxsimult is None:
        s += '#SBATCH --array=1-%d\n' % (nreps,)
    else:
        s += '#SBATCH --array=1-%d%%%d\n' % (nreps, maxsimult)
    s += '#SBATCH --job-name=beast\n'
    s += '#SBATCH -o beast-%a.out\n'
    s += '#SBATCH -e beast-%a.err\n'
    s += '#SBATCH --mem=150G\n'
    s += '\n'
    s += 'LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/lib"\n'
    s += 'export TIMEFORMAT="user-seconds %3U"\n'
    s += 'cd /home/%s/%s/rep${SLURM_ARRAY_TASK_ID}/beast\n' % (username, dirname)
    s += 'time $HOME/beast/bin/%s starbeast.xml\n' % beastprogname

    beastslurmf = open(beastslurmfn, 'w')
    beastslurmf.write(s)
    beastslurmf.close()

def createCopyDataPy():
    copydatafn =  os.path.join(dirname, 'copydata.py')
    shutil.copyfile('copydata_template.py', copydatafn)
    stuff = open(copydatafn, 'r').read()
    stuff, n = re.subn('__NLOCI__', '%d' % nloci, stuff, re.M | re.S)
    assert n == 1
    stuff, n = re.subn('__SEQLEN__', '%d' % seqlen, stuff, re.M | re.S)
    assert n == 1
    copydataf = open(copydatafn, 'w')
    copydataf.write(stuff)
    copydataf.close()

def createANOVAPy():
    anovafn =  os.path.join(dirname, 'anova-means.py')
    shutil.copyfile('anova-means-template.py', anovafn)
    stuff = open(anovafn, 'r').read()
    stuff, n = re.subn('__NREPS__', '%d' % nreps, stuff, re.M | re.S)
    assert n == 1
    stuff, n = re.subn('__SAMPLESIZE__', '%d' % nparticles, stuff, re.M | re.S)
    assert n == 1
    anovaf = open(anovafn, 'w')
    anovaf.write(stuff)
    anovaf.close()
    
    anovafn =  os.path.join(dirname, 'anova-full.py')
    shutil.copyfile('anova-full-template.py', anovafn)
    stuff = open(anovafn, 'r').read()
    stuff, n = re.subn('__NREPS__', '%d' % nreps, stuff, re.M | re.S)
    assert n == 1
    stuff, n = re.subn('__SAMPLESIZE__', '%d' % nparticles, stuff, re.M | re.S)
    assert n == 1
    anovaf = open(anovafn, 'w')
    anovaf.write(stuff)
    anovaf.close()

def createCrunch():
    # see https://blog.ronin.cloud/slurm-job-arrays/
    # see https://statsandr.com/blog/how-to-one-way-anova-by-hand/
    crunchfn = os.path.join(dirname, 'crunch.py')

    s   = ''
    s  += 'from math import sqrt\n'
    s  += '\n'
    s  += 'class DistSummary:\n'
    s  += '    def __init__(self):\n'
    s  += '        self.dists = {}\n'
    s  += '        self.count = {}\n'
    s  += '        self.total = {}\n'
    s  += '        self.sum = {}\n'
    s  += '        self.sumsq = {}\n'
    s  += '        self.cum = {}\n'
    s  += '        self.mean = {}\n'
    s  += '        self.var = {}\n'
    s  += '        self.stdev = {}\n'
    s  += '        self.group_mean = {}\n'
    s  += '        self.group_n = {}\n'
    s  += '    def zero(self, rep):\n'
    s  += '        self.dists[rep] = []\n'
    s  += '        self.count[rep] = 0\n'
    s  += '        self.total[rep] = 0\n'
    s  += '        self.sum[rep] = 0.0\n'
    s  += '        self.sumsq[rep] = 0.0\n'
    s  += '        self.cum[rep] = 0.0\n'
    s  += '        self.mean[rep] = 0.0\n'
    s  += '        self.var[rep] = 0.0\n'
    s  += '        self.stdev[rep] = 0.0\n'
    s  += '        self.group_mean[rep] = 0.0\n'
    s  += '        self.group_n[rep] = 0\n'
    s  += '\n'
    s  += 'def getKFDistances(fnprefix):\n'
    s  += '    d = DistSummary()\n'
    s  += '    for rep in range(%d):\n' % nreps
    s  += '        d.zero(rep)\n'
    s  += '        lines = open("%s%d.txt" % (fnprefix, rep+1,), "r").readlines()\n'
    s  += '        for line in lines[1:]:\n'
    s  += '            parts = line.strip().split()\n'
    s  += '            assert len(parts) == 3\n'
    s  += '            y = float(parts[1])\n'
    s  += '            d.dists[rep].append(y)\n'
    s  += '\n'
    s  += '        d.count[rep] = len(d.dists[rep])\n'
    s  += '        d.sum[rep] = sum(d.dists[rep])\n'
    s  += '        d.sumsq[rep] = sum([y*y for y in d.dists[rep]])\n'
    s  += '        d.cum[rep] += d.sum[rep]\n'
    s  += '        d.total[rep] += d.count[rep]\n'
    s  += '        d.mean[rep] = d.sum[rep]/d.count[rep]\n'
    s  += '        d.var[rep] = (d.sumsq[rep] - pow(d.mean[rep],2.)*d.count[rep])/(d.count[rep]-1)\n'
    s  += '        if d.var[rep] < 0.0:\n'
    s  += '            print("warning: variance negative (%g) for %s rep %d: mean = %g, sumsq = %g, count = %d" % (d.var[rep], fnprefix, rep, d.mean[rep], d.sumsq[rep], d.count[rep]))\n'
    s  += '            d.var[rep] = 0.0\n'
    s  += '        d.stdev[rep] = sqrt(d.var[rep])\n'
    s  += '        d.group_mean[rep] = d.mean[rep]\n'
    s  += '        d.group_n[rep] = d.count[rep]\n'
    s  += '    return d\n'
    s  += '\n'
    s  += '\n'
    s  += 'def getRFDistances(fnprefix):\n'
    s  += '    d = DistSummary()\n'
    s  += '    for rep in range(%d):\n' % nreps
    s  += '        d.zero(rep)\n'
    s  += '        lines = open("%s%d.txt" % (fnprefix, rep+1,), "r").readlines()\n'
    s  += '        for line in lines[1:]:\n'
    s  += '            parts = line.strip().split()\n'
    s  += '            assert len(parts) == 3\n'
    s  += '            y = float(parts[2])\n'
    s  += '            d.dists[rep].append(y)\n'
    s  += '\n'
    s  += '        d.count[rep] = len(d.dists[rep])\n'
    s  += '        d.sum[rep] = sum(d.dists[rep])\n'
    s  += '        d.sumsq[rep] = sum([y*y for y in d.dists[rep]])\n'
    s  += '        d.cum[rep] += d.sum[rep]\n'
    s  += '        d.total[rep] += d.count[rep]\n'
    s  += '        d.mean[rep] = d.sum[rep]/d.count[rep]\n'
    s  += '        d.var[rep] = (d.sumsq[rep] - pow(d.mean[rep],2.)*d.count[rep])/(d.count[rep]-1)\n'
    s  += '        if d.var[rep] < 0.0:\n'
    s  += '            print("warning: variance negative (%g) for %s rep %d: mean = %g, sumsq = %g, count = %d" % (d.var[rep], fnprefix, rep, d.mean[rep], d.sumsq[rep], d.count[rep]))\n'
    s  += '            d.var[rep] = 0.0\n'
    s  += '        d.stdev[rep] = sqrt(d.var[rep])\n'
    s  += '        d.group_mean[rep] = d.mean[rep]\n'
    s  += '        d.group_n[rep] = d.count[rep]\n'
    s  += '    return d\n'
    s  += 'kf = open("kf-summary.txt", "x")\n'
    s  += 'rf = open("rf-summary.txt", "x")\n'
    s  += 'dsmc_kf = getKFDistances("smcdists")\n'
    s  += 'dbeast_kf = getKFDistances("beastdists")\n'
    s  += 'dsmc_rf = getRFDistances("smcdists")\n'
    s  += 'dbeast_rf = getRFDistances("beastdists")\n'
    s  += 'print("%12s %38s %38s" % ("replicate", "----------------- SMC ----------------", "---------------- BEAST ---------------"))\n'
    s  += 'print("%12s %12s %12s %12s %12s %12s %12s" % ("replicate", "count", "mean", "stdev", "count", "mean", "stdev"))\n'
    s  += 'for rep in range(%d):\n' % nreps
    s  += '    print("kf: %12d %12d %12.5f %12.5f %12d %12.5f %12.5f" % (rep+1, dsmc_kf.count[rep], dsmc_kf.mean[rep], dsmc_kf.stdev[rep], dbeast_kf.count[rep], dbeast_kf.mean[rep], dbeast_kf.stdev[rep]))\n'
    s +=  '    kf.write("%12d %12d %12.5f %12.5f %12d %12.5f %12.5f \\n" % (rep+1, dsmc_kf.count[rep], dsmc_kf.mean[rep], dsmc_kf.stdev[rep], dbeast_kf.count[rep], dbeast_kf.mean[rep], dbeast_kf.stdev[rep]))\n'
    s +=  '    kf.close\n'

    s  += 'for rep in range(%d):\n' % nreps
    s  += '    print("rf: %12d %12d %12.5f %12.5f %12d %12.5f %12.5f" % (rep+1, dsmc_rf.count[rep], dsmc_rf.mean[rep], dsmc_rf.stdev[rep], dbeast_rf.count[rep], dbeast_rf.mean[rep], dbeast_rf.stdev[rep]))\n'
    s  += '    rf.write("%12d %12d %12.5f %12.5f %12d %12.5f %12.5f \\n" % (rep+1, dsmc_rf.count[rep], dsmc_rf.mean[rep], dsmc_rf.stdev[rep], dbeast_rf.count[rep], dbeast_rf.mean[rep], dbeast_rf.stdev[rep]))\n'
    s  += '    rf.close\n'

#    s  += 'dsvdq_kf = getKFDistances("svdqdists")\n'
#    s  += 'dsvdq_rf = getRFDistances("svdqdists")\n'
#    s  += 'print("%12s %38s %38s %38s" % ("replicate", "----------------- SMC ----------------", "---------------- BEAST ---------------", ---------------- SVDQ ---------------"))\n'
#    s  += 'print("%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s" % ("replicate", "count", "mean", "stdev", "count", "mean", "stdev", "count", "mean", "stdev"))\n'
#    s  += 'for rep in range(%d):\n' % nreps
#    s  += '    print("kf: %12d %12d %12.5f %12.5f %12d %12.5f %12.5f" % (rep+1, dsmc_kf.count[rep], dsmc_kf.mean[rep], dsmc_kf.stdev[rep], dbeast_kf.count[rep], dbeast_kf.mean[rep], dbeast_kf.stdev[rep]))\n'
#    s +=  '    kf.write("%12d %12d %12.5f %12.5f %12d %12.5f %12.5f \\n" % (rep+1, dsmc_kf.count[rep], dsmc_kf.mean[rep], dsmc_kf.stdev[rep], dbeast_kf.count[rep], dbeast_kf.mean[rep], dbeast_kf.stdev[rep]))\n'
#    s +=  '    kf.close\n'
#    s  += 'for rep in range(%d):\n' % nreps
#    s  += '    print("rf: %12d %12d %12.5f %12.5f %12d %12.5f %12.5f %12d %12.5f %12.5f" % (rep+1, dsmc_rf.count[rep], dsmc_rf.mean[rep], dsmc_rf.stdev[rep], dbeast_rf.count[rep], dbeast_rf.mean[rep], dbeast_rf.stdev[rep], dsvdq_rf.count[rep], dsvdq_rf.mean[rep], dsvdq_rf.stdev[rep]))\n'
#    s  += '    rf.write("%12d %12d %12.5f %12.5f %12d %12.5f %12.5f %12d %12.5f %12.5f \\n" % (rep+1, dsmc_rf.count[rep], dsmc_rf.mean[rep], dsmc_rf.stdev[rep], dbeast_rf.count[rep], dbeast_rf.mean[rep], dbeast_rf.stdev[rep], , dsvdq_rf.count[rep], dsvdq_rf.mean[rep], dsvdq_rf.stdev[rep]))\n'
#    s  += '    rf.close\n'
    s  += 'print(" ")\n'

    crunchf = open(crunchfn, 'w')
    crunchf.write(s)
    crunchf.close()

def createPAUP(pathname, fn, startat):
    # see https://blog.ronin.cloud/slurm-job-arrays/
    paupfn = os.path.join(dirname, '%spaup.nex' % pathname)

    s   = ''
    s  += '#nexus\n'
    s  += '\n'
    s  += 'begin paup;\n'
    s  += '    set maxtrees=%d;\n' % (nparticles+1)
    for rep in range(nreps):
        s  += '\n'
        s  += '    [### rep%d ###]\n' % (rep+1,)
        s  += '    gettrees file=rep%d/sim/true-species-tree.tre;\n' % (rep+1,)
        s  += '    gettrees file=rep%d/%s/%s mode=7 from=%d;\n' % (rep+1, pathname, fn, startat)
        s  += '    treedist reftree=1 measure=KF file=%sdists%d.txt replace;\n' % (pathname, rep+1)
        s  += '    cleartrees;\n'
    s  += '    quit;\n'
    s  += 'end;\n'

    paupf = open(paupfn, 'w')
    paupf.write(s)
    paupf.close()

def createTreeDist(pathname, fn, startat):
    # see https://blog.ronin.cloud/slurm-job-arrays/
    if (pathname == 'svdq') : #svdq tree files are in a different location
        tdfn = os.path.join(dirname, '%std.sh' % pathname)

        s   = '#!/bin/bash\n'
        for rep in range(nreps):
            s  += '\n\n### rep%d ###\n' % (rep+1,)
            s  += 'td --reffile rep%d/sim/true-species-tree.tre --treefile rep%d/sim/%s --skip %d --reftree 1 --outfile %sdists%d.txt\n' % (rep+1,rep+1, fn, startat, pathname, rep+1)

    
    else :
        tdfn = os.path.join(dirname, '%std.sh' % pathname)

        s   = '#!/bin/bash\n'
        for rep in range(nreps):
            s  += '\n\n### rep%d ###\n' % (rep+1,)
            s  += 'td --reffile rep%d/sim/true-species-tree.tre --treefile rep%d/%s/%s --skip %d --reftree 1 --outfile %sdists%d.txt\n' % (rep+1,rep+1, pathname, fn, startat, pathname, rep+1)

    tdf = open(tdfn, 'w')
    tdf.write(s)
    tdf.close()

def writeDeepCoalFile():
	deepcoalfn = os.path.join(dirname, 'get-deep-coal.py')
	s = "i = 1\n"
	s += "new_f = open('deep_coal.txt', 'x')\n"
	s += "for rep in range(%12d):\n" % nreps
	s += "	name = 'rep' + str(i)\n"
	s += "	lines = open(name+'/sim'+'/deep_coalescences.txt', 'r').readlines()[1:]\n"
	s += "	for line in lines[1:]:\n"
	s += "		parts = line.strip().split()\n"
	s += "		if (i == 1) :\n"
	s += "			new_f.write(parts[1])\n"
	s += "		else:\n"
	s += "			new_f.write(',' + parts[1])\n"
	s += "	i+=1\n"
	s += "\n"
	s += "new_f.close\n"
	deepcoalf = open(deepcoalfn, 'w')
	deepcoalf.write(s)
	deepcoalf.close()

def writeThetaFile():
	thetafn = os.path.join(dirname, 'get-average-theta.py')
	s = "import os\n"
	s += "import numpy\n"
	s += "import shutil\n"
	s += "import pandas as pd\n"
	s += "import csv\n"
	s +=  "nreps = %12d \n" % nreps
	s += "smc_theta_list = []\n"
	s += "beast_theta_list = []\n"
	s += "true_theta_list = []\n"
	s += "for i in range(%12d):" % nreps
	s += "	# get smc average theta\n"
	s += "	smc = pd.read_csv('rep' + str(i+1) + '/smc/params-beast-comparison.log', sep='\t')\n"
	s += "	theta_smc = smc['popMean '].mean()*4\n"
	s += "	smc_theta_list.append(theta_smc)\n"
	s += "	# get beast average theta\n"
	s += "	beast = pd.read_csv('rep' + str(i+1) + '/beast/starbeast3.log', sep = '\t', comment = '#')\n"
	s += "	theta_beast = beast['popMean'].mean()*4\n"
	s += "	beast_theta_list.append(theta_beast)\n"
	s += "	# get true theta\n"
	s += "	with  open('rep' + str(i+1) + '/sim/proj.conf') as fi:\n"
	s += "		id = []\n"
	s += "		for ln in fi:\n"
	s += "			if ln.startswith('theta'):\n"
	s += "				true_theta_list.append(float(ln[9:13]))\n"
	s += "rows = zip(smc_theta_list, beast_theta_list, true_theta_list)\n"
	s += "with open('theta-comparisons.txt', 'w') as f:\n"
	s += "	f.write('smc, beast, true\\n')\n"
	s += "	writer = csv.writer(f)\n"
	s += "	for row in rows:\n"
	s += "		writer.writerow(row)\n"
	thetaf = open(thetafn, 'w')
	thetaf.write(s)
	thetaf.close()

def writeTimeFile():
	timefn = os.path.join(dirname, 'get-times.py')
	s = "import os\n"
	s += "import numpy\n"
	s += "import shutil\n"
	s +=  "nreps = %12d \n" % nreps
	s += "timef = open ('times.txt', 'x')\n"
	s += "smc_time_list = []\n"
	s += "beast_time_list = []\n"
	s += "for i in range(nreps):\n"
	s += "	smc_file_name = 'smc-' + str(i+1) + '.err'\n"
	s += "	with open (smc_file_name, mode = 'r') as file:\n"
	s += "		for line in file:\n"
	s += "			pass\n"
	s += "		last_line = line\n"
	s += "	smc_time_list.append(float(last_line[13:]))\n"
	s +=  "	beast_file_name = 'beast-' + str(i+1) + '.err'\n"
	s += "	with open (beast_file_name, mode='r') as file:\n"
	s += "		for line in file:\n"
	s += "			pass\n"
	s += "		last_line = line\n"
	s += "	beast_time_list.append(float(last_line[13:]))\n"
	s += "smc_average = sum(smc_time_list) / len(smc_time_list)\n"
	s += "timef.write('smc average time: ' + str(smc_average))\n"
	s += "print('smc average time: ' + str(smc_average))\n"
	s += "beast_average = sum(beast_time_list) / len(beast_time_list)\n"
	s += "timef.write('beast average time: ' + str(beast_average))\n"
	s += "print('beast average time: ' + str(beast_average))\n"

	timef = open(timefn, 'w')
	timef.write(s)
	timef.close()

def creatergl3DPLOT():
	plotfn = os.path.join(dirname, 'rgl.R')
	s = "library(rgl)\n"
	if method == 'lognorm':
		Tstr = ['%g' % t for t in Tvect]
		s += 'T = c(%s)\n' % ','.join(Tstr)

		Rstr = ['%g' % r for r in Rvect]
		s += 'R = c(%s)\n' % ','.join(Rstr)

		thetastr = ['%g' % q for q in thetas]
		s += 'theta = c(%s)\n' % ','.join(thetastr)

		lambdastr = ['%g' % l for l in lambdas]
		s += 'lambda = c(%s)\n' % ','.join(lambdastr)

		s += 'plot(theta, lambda, type="p", pch=19, main="Simulation conditions", xlab="theta", ylab="lambda")\n'
	elif method == 'uniform':
		Tstr = ['%g' % t for t in Tvect]
		s += 'T = c(%s)\n' % ','.join(Tstr)

		thetastr = ['%g' % q for q in thetas]
		s += 'theta = c(%s)\n' % ','.join(thetastr)

		s += 'plot(theta/2, T, type="p", pch=19, main="Simulation conditions", xlab="theta/2", ylab="T")\n'
	elif method == 'grid':
		Tstr = ['%g' % t for t in Tvect]
		s += 'T = c(%s)\n' % ','.join(Tstr)

		thetastr = ['%g' % q for q in thetas]
		s += 'theta = c(%s)\n' % ','.join(thetastr)

		s += 'plot(theta/2, T, type="p", pch=19, main="Simulation conditions", xlab="theta/2", ylab="T")\n'
	else:
		assert False, 'method should be either "lognorm" or "uniform" but you specified "%s"' % method
	s += '\n'
	s += 'rf <- read.table(file="rf-summary.txt", header=FALSE)\n'
	s += 'df <- c("particle", "col1", "rf_smc", "col3", "col4", "rf_beast", "col6")\n'
	s += 'colnames(rf) <- df\n'
	s += 'RF_smc <- rf$V3\n'
	s += 'RF_beast <- rf$V6\n'
	s += '\n'
	s += 'kf <- read.table(file="kf-summary.txt")\n'
	s += 'KF_smc <- kf$V3\n'
	s += 'KF_beast <- kf$V6\n'
	s += '\n'
	s += 'theta_over_two <- theta/2\n'
	s += '\n'
	s += '#plot smc\n'
	s += 'open3d()\n'
	s += 'plot3d(x=theta_over_two, y=T, z=RF_smc, type="h", col="navy", size=1, lwd=1, xlab="theta/2", ylab="T", zlab="RF", box=FALSE, zlim = c(0,8))\n'
	s += 'spheres3d(x=theta_over_two, y=T, z=RF_smc, radius=.05, col="navy")\n'
	s += 'grd <- expand.grid(x=c(0, max(theta_over_two)), y=c(0,max(T)), z=0)\n'
	s += '\n'
	s += '# prediction surface\n'
	s += 'material3d(color = "red")\n'
	s += 'persp3d(x=unique(grd[[1]]), y=unique(grd[[2]]), z=matrix(grd[[3]],2,2), add=TRUE)\n'
	s += '\n'
	s += '# plot beast\n'
	s += 'open3d()\n'
	s += 'plot3d(x=theta_over_two, y=T, z=RF_beast, type="h", col="navy", size=1, lwd=1, xlab="theta/2", ylab="T", zlab="RF", box=FALSE, zlim = c(0,8))\n'
	s += 'spheres3d(x=theta_over_two, y=T, z=RF_beast, radius=.05, col="navy")\n'
	s += 'grd <- expand.grid(x=c(0, max(theta_over_two)), y=c(0,max(T)), z=0)\n'
	s += '\n'
	s += '# prediction surface\n'
	s += 'material3d(color = "red")\n'
	s += 'persp3d(x=unique(grd[[1]]), y=unique(grd[[2]]), z=matrix(grd[[3]],2,2), add=TRUE)\n'
	s += '\n'
	s += '# Replot SMC\n'
	s += 'open3d(zoom = zoom, userMatrix = userMatrix, windowRect=windowRect)\n'
	s += 'bg3d(color=bgcolor)\n'
	s += 'light3d(theta=45)\n'
	s += 'persp3d(x=T_vals, y=theta_vals, z=matrix(rf$rf_smc - rf$rf_beast,nrows,ncols), col=surfcolor, zlim=c(-8,8), xlab="T", ylab="theta", zlab="delta RF")\n'
	s += 'grd <- expand.grid(x=c(0, max(T_vals)), y=c(0,max(theta_vals)), z=0)\n'
	s += 'material3d(color = planecolor)\n'
	s += 'persp3d(x=unique(grd[[1]]), y=unique(grd[[2]]), z=matrix(grd[[3]],2,2), col=planecolor,add=TRUE)\n'
	s += '\n'
	s += '# Save to file\n'
	s += 'rgl.snapshot("rfdiff-smc-on-top.png", fmt="png", top=TRUE)\n'
	s += '\n'
	s += '# Plot BEAST\n'
	s += 'open3d(zoom = zoom, userMatrix = userMatrix, windowRect=windowRect)\n'
	s += 'bg3d(color=bgcolor)\n'
	s += 'light3d(theta=45)\n'
	s += 'persp3d(x=T_vals, y=theta_vals, z=matrix(rf$rf_beast - rf$rf_smc,nrows,ncols), col=surfcolor, zlim=c(-8,8), xlab="T", ylab="theta", zlab="delta RF")\n'
	s += 'grd <- expand.grid(x=c(0, max(T_vals)), y=c(0,max(theta_vals)), z=0)\n'
	s += 'material3d(color = planecolor)\n'
	s += 'persp3d(x=unique(grd[[1]]), y=unique(grd[[2]]), z=matrix(grd[[3]],2,2), col=planecolor,add=TRUE)\n'
	s += '\n'
	s += '# Save to file\n'
	s += 'rgl.snapshot("rfdiff-beast-on-top.png", fmt="png", top=TRUE)\n'
	plotf = open(plotfn, 'w')
	plotf.write(s)
	plotf.close()


def create3DRPlot():
	plotfn = os.path.join(dirname, 'plotdistances.py')
	s = "import plotly.express as px\n"
	s += "import plotly.figure_factory as ff\n"
	s += "import plotly.graph_objs as go\n"
	s += "import pandas as pd\n"
	s += "import math, random, sys, os\n"
	s += "import scipy,numpy\n"
	s += "import numpy as np\n"
	s += "from scipy.interpolate import griddata\n"
	s += "\n"
	s += "file = open('kf-summary.txt','r')\n"
	s += "lines = file.readlines()\n"
	s += "file.close()\n"
	s += "\n"
	s += "kf_smc = []\n"
	s += "kf_beast = []\n"
	s += "\n"
	s += "for line in lines:\n"
	s += "  parts = line.split() # split line into parts\n"
	s += "  column3 = parts[2]\n"
	s += "  kf_smc.append(column3)\n"
	s += "  column6 = parts[5]\n"
	s += "  kf_beast.append(column6)\n"
	s += "\n"
	s += "file = open('rf-summary.txt','r')\n"
	s += "lines = file.readlines()\n"
	s += "file.close()\n"
	s += "\n"
	s += "rf_smc = []\n"
	s += "rf_beast = []\n"
	s += "\n"
	s += "for line in lines:\n"
	s += "  parts = line.split() # split line into parts\n"
	s += "  column3 = parts[2]\n"
	s += "  rf_smc.append(column3)\n"
	s += "  column6 = parts[5]\n"
	s += "  rf_beast.append(column6)\n"
	s += "\n"
	s += "plot_type = 'scatterplot'\n"
	s += "#plot_type = 'meshplot'\n"
	s += "plot_type = 'surface'\n"
	s += "\n"
	s += "plot_saved_to_file = False\n"
	s += "plotting_filename = '%s.png' % plot_type\n"
	s += "plotting_height = 1000\n"
	s += "plotting_scale = 1\n"
	s += "\n"
	if method == 'lognorm':
		Tstr = ['%g' % t for t in Tvect]
		s += 'T = [%s]\n' % ','.join(Tstr)
		Rstr = ['%g' % r for r in Rvect]
		s += 'R = [%s]\n' % ','.join(Rstr)
		thetastr = ['%g' % q for q in thetas]
		s += 'theta = [%s]\n' % ','.join(thetastr)
		lambdastr = ['%g' % l for l in lambdas]
		s += 'lambda = c(%s)\n' % ','.join(lambdastr)
	elif method == 'uniform':
		Tstr = ['%g' % t for t in Tvect]
		s += 'T = [%s]\n' % ','.join(Tstr)
		thetastr = ['%g' % q for q in thetas]
		s += 'theta = [%s]\n' % ','.join(thetastr)
	elif method == 'grid':
		Tstr = ['%g' % t for t in Tvect]
		s += 'T = [%s]\n' % ','.join(Tstr)
		thetastr = ['%g' % q for q in thetas]
		s += 'theta = [%s]\n' % ','.join(thetastr)
	else:
		assert False, 'method should be either "lognorm" or "uniform" or "grid" but you specified "%s"' % method
	s += 'RF_smc = [float(i) for i in rf_smc]\n'
	s += 'RF_beast = [float(i) for i in rf_beast]\n'
	s += 'KF_smc = [float(i) for i in kf_smc]\n'
	s += 'KF_beast = [float(i) for i in kf_beast]\n'
	s += 'halfT     = [0.5*t for t in T]\n'
	s += 'KFdiff    = [kfs - kfb for (kfs,kfb) in zip(KF_smc, KF_beast)]\n'
	s += 'KFabsdiff = [math.fabs(kfs - kfb) for (kfs,kfb) in zip(KF_smc, KF_beast)]\n'
	s += 'RFdiff    = [rfs - rfb for (rfs,rfb) in zip(RF_smc, RF_beast)]\n'
	s += '#RFdiff    = [rfs - rfb for (rfs,rfb) in zip(RF_beast, RF_smc)]\n'
	s += 'RFabsdiff = [math.fabs(rfs - rfb) for (rfs,rfb) in zip(RF_smc, RF_beast)]\n'
	s += '\n'
	s += 'if plot_type == "scatterplot":\n'
	s += '	fig = go.Figure(data=[go.Scatter3d(\n'
	s += '		x=halfT,\n'
	s += '		y=theta,\n'
	s += '		z=RFdiff,\n'
	s += '		mode="markers",\n'
	s += '		marker=dict(\n'
	s += '			size=12,\n'
	s += '			color=RFdiff,                # set color to an array/list of desired values\n'
	s += '			colorscale="Rainbow",   # choose a colorscale\n'
	s += '			opacity=0.8\n'
	s += '		)\n'
	s += '	)])\n'
	s += 'elif plot_type == "meshplot":\n'
	s += '	fig = go.Figure(data=[go.Mesh3d(\n'
	s += '		x=halfT,\n'
	s += '		y=theta,\n'
	s += '		z=RFdiff,\n'
	s += '		flatshading=True,\n'
	s += '		intensity=RFabsdiff,\n'
	s += '		color="firebrick"\n'
	s += '	)])\n'
	s += 'elif plot_type == "surface":\n'
	s += '	x = np.linspace(min(halfT), max(halfT), num=100)\n'
	s += '	y = np.linspace(min(theta), max(theta), num=100)\n'
	s += '	x_grid, y_grid = np.meshgrid(x,y)\n'
	s += '	z_grid = griddata((halfT,theta),RFdiff,(x_grid,y_grid),method="cubic")\n'
	s += '\n'
	s += '	fig = go.Figure(data=[go.Surface(\n'
	s += '		x=x_grid,\n'
	s += '		y=y_grid,\n'
	s += '		z = z_grid\n'
	s += '	)])\n'
	s += '\n'
	s += '# tight layout\n'
	s += 'fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))\n'
	s += 'if plot_saved_to_file:\n'
	s += '	fig.write_image(plotting_filename, height=plotting_height, scale=plotting_scale) # this creates a static file    \n'
	s += 'else:\n'
	s += '	fig.show(renderer="browser") # this opens your browser to show you the plot now\n'
	s += '# Named colors\n'
	s += '# aliceblue, antiquewhite, aqua, aquamarine, azure, beige, bisque, black, blanchedalmond, blue, blueviolet, brown, burlywood, cadetblue, chartreuse, chocolate, coral, cornflowerblue, cornsilk, crimson, cyan, darkblue, darkcyan, darkgoldenrod, darkgray, darkgrey, darkgreen, darkkhaki, darkmagenta, darkolivegreen, darkorange, darkorchid, darkred, darksalmon, darkseagreen, darkslateblue, darkslategray, darkslategrey, darkturquoise, darkviolet, deeppink, deepskyblue, dimgray, dimgrey, dodgerblue, firebrick, floralwhite, forestgreen, fuchsia, gainsboro, ghostwhite, gold, goldenrod, gray, grey, green, greenyellow, honeydew, hotpink, indianred, indigo, ivory, khaki, lavender, lavenderblush, lawngreen, lemonchiffon, lightblue, lightcoral, lightcyan, lightgoldenrodyellow, lightgray, lightgrey, lightgreen, lightpink, lightsalmon, lightseagreen, lightskyblue, lightslategray, lightslategrey, lightsteelblue, lightyellow, lime, limegreen, linen, magenta, maroon, mediumaquamarine, mediumblue, mediumorchid, mediumpurple, mediumseagreen, mediumslateblue, mediumspringgreen, mediumturquoise, mediumvioletred, midnightblue, mintcream, mistyrose, moccasin, navajowhite, navy, oldlace, olive, olivedrab, orange, orangered, orchid, palegoldenrod, palegreen, paleturquoise, palevioletred, papayawhip, peachpuff, peru, pink, plum, powderblue, purple, red, rosybrown, royalblue, rebeccapurple, saddlebrown, salmon, sandybrown, seagreen, seashell, sienna, silver, skyblue, slateblue, slategray, slategrey, snow, springgreen, steelblue, tan, teal, thistle, tomato, turquoise, violet, wheat, white, whitesmoke, yellow, yellowgreen\n'
	plotf = open(plotfn, 'w')
	plotf.write(s)
	plotf.close()

if __name__ == '__main__':
    createMainDir()
    for rep in range(nreps):
        createRepDir(rep)
        
        createSimDir(rep)
        createSMCDir(rep)
        createBeastDir(rep)
        
        createSimConf(rep)
        createSMCConf(rep)
        createBeastXML(rep)
        
    createREADME()
    createRplot()
    createSimBash()
    createSMCSlurm()
    createBeastSlurm()
    createCopyDataPy()
    createCrunch()
    writeDeepCoalFile()
    createGetSVDQHeight()
    createReplacePy()
    #createPAUP('smc', smctreefname, 1)
    #createPAUP('beast', beasttreefname, 2)
    createTreeDist('smc', smctreefname, 1)
    createTreeDist('beast', beasttreefname, 2)
    createTreeDist('svdq', svdqtreefname, 0)
    createANOVAPy()
    writeTimeFile()
    writeThetaFile()
    create3DRPlot()
    creatergl3DPLOT()
