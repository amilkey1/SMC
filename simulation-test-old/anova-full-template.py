from math import sqrt

row_col_table = False
verbose    = False
anova_table = True

usescipy = True
try:
    from scipy.stats import f
except ImportError:
    usescipy = False
    print('Could not import scipy.stats.f')
    print('You will need to use R to compute P value')
else:
    print('Found scipy.stats.f and will use it to compute P value')

# Dimensions
a = 2
n = __NREPS__
p = __SAMPLESIZE__

def getDistances(fnprefix):
    d = {}
    for rep in range(n):
        d[rep] = []
        lines = open("%s%d.txt" % (fnprefix, rep+1,), "r").readlines()
        for line in lines[1:]:
            parts = line.strip().split()
            assert len(parts) == 2
            y = float(parts[1])
            d[rep].append(y)
    return d

dsmc = getDistances("smcdists")
dbeast = getDistances("beastdists")

nrows_smc = len(dsmc)
ncols_smc = len(dsmc[0])
nrows_beast = len(dbeast)
ncols_beast = len(dbeast[0])
if row_col_table:
    print('%12s %12s %12s' % ('Method', 'Rows', 'Columns'))
    print('%12s %12d %12d' % ('SMC', nrows_smc, ncols_smc))
    print('%12s %12d %12d' % ('BEAST', nrows_beast, ncols_beast))
assert nrows_smc == nrows_beast
assert ncols_smc == ncols_beast
nrows = nrows_smc
ncols = ncols_smc

#print('\ndsmc:')
#print('%6s %12s %12s ... %12s' % ('key', '0', '1', ncols-1))
#for k in dsmc.keys():
#    print('%6d %12.5f %12.5f ... %12.5f' % (k,dsmc[k][0],dsmc[k][1],dsmc[k][-1]))
#    
#print('\ndbeast:')
#print('%6s %12s %12s ... %12s' % ('key', '0', '1', ncols-1))
#for k in dbeast.keys():
#    print('%6d %12.5f %12.5f ... %12.5f' % (k,dbeast[k][0],dbeast[k][1],dbeast[k][-1]))

# Calculate overall mean
M = 0.0
ncheck = 0
for i in range(n):
    ncheck += len(dsmc[i])
    M += sum(dsmc[i])
    ncheck += len(dbeast[i])
    M += sum(dbeast[i])
assert ncheck == a*n*p
M /= a*n*p

# Calculate total sum-of-squares SST
SST = 0.0
for i in range(n):
    SST += sum([pow(x - M, 2.) for x in dsmc[i]])
    SST += sum([pow(x - M, 2.) for x in dbeast[i]])
dfT = a*n*p - 1
MST = SST/dfT

#print('dfT = a n p - 1 = %d %d %d - 1 = %d' % (a, n, p, a*n*p-1))

# Calculate within-run sum-of-squares SSWS
SSWS = 0.0
meanrowcol = {}
for i in range(n):
    # Calculate mean for replicate i, treatment SMC
    m = sum(dsmc[i])/len(dsmc[i])
    meanrowcol[(i,0)] = m
    SSWS += sum([pow(x - m, 2.) for x in dsmc[i]])
for i in range(n):
    # Calculate mean for replicate i, treatment BEAST
    m = sum(dbeast[i])/len(dbeast[i])
    meanrowcol[(i,1)] = m
    SSWS += sum([pow(x - m, 2.) for x in dbeast[i]])
dfWS = a*n*(p-1)    
MSWS = SSWS/dfWS

#print('dfWS = a n (p - 1) = %d %d (%d - 1) = %d' % (a, n, p, a*n*(p-1)))
 
# Calculate treatment sum-of-squares SSA
meancol = {}
meancol[0] = 0.0
meancol[1] = 0.0
for i in range(n):
    meancol[0] += sum(dsmc[i])    
    meancol[1] += sum(dbeast[i])
meancol[0] /= n*p
meancol[1] /= n*p
SSA = n*p*pow(meancol[0] - M, 2.) + n*p*pow(meancol[1] - M, 2.)
dfA = a - 1
MSA = SSA/dfA

#print('dfA = a - 1 = %d - 1 = %d' % (a,a-1))

# Calculate subject sum-of-squares SSS
SSS = 0.0
meanrow = {}
for i in range(n):
    meanrow[i] = sum(dsmc[i]) + sum(dbeast[i])
    meanrow[i] /= p*a
    SSS += p*a*pow(meanrow[i] - M, 2.)
dfS = n - 1
MSS = SSS/dfS

#print('dfS = n - 1 = %d - 1 = %d' % (n,n-1))

# Calculate treatment X subject sum-of-squares SSAxS
SSAxS = 0.0
for i in range(n):
    for j in range(2):
        SSAxS += p*pow(meanrowcol[i,j] - meanrow[i] - meancol[j] + M, 2.)
dfAxS = (a-1)*(n-1)
MSAxS = SSAxS/dfAxS

#print('dfAxS = (a - 1) (n - 1) = (%d - 1) (%d - 1) = %d' % (a,n,(a-1)*(n-1)))

FA = MSA/MSAxS
if usescipy:
    pA = 1.0 - f.cdf(FA, dfA, dfAxS)
if verbose:
    print('\nTest of whether SMC and BEAST differ:')
    print('  MSA   = %.5f' % MSA)
    print('  MSAxS = %.5f' % MSAxS)
    print('  F     = %.5f' % FA)
    print('  dfA   = %.5f' % dfA)
    print('  dfAxS = %.5f' % dfAxS)
    if usescipy:
        pA = 1.0 - f.cdf(FA, dfA, dfAxS)
        print('  P-value = %.5f' % pA)
        if pA < 0.05:
            print('  Result: significant difference between SMC and BEAST at the 5% level')
        else:
            print('  Result: no significant difference between SMC and BEAST at the 5% level')
    else:
        print('In R, issue this command to compute P value:')
        print('1 - pf(%.5f, %g, %g)' % (FA, dfA, dfAxS))

FS = MSS/MSAxS
if usescipy:
    pS = 1.0 - f.cdf(FS, dfS, dfAxS)
if verbose:
    print('\nTest of whether repeated measures was needed:')
    print('  MSS   = %.5f' % MSS)
    print('  MSAxS = %.5f' % MSAxS)
    print('  F     = %.5f' % FS)
    print('  dfS   = %.5f' % dfS)
    print('  dfAxS = %.5f' % dfAxS)
    if usescipy:
        pS = 1.0 - f.cdf(FS, dfS, dfAxS)
        print('  P-value = %.5f' % pS)
        if pS < 0.05:
            print('  Result: repeated measures was useful')
        else:
            print('  Result: repeated measures was not necessary')
    else:
        print('In R, issue this command to compute P value:')
        print('1 - pf(%.5f, %g, %g)' % (FS, dfS, dfAxS))

if anova_table:    
    if usescipy:
        print('\nRepeated-measures ANOVA (including within-run variation)\n')
        print('%12s %12s %12s %12s %12s %12s'         % ('Source',       'SS',           'df',           'MS',           'F',            'P-value'     ))
        print('%12s %12s %12s %12s %12s %12s'         % ('------------', '------------', '------------', '------------', '------------', '------------'))
        print('%12s %12.1f %12d %12.1f %12.2f %12.3f' % ('A',             SSA,            dfA,            MSA,            FA,             pA           ))
        print('%12s %12.1f %12d %12.1f %12.2f %12.3f' % ('S',             SSS,            dfS,            MSS,            FS,             pS           ))
        print('%12s %12.1f %12d %12.1f'               % ('AxS',           SSAxS,          dfAxS,          MSAxS                                        ))
        print('%12s %12.1f %12d %12.1f'               % ('WS',            SSWS,           dfWS,           MSWS                                         ))
        print('%12s %12.1f %12d %12.1f'               % ('Total',         SST,            dfT,            MST                                          ))
        print('%12s %12s %12s %12s %12s %12s'         % ('------------', '------------', '------------', '------------', '------------', '------------'))
    else:    
        print('\nRepeated-measures ANOVA (including within-run variation)\n')
        print('%12s %12s %12s %12s %12s'       % ('Source',       'SS',           'df',           'MS',           'F'))
        print('%12s %12s %12s %12s %12s'       % ('------------', '------------', '------------', '------------', '------------'))
        print('%12s %12.1f %12d %12.1f %12.2f' % ('A',             SSA,            dfA,            MSA,            FA           ))
        print('%12s %12.1f %12d %12.1f %12.2f' % ('S',             SSS,            dfS,            MSS,            FS           ))
        print('%12s %12.1f %12d %12.1f'        % ('AxS',           SSAxS,          dfAxS,          MSAxS                        ))
        print('%12s %12.1f %12d %12.1f'        % ('WS',            SSWS,           dfWS,           MSWS                         ))
        print('%12s %12.1f %12d %12.1f'        % ('Total',         SST,            dfT,            MST                          ))
        print('%12s %12s %12s %12s %12s'       % ('------------', '------------', '------------', '------------', '------------'))
