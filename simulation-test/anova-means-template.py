from math import sqrt

row_col_table = False
verbose       = False
anova_table   = True
data_to_file  = True
R_file        = True

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
n = __NLOCI__
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

# subject is simulation replicate
# treatment is either SMC or BEAST

# Calculate meanrowcol for each subject and treatment combination
meanrowcol = {}
for i in range(n):
    # Calculate mean for replicate i, treatment SMC
    meanrowcol[(i,0)] = sum(dsmc[i])/len(dsmc[i])

    # Calculate mean for replicate i, treatment BEAST
    meanrowcol[(i,1)] = sum(dbeast[i])/len(dbeast[i])
    
# Calculate meanrow for each subject
meanrow = {}
for i in range(n):
    meanrow[i] = 0.0
    for j in range(a):
        meanrow[i] += meanrowcol[(i,j)]
    meanrow[i] /= a

# Calculate meancol for each treatment
meancol = {}
for j in range(a):
    meancol[j] = 0.0
    for i in range(n):
        meancol[j] += meanrowcol[(i,j)]
    meancol[j] /= n

# Calculate overall mean
M = 0.0
for i in range(n):
    for j in range(a):
        M += meanrowcol[(i,j)]
M /= a*n

# Calculate total sum-of-squares SST
SST = 0.0
for i in range(n):
    for j in range(a):
        SST += pow(meanrowcol[(i,j)] - M, 2.)
dfT = a*n - 1
MST = SST/dfT

# print('dfT = a n - 1 = %d %d - 1 = %d' % (a, n, a*n-1))
 
# Calculate treatment sum-of-squares SSA
SSA = 0.0
for j in range(a):
    SSA += pow(meancol[j] - M, 2.)
SSA *= n
dfA = a - 1
MSA = SSA/dfA

# print('dfA = a - 1 = %d - 1 = %d' % (a,a-1))

# Calculate subject sum-of-squares SSS
SSS = 0.0
for i in range(n):
    SSS += pow(meanrow[i] - M, 2.)
SSS *= a
dfS = n - 1
MSS = SSS/dfS

# print('dfS = n - 1 = %d - 1 = %d' % (n,n-1))

# Calculate treatment X subject sum-of-squares SSAxS
SSAxS = 0.0
for i in range(n):
    for j in range(2):
        SSAxS += pow(meanrowcol[i,j] - meanrow[i] - meancol[j] + M, 2.)
dfAxS = (a-1)*(n-1)
MSAxS = SSAxS/dfAxS

# print('dfAxS = (a - 1) (n - 1) = (%d - 1) (%d - 1) = %d' % (a,n,(a-1)*(n-1)))

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
        print('\nRepeated-measures ANOVA (no within-run variation)\n')
        print('%12s %12s %12s %12s %12s %12s'         % ('Source',       'SS',           'df',           'MS',           'F',            'P-value'     ))
        print('%12s %12s %12s %12s %12s %12s'         % ('------------', '------------', '------------', '------------', '------------', '------------'))
        print('%12s %12.1f %12d %12.1f %12.2f %12.3f' % ('A',             SSA,            dfA,            MSA,            FA,             pA           ))
        print('%12s %12.1f %12d %12.1f %12.2f %12.3f' % ('S',             SSS,            dfS,            MSS,            FS,             pS           ))
        print('%12s %12.1f %12d %12.1f'               % ('AxS',           SSAxS,          dfAxS,          MSAxS                                        ))
        print('%12s %12.1f %12d %12.1f'               % ('Total',         SST,            dfT,            MST                                          ))
        print('%12s %12s %12s %12s %12s %12s'         % ('------------', '------------', '------------', '------------', '------------', '------------'))
    else:    
        print('\nRepeated-measures ANOVA (no within-run variation)\n')
        print('%12s %12s %12s %12s %12s'       % ('Source',       'SS',           'df',           'MS',           'F'))
        print('%12s %12s %12s %12s %12s'       % ('------------', '------------', '------------', '------------', '------------'))
        print('%12s %12.1f %12d %12.1f %12.2f' % ('A',             SSA,            dfA,            MSA,            FA           ))
        print('%12s %12.1f %12d %12.1f %12.2f' % ('S',             SSS,            dfS,            MSS,            FS           ))
        print('%12s %12.1f %12d %12.1f'        % ('AxS',           SSAxS,          dfAxS,          MSAxS                        ))
        print('%12s %12.1f %12d %12.1f'        % ('Total',         SST,            dfT,            MST                          ))
        print('%12s %12s %12s %12s %12s'       % ('------------', '------------', '------------', '------------', '------------'))

if data_to_file:
    outf = open('anova-means-data.txt', 'w')
    outf.write('%12s\t%12s\t%12s\t%12s\n' % ('Replicate', 'SMC', 'BEAST', 'Mean'))
    for i in range(n):
        outf.write('%12d\t%12.5f\t%12.5f\t%12.5f\n' % (i+1, meanrowcol[(i,0)], meanrowcol[(i,1)], meanrow[i]))
    outf.write('%12s\t%12.5f\t%12.5f\t%12.5f\n' % ('mean', meancol[0], meancol[1], M))
    outf.close()
    
if R_file:
    Rfile = open('anova-means.R', 'w')
    subject = []
    treatment = []
    measure = []
    for i in range(n):
        for j in range(a):
            subject.append('%d' % i)
            treatment.append('%d' % j)
            measure.append('%g' % (meanrowcol[(i,j)],))
    total = n*a
    assert len(subject) == total
    assert len(treatment) == total
    assert len(measure) == total
    Rfile.write('library(lme4)\n')
    Rfile.write('subject <- c(%s)\n' % ','.join([s for s in subject]))
    Rfile.write('treatment <- c(%s)\n' % ','.join([t for t in treatment]))
    Rfile.write('measure <- c(%s)\n' % ','.join([m for m in measure]))
    Rfile.write('df <- data.frame(subject, treatment, measure)\n')
    Rfile.write('df$subject <- factor(df$subject)\n')
    Rfile.write('df$treatment <- factor(df$treatment)\n')
    Rfile.write('#my.aov <- with(df, aov(measure ~ treatment + Error(subject/treatment), data = df))\n')
    Rfile.write('linmod <- lmer(measure ~ treatment + (1|replicate, data=df)\n')
    Rfile.close()    

    
