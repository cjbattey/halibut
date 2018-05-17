import numpy;import scipy; import pandas
import sklearn;import allel;
import tensorflow;import keras;
import os;import math;import numpy as np;
import msprime
os.chdir("/Users/cj/Dropbox/halibut/")

#functions for window boundaries
def getSubWinBounds(subWinLen, totalPhysLen): # get inclusive subwin bounds
    subWinStart = 1
    subWinEnd = subWinStart + subWinLen - 1
    subWinBounds = [(subWinStart, subWinEnd)]
    numSubWins = totalPhysLen//subWinLen
    for i in range(1, numSubWins-1):
        subWinStart += subWinLen
        subWinEnd += subWinLen
        subWinBounds.append((subWinStart, subWinEnd))
    subWinStart += subWinLen
    # if our subwindows are 1 bp too short due to rounding error, the last window picks up all of the slack
    subWinEnd = totalPhysLen
    subWinBounds.append((subWinStart, subWinEnd))
    return subWinBounds
def getSnpIndicesInSubWins(subWinBounds, snpLocs):
    snpIndicesInSubWins = []
    for subWinIndex in range(len(subWinBounds)):
        snpIndicesInSubWins.append([])
    subWinIndex = 0
    for i in range(len(snpLocs)):
        while not (snpLocs[i] >= subWinBounds[subWinIndex][0] and snpLocs[i] <= subWinBounds[subWinIndex][1]):
            subWinIndex += 1
        snpIndicesInSubWins[subWinIndex].append(i)
    return snpIndicesInSubWins

#function for windowed summary stats
def getContigStats(vcf,window_size,outfile,subpops):
    #load data
    vcf=allel.read_vcf(vcf)
    snps=allel.GenotypeArray(vcf['calldata/GT'])
    positions=vcf['variants/POS']
    sample_indices=dict()
    for i in range(len(vcf['samples'])): sample_indices[vcf['samples'][i]]=i

    #prep output file
    outfile=open(str(outfile),'w')
    outfile.write('chrom\tchromStart\tchromEnd\tnumSites\tfst\ttajD1\ttajD2\tthetaW1\tthetaW2\tdxy_bw\tpi\tdfd\n')
    outfile.close()

    #get window bounds
    window_bounds=getSubWinBounds(window_size,max(positions))
    window_bound_indices=getSnpIndicesInSubWins(window_bounds,positions)
    nwindows=max(positions)//window_size - 1

    #loop over windows and print summary stats to file
    for i in range(nwindows):
        if(len(window_bound_indices[i])<10): #if <n snps in the window
            outfile=open(str(outfile),'a')
            sumstats=[vcf['variants/CHROM'][0],str(window_bounds[i][0]),str(window_bounds[i][1]),str(0),
                      "NA","NA","NA","NA","NA","NA","NA","NA"]
            sumstats='\t'.join(sumstats)+'\n'
            outfile.write(sumstats)
            outfile.close()
        else:
            window_snp_positions=positions[window_bound_indices[i]]
            window_snps=snps.subset(window_bound_indices[i])
            window_ac_all=window_snps.count_alleles()
            window_ac_subpop=window_snps.count_alleles_subpops(subpops=subpops)
            window_ac_per_ind=window_snps.to_allele_counts()

            #summary stats
            a,b,c=allel.stats.fst.weir_cockerham_fst(window_snps,[subpops['rufus'],subpops['sasin']])
            fst=np.sum(a) / (np.sum(a) + np.sum(b) + np.sum(c))
            tajD1=allel.stats.diversity.tajima_d(window_ac_subpop['rufus'])
            tajD2=allel.stats.diversity.tajima_d(window_ac_subpop['sasin'])
            thetaW1=allel.stats.diversity.watterson_theta(window_snp_positions,window_ac_subpop['rufus'])
            thetaW2=allel.stats.diversity.watterson_theta(window_snp_positions,window_ac_subpop['sasin'])
            dxy_bw=allel.stats.diversity.sequence_divergence(window_snp_positions,window_ac_subpop['rufus'],window_ac_subpop['sasin'])
            pi=allel.stats.diversity.sequence_diversity(window_snp_positions,window_ac_all)
            dfd=allel.stats.diversity.windowed_df(window_snp_positions,window_ac_subpop['rufus'],window_ac_subpop['sasin'],size=window_size)[0][0]
            # pdxy=allel.stats.distance.pairwise_dxy(window_snp_positions,window_ac_per_ind)
            # dmax=pdxy.max()
            # dmin=pdxy.min()
            # f2=allel.stats.admixture.patterson_f2(window_ac_subpop['rufus'],window_ac_subpop['sasin']) #need to drop non-biallelic sites for this

            #write a vector of summary stats to file
            outfile=open(str(outfile),'a')
            sumstats=[vcf['variants/CHROM'][0],str(window_bounds[i][0]),str(window_bounds[i][1]),str(window_snps.shape[0]),
                      str(round(fst,6)),str(round(tajD1,6)),str(round(tajD2,6)),
                      str(round(thetaW1,6)),str(round(thetaW2,6)),str(round(dxy_bw,6)),
                      str(round(pi,6)),str(round(dfd,6))]
            sumstats='\t'.join(sumstats)+'\n'
            outfile.write(sumstats)
            outfile.close()

def getSimulationStats(vcf,outfile,subpops,length,append,mig):
    #load data
    vcf=allel.read_vcf(vcf)
    window_snps=allel.GenotypeArray(vcf['calldata/GT'])
    window_snp_positions=vcf['variants/POS']
    window_ac_all=window_snps.count_alleles()
    window_ac_subpop=window_snps.count_alleles_subpops(subpops=subpops)
    window_ac_per_ind=window_snps.to_allele_counts()

    #prep output file
    if append==0:
        out=open(str(outfile),'w')
        out.write('chrom\tchromStart\tchromEnd\tnumSites\tfst\ttajD1\ttajD2\tthetaW1\tthetaW2\tdxy_bw\tpi\tdfd\tmig\n')
        out.close()

    #summary stats
    a,b,c=allel.stats.fst.weir_cockerham_fst(window_snps,[subpops['rufus'],subpops['sasin']])
    fst=np.sum(a) / (np.sum(a) + np.sum(b) + np.sum(c))
    tajD1=allel.stats.diversity.tajima_d(window_ac_subpop['rufus'])
    tajD2=allel.stats.diversity.tajima_d(window_ac_subpop['sasin'])
    thetaW1=allel.stats.diversity.watterson_theta(window_snp_positions,window_ac_subpop['rufus'])
    thetaW2=allel.stats.diversity.watterson_theta(window_snp_positions,window_ac_subpop['sasin'])
    dxy_bw=allel.stats.diversity.sequence_divergence(window_snp_positions,window_ac_subpop['rufus'],window_ac_subpop['sasin'])
    pi=allel.stats.diversity.sequence_diversity(window_snp_positions,window_ac_all)
    dfd=allel.stats.diversity.windowed_df(window_snp_positions,window_ac_subpop['rufus'],window_ac_subpop['sasin'],size=length)[0][0]
    # pdxy=allel.stats.distance.pairwise_dxy(window_snp_positions,window_ac_per_ind)
    # dmax=pdxy.max()
    # dmin=pdxy.min()
    # f2=allel.stats.admixture.patterson_f2(window_ac_subpop['rufus'],window_ac_subpop['sasin']) #need to drop non-biallelic sites for this

    #write a vector of summary stats to file
    out=open(str(outfile),'a')
    sumstats=[vcf['variants/CHROM'][0],str(1),str(length),str(window_snps.shape[0]),
              str(round(fst,6)),str(round(tajD1,6)),str(round(tajD2,6)),
              str(round(thetaW1,6)),str(round(thetaW2,6)),str(round(dxy_bw,6)),
              str(round(pi,6)),str(round(dfd,6)),str(mig)]
    sumstats='\t'.join(sumstats)+'\n'
    out.write(sumstats)
    out.close()

############# get summary stats for empirical data #################
#indices of samples in each population (ie 0-indexed columns in vcf genotype region)
subpops={'rufus':[0,1,2,3,12,13,14,15],
         'sasin':[4,5,6,7,8,9,10,11]}

getContigStats("data_to_classify/rs_chr1peak.vcf",20000,"sumstats.txt",subpops)

############# get summary stats for simulated data ################
#simulate data with and without migration, get summary stats, print to file
#TODO: keep simulations in memory, swap rows in genotypes to erase phasing
for i in range(10000): #total simulations is n iterations * 2 (mig/nomig)
    #draw parameter values from uniform distributions
    tdiv=random.uniform(1e5,6e5)
    tmig=random.uniform(1e4,2e4)
    m=random.uniform(1e-5,1e-3)
    m=[[0,m],
       [m,0]]
    Ne_sasin=random.uniform(1e4,1e5)
    Ne_rufus=random.uniform(5e4,5e6)

    #simulate with migration
    pops={'rufus':[0,1,2,3,4,5,6,7],
          'sasin':[8,9,10,11,12,13,14,15]}
    pop_configs=[msprime.PopulationConfiguration(
                 sample_size=16,initial_size=Ne_sasin),
                 msprime.PopulationConfiguration(
                 sample_size=16,initial_size=Ne_rufus),
                 ]
    dem_events=[msprime.MigrationRateChange(time=tmig,rate=0),
                msprime.MassMigration(time=tdiv,source=1,destination=0,proportion=1)]
    trees=msprime.simulate(population_configurations=pop_configs,
                           demographic_events=dem_events,
                           migration_matrix=m,
                           mutation_rate=1e-8,
                           recombination_rate=1.4e-9,
                           length=20000)
    vcf_file = open("/Users/cj/Dropbox/halibut/sim.vcf","w")
    trees.write_vcf(vcf_file,2)
    vcf_file.close()
    if i==1:
        getSimulationStats(vcf="sim.vcf",outfile="sim_sumstats.txt",subpops=pops,length=20000,append=0,mig=1)
    else:
        getSimulationStats(vcf="sim.vcf",outfile="sim_sumstats.txt",subpops=pops,length=20000,append=1,mig=1)

    #simulate without migration
    dem_events=[msprime.MassMigration(time=tdiv,source=1,destination=0,proportion=1)]
    trees=msprime.simulate(population_configurations=pop_configs,
                           demographic_events=dem_events,
                           mutation_rate=1e-8,
                           recombination_rate=1.4e-9,
                           length=20000)
    vcf_file = open("/Users/cj/Dropbox/halibut/sim.vcf","w")
    trees.write_vcf(vcf_file,2)
    vcf_file.close()
    getSimulationStats(vcf="sim.vcf",outfile="sim_sumstats.txt",subpops=pops,length=20000,append=1,mig=0)

############################ train extra trees classifier ###########################
#load summary stats
ss=np.loadtxt("sim_sumstats.txt",skiprows=1)
statnames=open("sim_sumstats.txt","r").readline().split("\t")

#get summary stats and true stats NOTE: fix this...
out=[[0 for x in range(9)]for y in range(len(ss))];migstate=[0 for x in range(len(ss))]
for i in range(ss.shape[0]):
    line=ss[i]
    for j in range(3,len(line)-1):
        out[i][j-3]=line[j]
        migstate[i]=line[12]

#will need to figure out wtf this is
param_grid_forest = {"max_depth": [3, 10, None],
              "min_samples_split": [2, 3, 10],
              "min_samples_leaf": [1, 3, 10],
              "bootstrap": [True, False],
              "criterion": ["gini", "entropy"]}

clf = ExtraTreesClassifier(n_estimators=100, random_state=0)
grid_search = GridSearchCV(clf,param_grid=param_grid_forest,cv=10,n_jobs=10)
grid_search.fit(out, migstate)
halibut=grid_search.best_estimator_



############### generate predictions for empirical data #################
ss=np.loadtxt("sumstats.txt",skiprows=1)
statnames=open("sumstats.txt","r").readline().split("\t")

#get summary stats and true stats NOTE: fix formatting for predictions...
out=[[0 for x in range(9)]for y in range(len(ss))]
for i in range(ss.shape[0]):
    line=ss[i]
    for j in range(3,len(line)-1):
        out[i][j-3]=line[j]

pred=halibut.predict(out)
pred_out=open("chr4_preds.txt","w")
pred_out.write(str(pred))
pred_out.close()

def report(grid_scores, n_top=3):
    top_scores = sorted(grid_scores, key=itemgetter(1), reverse=True)[:n_top]
    for i, score in enumerate(top_scores):
        print("Model with rank: {0}".format(i + 1))
        print("Mean validation score: {0:.3f} (std: {1:.3f})".format(
              score.mean_validation_score,
              np.std(score.cv_validation_scores)))
        print("Parameters: {0}".format(score.parameters))
        print("")

report(grid_search.grid_scores_)


#classify the empirical data
