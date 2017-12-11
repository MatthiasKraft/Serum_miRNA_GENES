import numpy as np
import pandas as pd
import matplotlib

f = path # path to DESeq normalized read counts.


def filterExprResults(df, rc):
	df = df.loc[df.mean(axis=1)>=rc,:]# keeps only rows (i.e. genes) with mean read count >= mean read count cutoff
	#debugging
	#print('cutoff is:')
	#print(rc)
	#print(len(df.index))
	return(df)
	
def ks_test(df, ids): # no need for pairs of ids, i will create them below.
	''' 
		Computes KS statistic, returns d-statistic, and pvaue for each pair of individuals.
		df is a pandas dataframe of the gene expression file
		ids is a list of all possible pairs of individual names
	'''
	global outids

	def getMin(x):
		r=0
		x=list(x)
		m=x[0]
		print(x)
		print(len(x))
		for i in range(len(x)-1):
			if x[i]<m and x[i+1]>x[i] and abs(m-x[i])>0.015:
				r=readcuts[i]
				break
		return r
	
	from scipy.stats import ks_2samp
	import pylab as plt
	
	readcuts = np.arange(0,100,1) # make array of cuttoff values
	result = []
	count = 0
	for i in ids:
		stats = []
		s1, s2 = i
		for rc in readcuts:
			d = filterExprResults(df, rc)
			d=np.log(d+1)
			x = d[s1]
			y = d[s2]
			val = ks_2samp(x, y) # val is D-statistic, and pvalue
			stats.append(val[0]) # grab only the D-statistic
			#print s1, s2, str(rc)
			#print val
		result.append(stats)
		count += 1
		if count % 100 == 0:
			print(count)
	result = pd.DataFrame(np.array(result).T)
	result=result.set_index(readcuts)
	outids =["|".join(i) for i in ids]
	result.columns=outids
	result.to_csv('Dstats.txt', index=True, header=True, sep='\t')
	tmins = result.apply(getMin) # sends one columns at a time, one column is one pair of samples with a D-statistic for each read count cutoff.
	#print 'These are the cutoff values:'
	#print (tmins)
	tmins.index = outids
	print(type(tmins))
	tmins2 = tmins[tmins != 0]
	print(tmins2)
	print(type(tmins2))

	mean=tmins2.mean()
	median=tmins2.median()
	print('Mean read cutoff value for all pairs of samples: %d' % (mean))
	print('Median read cutoff value for all pairs of samples: %d' % (median))
	tmins2.to_csv('MeanReadCutoffValues.txt', index=True, header=True, sep='\t')
	tmins.to_csv('AllMeanReadCutoffValues.txt', index=True, header=True, sep='\t')
	std=tmins.std()
	print("Dimension of result df:")
	print(result.shape)
	result2 = result.loc[:,tmins2.index]
	print("Dimension of result2 df:")
	print(result2.shape)
	result2.plot(lw=1.5,colormap='Set1',legend=False)

	plt.axvspan(mean-std/2,mean+std/2,color='g',alpha=0.3)
	plt.xlabel('read count threshold')
	plt.ylabel('KS')
	plt.savefig('KS_test.AllSamples.pdf')
	
	randcols = result2.sample(5,axis=1)
	randcols.plot(lw=1.5,colormap='Set1',legend=False)

	# Make plot readable
	plt.axvspan(mean-std/2,mean+std/2,color='g',alpha=0.3)
	plt.xlabel('read count threshold')
	plt.ylabel('KS')
	plt.savefig('KS_test.5Samples.pdf')
	
def CreateListOfAllIDPairs(names):
	'''
	Example:
	for i in itertools.combinations([1,2,3], 2):
		print i
	'''
	import itertools
	idpairs = []
	for i in itertools.combinations(names, 2):
		idpairs.append(i)
	return(idpairs)
	
def ReadNamesFromSampleInformation(f):
	FILE = open(f, "r")
	names = [name for name in FILE.readline().strip().split('\t')] # read the first header line
	print(names)
	FILE.close()
	return(names)
	
def main():
	print('Reading sample info and getting sample names')
	names = ReadNamesFromSampleInformation(f)
	
	print('Creating all pairs of sample IDs')
	ids = CreateListOfAllIDPairs(names)
	print('Number of sample IDs; %d\n Number of pairs of samples: %d' % (len(names), len(ids)))
	
	print('Reading gene expression file')
	df = pd.read_csv(f, sep='\t') # pandas dataframe of expression values
	
	print('Starting KS-Test iteration')
	ks_test(df, ids)
	
main()
