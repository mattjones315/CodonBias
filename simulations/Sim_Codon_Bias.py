from __future__ import print_function, division
import argparse 
import os
import scipy.stats as st
import numpy as np
import matplotlib.pyplot as plt
import pickle as pic
import statsmodels.stats.multitest as smm
import pandas as pd


class ContextDistribution:

	def __init__(self, mut_codon, codons, kmer_len=7, _type="greedy"):

		## for now just assume that kmer length is 
		self.codons = codons
		self.mut_codon = mut_codon
		self.kmer_len = kmer_len
		self.type = _type
		self.expected_probs = np.zeros(1)
		self.context_codon = ""
		self.null_expected = np.zeros(1)
		self.null_freqs = np.zeros(1)
		self.null_observed = np.zeros(1)

	def createDistribution(self):
		#print(codons)

		self.createContext()

		if self.type == "greedy":

			self.drawProbsGreedy()


	def createContext(self):

		bases = np.array(["A", "C", "G", "T"])
		
		context = np.random.choice(bases, 4)

		# Create new codon of cBBBccc
		self.context_codon = context[0] + self.mut_codon + context[1] + context[2] + context[3]

	def find_mutation(self, c1, c2):

		muts = []
		for i in range(len(c1)):
			if c1[i] != c2[i]:
				return i

	def drawProbsGreedy(self):
		'''
		Simple probability model, where the GC content of the resulting context gives an 
		added fold increase in probability of mutation of 1.2. Every G/C awards this increase
		in probability. 

		We scale at the end by dividing each probability by the total probability of a mutation occurring
		from a given context to get probabilities summing to 1. 
		'''

		C = len(self.codons)

		# initialize the transition matrix
		probs = pd.DataFrame(np.ones((C)), index=self.codons).T

		# Add to probabilities according to GC content
		for j in np.arange(C):
			_to = probs.columns[j]
			mut_ii = self.find_mutation(self.mut_codon, _to)
			mut = _to[mut_ii]
			factor = (1 + self.context_codon.count(mut))
			probs.iloc[0,j] += factor

		## Now normalize probabilities
		probs.iloc[0,:] /= np.sum(probs.iloc[0,:])

		self.expected_probs = probs.iloc[0,:]

	def simulateNull(self, N=1000):
		C = len(self.codons)

		self.null_expected = pd.DataFrame(np.zeros((C)), index=self.codons).iloc[:,0]
		self.null_observed = pd.DataFrame(np.zeros((C)), index=self.codons).iloc[:,0]

		for i in range(N):
			self.createDistribution()
			self.null_expected += self.expected_probs
			mut_i = np.random.choice(self.codons, p=self.expected_probs)
			self.null_observed[mut_i] += 1

		freqs = self.null_observed / N 
		self.null_freqs = freqs


class Statistics:

	def __init__(self, counts, freqs, test, expected=None):

		self.counts = counts
		self.freqs = freqs
		self.test = test
		self.expected = expected

	def significance(self, N, cdist=None):

		if self.test == "simple":

			return self.singleSignificance(N)

		elif self.test == "mult":

			return self.multipleSignificance(N)

		elif self.test == "context":

			return self.contextSignificance(N, cdist)

		else:

			return 1.0

	def singleSignificance(self, N):
		'''
		Runs a Chi-squared test under the most basic null - expecting equal probability transitioning a single
		codon to all other possible codons.
		'''

		C = len(self.counts)
		expected = np.repeat(float(sum(self.counts)) / len(self.counts), len(self.counts))

		chisq = sum(((expected  - self.counts)**2 / expected))
		pv = 1 - st.chi2.cdf(chisq, C-1)

		return pv, chisq

	def multipleSignificance(self, N):
		'''
		Runs a Chi Squared test under a basic null, expecting uniform transition probability from each codon to 
		all other codons.
		'''

		C = len(self.counts)
		codon_muts = np.sum(self.counts, 1)
		expected = np.full((C,C), codon_muts / (C-1), dtype=float)
		np.fill_diagonal(expected, 0)

		chisq = ((expected - self.counts)**2 / expected) 
		chisq = np.sum(np.nan_to_num(chisq), 1)
		pv = 1 - st.chi2.cdf(chisq, C-2)
		return pv, chisq

	def contextSignificance(self, N, cdist):

		## Create null transition matrix
		C = len(cdist.codons)
		print(cdist.null_expected)
		print(cdist.null_observed)
		chisq = np.sum((cdist.null_expected - cdist.null_observed)**2 / cdist.null_expected)
		pv = 1 - st.chi2.cdf(chisq, C-1)
		return pv, chisq

def parseArgs():
	'''
	Defines command line arguments for the codon bias simulation
	'''

	parser = argparse.ArgumentParser(prog="SimCodonBias", description = "Simulate Codon Bias")

	parser.add_argument("--input", help="Input Codon File")
	parser.add_argument("--context", help="Run simulation with input file", action="store_true")
	parser.add_argument("--single", help="Run Basic Simulation", action="store_true")
	parser.add_argument("--mult", help="Run Simulation with Many Mutation Events", action="store_true")
	parser.add_argument("-n", help="Number of mutation events to simulate", default=1e8)
	parser.add_argument("--bias", help="Amount of bias to simulate, ranging from 0 to 1", default=0)
	parser.add_argument("-m", "--mut_rate", help="Mutation rate", default=1.1e-8)
	parser.add_argument("-c", "--num_codons", help="Number of codons to simulate", default=4)
	args = parser.parse_args()
	args = vars(args)

	return args

def runSingleSim(N=1e8, bias=0, num_codons=4, mut_rate=1.1e-8):
	'''
	Run Simple Simulation, where for now we simulate a single amino acide which is coded by NUM_CODONS.

	We'll draw a vector of different transition probabilities, corresponding to the relative probability of
	a particular codon transitioning to another synonymous mutation. BIAS controls the amount of bias 
	introduced to the system, manifested as a higher probability of mutation to a randomly selected codon.

	We'll simulate over N generations and then return the frequency table.
	'''

	counts = np.zeros(num_codons-1)
	bcodon = np.random.randint(0, num_codons-1, 1)

	bias_prob = float((1 + bias*(num_codons - 2))) / (num_codons-1)
	probs = np.repeat((1.0 - bias_prob)/(num_codons-2), num_codons-1, axis=0)
	probs[bcodon] = bias_prob

	for i in range(int(N)):
		mut_i = np.random.choice(np.array(range(num_codons-1)), p=probs)
		counts[mut_i] += 1

	freqs = counts / sum(counts)
	return counts, freqs


def runMultSim(N=1e8, bias=0, C=4, mut_rate = 1.1e-8):
	'''
	Run Mult simulation, where we have C codons which can all mutate to one another with some 
	probability influenced by BIAS.  

	We'll randomly draw a vector different transition probabilities, corresponding to the mutation rate from
	a single codon to the other NUM_CODONS-1 codons. Bias is a float number, between 0 and 1, which quantifies
	the degree of bias to introduce when drawing these transition probabilities. 0 corresponds to 0 bias (i.e. uniform
	probability of transition to each codon) and 1 represents only complete bias (i.e. only transition to one particular
	codon). The codon to be treated as the "biased" codon will be chosen randomly. 

	We'll then simulate over N mutation events, and then return the frequency table. We'll assume that all codons
	are present in the population with equal frequency at the beginning (equal to N), and each codon will be allowed
	to mutate according to the transition probabilities. 
	'''

	counts = np.zeros((C,C)) # instantate count matrix
	bcodon = np.random.randint(0, C, 1) # choose random codon to be biased

	# Create probability vectors for each codon

	if C == 2:
		probs = np.ones((C, C))
		probs[1,1], probs[0,0] = 0.0, 0.0

	else:
		probs = np.zeros((C, C))
		bias_prob = float((1 + bias*(C - 1))) / (C - 1)
		unif_prob = float(1 / (C-1))
		for c in np.arange(C):

			if bcodon != c:
				p = np.repeat((1.0 - bias_prob)/(C-2), C, axis=0)
				p[bcodon] = bias_prob
			else:
				p = np.repeat(unif_prob, C, axis=0)

			p[c] = 0
			probs[c,:] = p


	for c in range(C):
		for i in range(int(N)):
			mut_i = np.random.choice(C, p=probs[c,:])
			counts[c,mut_i] += 1

	freqs = np.apply_along_axis(lambda x: x / np.sum(x), 1, counts)
	return counts, freqs


def runContextSim(names, _freqs, N=1e8):
	'''
	Run a Context Simulation, where we have already read in a set of codon names and frequencies (which are
	treated as probabilites). 

	For N iterations, we'll compare the results of mutating to a codon proportional to the probabilities
	read in and the null model, here defined in relation to a mutated codon's 7 nucleotide context. 
	'''

	C = len(names)
	counts = np.zeros((C,C))

	## Create transition matrix
	probs = np.tile(_freqs, (C, 1))
	np.fill_diagonal(probs,0)
	# renormalize probabilities
	probs = np.apply_along_axis(lambda x: x / np.sum(x), 1, probs)

	counts = pd.DataFrame(counts, index=names, columns=names)

	freqs = np.apply_along_axis(lambda x: x / np.sum(x), 1, counts)
	return counts, freqs


def parseInputFile(inp):

	
	data = pd.read_csv(inp, delimiter="\t", index_col=0)
	names = data.index
	probs = data[data.columns[0]].tolist()

	return names, probs

def plot_results(tvec, pvec, power, bias, plot):

	if plot == "chi-squared":
		plt.plot(tvec, pvec, '.')
		plt.xlabel("Chi Squared Statistic")
		plt.ylabel("P value")
		plt.title("Chi Squared Distribution, Bias = " + str(bias))
		plt.show()	

	elif plot == "p_dist":

		plt.hist(pvec)
		plt.xlabel("P value")
		plt.ylabel("Frequency")
		plt.title("P Value Distribution, Bias = " + str(bias))
		#savefig('pdist.png', bbox_inches='tight')
		plt.show()

	elif plot == "qqplot":
		opv = -1.0 * np.log10(np.sort(pvec))
		epv = -1.0 * np.log10(np.arange(1, 1+len(opv))/len(opv))
		plt.plot(epv, opv, "r.")
		plt.plot(epv, epv, "k-")
		plt.xlabel("Expected P Values")
		plt.ylabel("Observed P Values")
		plt.title("QQ-Plot of P Values, Bias = " + str(bias))
		#savefig('qq.png', bbox_inches='tight')
		plt.show()

	elif plot == "power":
		plt.plot(power.keys(), power.values())
		plt.xlabel("Number of Codons")
		plt.ylabel("Power")
		plt.title("Power Simulation for Varying C, N=1000 Mutations")
		plt.show()

	else:
		print("Type of plot not recognized")


def main():

	args = parseArgs()
	print(args)

	N = int(args["n"])
	bias = float(args["bias"])
	mut_rate = float(args["mut_rate"])
	C = int(args["num_codons"])

	if args["single"]:
		pvalues = {}
		tvalues = {}
		power = {}
		pvec = list()
		tvec = list()

		for i in range(1000):
			print(i)
			counts, freqs = runSingleSim(N=N, bias=bias, mut_rate=mut_rate, num_codons=C)
			test = Statistics(counts, freqs, "simple")
			pv, tstat = test.significance(N)
			pvec.append(pv)
			tvec.append(tstat)
		pvec = np.array(pvec)
		## Correct for mutliple hypothesis testing
		#pvec = smm.multipletests(pvec, alpha=0.05, method="fdr_bh")[1]
		tvec = np.array(tvec)
		#print(np.mean(pvec))
		#print(np.mean(tvec))

		#pvalues[b] = pvec 
		#tvalues[b] = tvec
		#power[C] = len(np.where(pvec < 0.05)[0]) / len(pvec)

		#plot_results(tvec, pvec, power, bias, "chi-squared")
		plot_results(tvec, pvec, power, bias, "p_dist")
		plot_results(tvec, pvec, power, bias, "qqplot")

		#plot_results(None, None, power, "power")



	elif args["mult"]:
		pvalues = {}
		tvalues = {}
		power = {}
		
		pvec = np.array([])
		tvec = np.array([])
		for i in np.arange(500):
			print(i)
			counts, freqs = runMultSim(N=N, bias=bias, mut_rate=mut_rate, C=C)
			test = Statistics(counts, freqs, "mult")
			pv, tstat = test.significance(N)
			pvec = np.append(pvec, pv)
			tvec = np.append(tvec, tstat)
		
		#tvalues[C] = tvec
		#np.savetxt("teststat.txt", tvec, delimiter="\t")
		#pic.dump(tvalues, open("tvalues.pkl", "wb"))
		plot_results(tvec, pvec, power, bias, "chi-squared")
		plot_results(tvec, pvec, power, bias, "p_dist")
		plot_results(tvec, pvec, power, bias, "qqplot")

	elif args["context"]:
		names, probs = parseInputFile(args["input"])

		pvalues = {}
		tvalues = {}
		power = {}

		pvec = list()
		tvec = list()
		counts = None
		freqs = None

		for i in np.arange(1000):
			print(i)
			### Arbitrarily select a codon to be mut_codon, and let the rest be possible mutations
			cdist = ContextDistribution(mut_codon=names[0], codons=names[1:], _type="greedy")
			cdist.simulateNull(N)
			#counts, freqs = runContextSim(names, probs, N=N)
			test = Statistics(counts, freqs, "context")
			pv, tstat = test.significance(N, cdist)
			pvec = np.append(pvec, pv)
			tvec = np.append(tvec, tstat)

		#np.savetxt("teststat.txt", tvec, delimiter="\t")
		plot_results(tvec, pvec, power, bias, "chi-squared")
		plot_results(tvec, pvec, power, bias, "p_dist")
		plot_results(tvec, pvec, power, bias, "qqplot")



	else:
		print("Don't recognize simulation type")
	
	#plt.plot(pvalues.keys(), pvalues.values())
	#plt.xlabel("Bias")
	#plt.ylabel("Power")
	#plt.title("Power Simulation for 4 Codons, N=2000 Mutations")
	#plt.show()

if __name__ == "__main__":
	main();

				









	

