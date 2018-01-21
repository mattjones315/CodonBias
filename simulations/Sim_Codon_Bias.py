from __future__ import print_function, division
import argparse 
import os
import scipy.stats as st
import numpy as np
import matplotlib.pyplot as plt
import pickle as pic
import statsmodels.stats.multitest as smm


class Statistics:

	def __init__(self, counts, freqs, test, expected=None):

		self.counts = counts
		self.freqs = freqs
		self.test = test
		self.expected = expected

	def significance(self, N):

		if self.test == "simple":

			return self.singleSignificance()

		elif self.test == "mult":

			return self.multipleSignificance(N)

		else:

			return 1.0

	def singleSignificance(self):
		"""
		Runs a Chi squared test under the most basic null - expecting equal probability transitioning a single
		codon to all other possible codons.
		"""

		expected = np.repeat(float(sum(self.counts)) / len(self.counts), len(self.counts))

		chisq = sum(((expected  - self.counts)**2 / expected))
		
		pv = st.chisquare(self.counts, f_exp=expected)

		return pv[1], chisq

	def multipleSignificance(self, N):
		"""
		Runs a Chi Squared test under a basic null, expecting uniform transition probability from each codon to 
		all other codons.
		"""

		C = len(self.counts)
		expected = np.repeat(N, C)

		chisq = sum( ( (expected - self.counts)**2 / expected) )
		pv = st.chisquare(self.counts, f_exp=expected, ddof=1)
		print(self.counts)
		print(expected)
		return pv[1], chisq


def parseArgs():
	"""
	Defines command line arguments for the codon bias simulation
	"""

	parser = argparse.ArgumentParser(prog="SimCodonBias", description = "Simulate Codon Bias")

	parser.add_argument("--input", help="Input Codon File")
	parser.add_argument("--single", help="Run Basic Simulation", action="store_true")
	parser.add_argument("--mult", help="Run Simulation with Many Mutation Events", action="store_true")
	parser.add_argument("-n", help="Number of mutation events to simulate", default=1e8)
	parser.add_argument("--bias", help="Amount of bias to simulate, ranging from 0 to 1", default=0)
	parser.add_argument("-m", "--mut_rate", help="Mutation rate", default=1.1e-8)
	parser.add_argument("-c", "--num_codons", help="Number of codons to simulate", default=4)
	args = parser.parse_args()
	args = vars(args)

	return args



def runMultSim(N=1e8, bias=0, num_codons=4, mut_rate = 1.1e-8):
	"""
	Run Mult simulation, where we have NUM_CODONS codons which can all mutate to one another with some 
	probability influenced by BIAS.  

	We'll randomly draw a vector different transition probabilities, corresponding to the mutation rate from
	a single codon to the other NUM_CODONS-1 codons. Bias is a float number, between 0 and 1, which quantifies
	the degree of bias to introduce when drawing these transition probabilities. 0 corresponds to 0 bias (i.e. uniform
	probability of transition to each codon) and 1 represents only complete bias (i.e. only transition to one particular
	codon). The codon to be treated as the "biased" codon will be chosen randomly. 

	We'll then simulate over N mutation events, and then return the frequency table. We'll assume that all codons
	are present in the population with equal frequency at the beginning (equal to N), and each codon will be allowed
	to mutate according to the transition probabilities. 
	"""

	
	counts = np.zeros(num_codons) # instantate count matrix
	bcodon = np.random.randint(0, num_codons, 1) # choose random codon to be biased

	# Create probability vectors for each codon
	probs = np.zeros((num_codons, num_codons))
	for c in np.arange(num_codons):
		bias_prob = float((1 + bias*(num_codons - 1))) / (num_codons - 1)
		p = np.repeat((1.0 - bias_prob)/(num_codons-2), num_codons, axis=0)
		p[bcodon] = bias_prob
		p[c] = 0
		probs[c,:] = p

	for i in range(int(N)):
		for j in range(num_codons):
			mut_i = np.random.choice(num_codons, p=probs[j,:])
			counts[mut_i] += 1

	return counts

def runSingleSim(N=1e8, bias=0, num_codons=4, mut_rate=1.1e-8):
	"""
	Run Simple Simulation, where for now we simulate a single amino acide which is coded by NUM_CODONS.

	We'll draw a vector of different transition probabilities, corresponding to the relative probability of
	a particular codon transitioning to another synonymous mutation. BIAS controls the amount of bias 
	introduced to the system, manifested as a higher probability of mutation to a randomly selected codon.

	We'll simulate over N generations and then return the frequency table.
	"""

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

def runCustomSim(names, probs, N=1e8):
	"""
	Run a Custom Simulation, where we have already read in a set of codon names and probabilities. 

	For N iterations, we'll compare the results of mutating to a codon proportional to the probabilities
	read in and the null model, here defined in relation to a mutated codon's 7 nucleotide context. 
	"""

	counts = np.zeros(num_codons)
	
	return None


def parseInputFile(inp):

	probs = np.loadtxt(inp, skiprows=1, usecols=1)
	names = np.loadtxt(inp, skiprows=1, usecols=0, dtype="str")

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
		plt.show()

	elif plot == "qqplot":
		opv = -1.0 * np.log10(np.sort(pvec))
		print(np.arange(len(opv)) / len(opv))
		epv = -1.0 * np.log10(np.arange(1, 1+len(opv))/len(opv))
		plt.plot(epv, opv, "r.")
		plt.plot(epv, epv, "k-")
		plt.xlabel("Expected P Values")
		plt.ylabel("Observed P Values")
		plt.title("QQ-Plot of P Values, Bias = " + str(bias))
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
		for C in np.arange(3, 42, step=3):
			pvec = list()
			tvec = list()

			for i in range(1000):
				print(i)
				counts, freqs = runSingleSim(N=N, bias=bias, mut_rate=mut_rate, num_codons=C)
				test = Statistics(counts, freqs, "simple")
				pv, tstat = test.significance()
				pvec.append(pv)
				tvec.append(tstat)
			pvec = np.array(pvec)
			## Correct for mutliple hypothesis testing
			#pvec = smm.multipletests(pvec, alpha=0.05, method="fdr_bh")[1]
			tvec = np.array(tvec)
			print(np.mean(pvec))
			print(np.mean(tvec))

			#pvalues[b] = pvec 
			#tvalues[b] = tvec
			power[C] = len(np.where(pvec < 0.05)[0]) / len(pvec)

			#plot_results(tvec, pvec, power, bias, "chi-squared")
			#plot_results(tvec, pvec, power, bias, "p_dist")
			#plot_results(tvec, pvec, power, bias, "qqplot")

		#plot_results(None, None, power, "power")



	elif args["mult"]:
		pvalues = {}
		tvalues = {}
		power = {}
		pvec = list()
		tvec = list()
		for i in np.arange(1000):
			print(i)
			results = runMultSim(N=N, bias=bias, mut_rate=mut_rate, num_codons=C)
			test = Statistics(results, None, "mult")
			pv, tstat = test.significance(N)
			pvec.append(pv)
			tvec.append(tstat)
		
		pvec = np.array(pvec)
		tvec = np.array(tvec)

		plot_results(tvec, pvec, power, bias, "chi-squared")
		plot_results(tvec, pvec, power, bias, "p_dist")
		plot_results(tvec, pvec, power, bias, "qqplot")

	elif args["custom"]:
		names, probs = parseInputFile(args["input_file"])

		results = runCustomSim(names, probs, N=N)
		test = Statistics(results, None, "context")
		pv, tstat = test.significance()



	else:
		print("Don't recognize simulation type")


	#plt.plot(pvalues.keys(), pvalues.values())
	#plt.xlabel("Bias")
	#plt.ylabel("Power")
	#plt.title("Power Simulation for 4 Codons, N=2000 Mutations")
	#plt.show()

if __name__ == "__main__":
	main();

				









	

