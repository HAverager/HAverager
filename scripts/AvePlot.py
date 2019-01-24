#!/usr/bin/env python

# Import libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from operator import truediv
import argparse
import AveUtils
from matplotlib.colors import LogNorm


def PlotMatrix(matrix, fname, vmin=-1, vmax=1, xtitle='bin number', ytitle='bin number', label='', title='',
			   log = False, ytickspos=None, yticks=None, xtickspos=None, xticks=None, cmap=None):
	plt.figure()
	lognorm = None
	if log:
		lognorm = LogNorm()

	im = plt.imshow(matrix, interpolation='none', alpha=None, vmin=vmin, vmax=vmax,
					norm=lognorm, cmap=cmap)
	plt.xlabel(xtitle)
	plt.ylabel(ytitle)
	plt.title(title)

	if ytickspos is not None:
		plt.yticks(ytickspos, yticks)

	if xtickspos is not None:
		plt.xticks(xtickspos, xticks)

	clb = plt.colorbar()
	clb.set_label(label, labelpad=-40, y=1.05, rotation=0)
	plt.savefig(fname)

def PlotHist(data, fname, xtitle=''):
	plt.figure()
	plt.hist(data)
	plt.xlabel(xtitle)
	plt.plot()
	plt.savefig(fname)
	plt.close()

def PlotHistDecompose(data, fname, histtitle='entry/bin', xtitle='', ytitle=''):
	# definitions for the axes
	left, width = 0.1, 0.65
	bottom, height = 0.1, 0.85
	bottom_h = left_h = left + width + 0.01

	rect_scatter = [left, bottom, width, height]
	rect_histx = [left, bottom_h, width, 0.2]
	rect_histy = [left_h, bottom, 0.2, height]

	plt.figure()
	axScatter = plt.axes(rect_scatter)
	axScatter.errorbar(np.arange(len(data)),
					   np.nan_to_num(data), yerr=0, ecolor='black', marker='o', ls='')
	plt.xlim(-0.5, len(data))
	# plt.ylim(-2,2)
	plt.ylabel(ytitle)
	plt.xlabel(xtitle)

	axHisty = plt.axes(rect_histy)
	axHisty.yaxis.set_major_formatter(NullFormatter())
	axHisty.hist(np.nan_to_num(data), orientation='horizontal')
	plt.xlabel(histtitle)
	plt.plot()
	plt.savefig(fname)
	plt.close()

def PlotProfile(data, stat, total, labels, fname, ratio=2, xtitle='', ytitle='', rlims=[0.1, 1.9]):
	bins = np.arange(len(data[0]))
	fig = plt.figure()
	ax = fig.add_subplot(111)

	if ratio==1:
		fig, axs = plt.subplots(2, sharex=True)

		# Loop over all data samples
		for i in range(len(data)):
			axs[0].errorbar(bins, data[i], yerr=stat, marker='o', ls='', label=labels[i])
		axs[0].legend(numpoints=1, loc=0)
		axs[0].set_ylabel('data')
		#axs[0].set_ylim(1e-5, 1e+2)
		axs[0].set_yscale('log')

		ax = axs[1]

	# plot ratio plot
	ax.fill_between(bins, 1 + (total / data[0]), 1 - (total / data[0]),
					 alpha=0.2, edgecolor='#1B2AFF', facecolor='#089FCC')


	ax.fill_between(bins, 1 + (stat / data[0]), 1 - (stat / data[0]),
					 alpha=0.2, edgecolor='#1B2ACC', facecolor='#089FFF')

	for i in range(len(data)-1):
		ax.errorbar(bins, map(truediv, data[i+1], data[0]), yerr=map(truediv, stat, data[0]),
						marker='o', ls='', label=labels[i+1])


	plt.legend(numpoints=1, loc=0)
	plt.ylim(rlims)
	plt.xlabel(xtitle)
	plt.ylabel(ytitle)
	plt.plot()
	plt.savefig(fname)
	plt.close()


def plotAll(folder):
	systpull = np.loadtxt(folder+'sys.txt', usecols=[4])
	systnames = np.loadtxt(folder+'sys.txt', dtype='str', usecols=[1])
	data = np.loadtxt(folder+'tab.dat')
	ndata = len(data)
	pulldata = np.loadtxt(folder+'chi2map.dat', usecols=[-2])

	print ('Number of data points: ', ndata)
	print ('Number of measurements: ', len(pulldata)/ndata)
	print ('Number of systematics: ', len(systnames))
	print ('List of systematics: ', systnames)


	# plotting

	# hist pulls of systematics
	PlotHist(systpull, 'AvPullHist.pdf', 'Systematic pulls')
	PlotHistDecompose(systpull, 'AvPullHist2.pdf',
					  ytitle='pull of systematics', xtitle='index of systematic source')

	# Correlation matrix of systematic uncertainties
	Corr = AveUtils.GetCorrMatrix(np.loadtxt(folder+'tab.dat')[:,1:])
	PlotMatrix(Corr, 'Corr.pdf')

	# hist pulls of data
	for i in range(len(pulldata)/ndata):
		cutpulldata = pulldata[ndata*i:ndata*(i+1)]
		PlotHist(cutpulldata, 'AvDataPullHist_'+str(i+1)+'.pdf', 'Data pulls')
		PlotHistDecompose(cutpulldata, 'AvDataPullHist2_'+str(i+1)+'.pdf',
						  ytitle='pull of data', xtitle='index of measured point')

def run():

	########################
	# setup input parameters
	########################
	parser = argparse.ArgumentParser(description="Plot results of the averager stored in the given folder.")
	parser.add_argument('folder', metavar='folder name', type=str,
                    help='folder path to result of averager.')

	arguments = parser.parse_args()

	# Number of data points, measurements, systematic
	plotAll(arguments.folder+'/')


if __name__ == "__main__":
    run()
