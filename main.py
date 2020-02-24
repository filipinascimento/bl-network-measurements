#!/usr/bin/env python

import sys
import os.path
import re
import numpy as np
from tqdm import tqdm

def loadEdgesList(filename):
	edgesList = []
	weights = [] # Will be weighted if len(edgesList)==len(weights);
	with open(filename,"r") as fd:
		for line in fd:
			line = line.strip()
			if line:
				entries = [int(entry) for entry in re.split(r'\W+',line)]
				fromIndex = entries[0]
				toIndex = entries[1]
				edgesList.append((fromIndex,toIndex));
				if(len(entries)>2):
					weights.append(float(entries[2]));
				edgesList.append((fromIndex,toIndex));
	if(len(edgesList)==len(weights)):
		edgesList = [(fromIndex,toIndex,weights[index]) for index,(fromIndex,toIndex) in enumerate(edgesList)];
	return edgesList

def loadCSV(filename,threshold=0.5,directed=False):
	edgesList = []
	with open(filename,"r") as fd:
		for fromIndex,line in enumerate(fd):
			line = line.strip()
			if line:
				entries = line.split(",");
				if(not directed):
					entries=entries[:fromIndex];
				edgesList+=[(fromIndex,toIndex,float(weightText)) for toIndex,weightText in enumerate(entries) if float(weightText)>threshold]
				# edgesList.append((fromIndex,toIndex,weight))
	return edgesList

def saveEdgesListTo(edgesList,filename):
	if(edgesList and len(edgesList[0])>2):
		#weighted
		with open(filename,"w") as fd:
			fd.write("\n".join(["%s\t%s\t%g"%edge for edge in edgesList]))
	else:
		#non-weighted
		with open(filename,"w") as fd:
			fd.write("\n".join(["%s\t%s"%edge for edge in edgesList]))


networkName = "data/cm.csv"
realizations = 100
alpha = 1.0
rewireMode = "basic"

argCount = len(sys.argv)
if(argCount>1):
	networkName = sys.argv[1]
	if(argCount>2):
		realizations = int(sys.argv[2])
	if(argCount>3):
		alpha = float(sys.argv[3])
	if(argCount>4):
		alpha = float(sys.argv[4])


networkBaseName = os.path.splitext(os.path.basename(networkName))[0]
outputDirectory = "output"
if(not os.path.exists(outputDirectory)):
	os.makedirs(outputDirectory)

originalEdges = loadCSVList(networkName)
for realization in tqdm(range(realizations)):
	if(rewireMode=="configurational"):
		edgesCount = len(originalEdges)
		selectedEdges = list(np.where(np.random.random(edgesCount)<alpha)[0])
		np.random.shuffle(selectedEdges)
		newEdges = originalEdges.copy()
		if(len(selectedEdges)>2):
			for edgeIndex in selectedEdges:
				while(True):
					selectedEdge = newEdges[edgeIndex]
					crossIndex = np.random.choice(selectedEdges); #selected from selected randomly
					crossEdge = newEdges[crossIndex]

					selectedSource = min(selectedEdge[0:2])
					selectedTarget = max(crossEdge[0:2])
					crossSource = min(crossEdge[0:2])
					crossTarget = max(selectedEdge[0:2])

					if(selectedSource!=selectedTarget and crossSource!=crossTarget):
						newEdges[edgeIndex] = (selectedSource,selectedTarget)
						newEdges[crossIndex] = (crossSource,crossTarget)
						break
	else:
		edgesCount = len(originalEdges);
		newEdges = np.array(originalEdges);
		selectedEdgesIndices = np.where(np.random.random(edgesCount)<alpha)[0];
		generatedSelectedEdges = np.array(np.random.randint(0,vertexCount,(len(selectedEdgesIndices),2)));
		newEdges[selectedEdgesIndices] = generatedSelectedEdges;
		newEdges = [(fromIndex,toIndex) for fromIndex,toIndex in newEdges];
	saveEdgesListTo(newEdges,os.path.join(outputDirectory, "%s-nullA%.3f-%d.edgeslist"%(networkBaseName,alpha,realization)))
