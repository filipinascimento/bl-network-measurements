#!/usr/bin/env python

import sys
import os.path
from os.path import join as PJ
import re
import json
import numpy as np
from tqdm import tqdm
import igraph as ig

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt



def calcModularity(g):
	if("Community" in g.vertex_attributes()):
		Ci = reindexList(g.vs["Community"])
	else:
		return (None,None)
	if("weight" in g.edge_attributes()):
		return None, g.modularity(Ci, weights="weight");
	else:
		return None, g.modularity(Ci, weights=None);



def calcDegree(g):
	results = np.array(g.degree(mode="ALL"))
	return results, np.average(results)


def calcInDegree(g):
	if(not g.is_directed()):
		return (None,None)
	results = np.array(g.indegree())
	return results, np.average(results)

def calcOutDegree(g):
	if(not g.is_directed()):
		return (None,None)
	results = np.array(g.outdegree())
	return results, np.average(results)

def calcStrength(g):
	if("weight" not in g.edge_attributes()):
		return (None,None)
	results = np.array(g.strength(mode="ALL", weights = "weight"))
	return results, np.average(results)

def calcInStrength(g):
	if("weight" not in g.edge_attributes() or not g.is_directed()):
		return (None,None)
	results = np.array(g.strength(mode="IN", weights = "weight"))
	return results, np.average(results)

def calcOutStrength(g):
	if("weight" not in g.edge_attributes() or not g.is_directed()):
		return (None,None)
	results = np.array(g.strength(mode="OUT", weights = "weight"))
	return results, np.average(results)

def calcClusteringCoefficient(g):
	# if("weight" in g.edge_attributes()):
	results = g.transitivity_local_undirected(weights=None)
	# else:
	# 	results = g.transitivity_local_undirected(weights="weight")
	return np.nan_to_num(results,0), np.nanmean(results)

def calcCoreness(g):
	results = np.array(g.coreness(mode="ALL"))
	return results, None

def calcMatchIndex(g):
	degree = np.array(g.degree())
	matchIndex = np.zeros(g.ecount())
	for id,e in enumerate(g.es):
		node1,node2 = e.tuple
		viz1 = g.neighbors(node1)
		viz2 = g.neighbors(node2)
		sharedNei = set(viz1) & set(viz2)
		if ((degree[node1]+degree[node2]) > 2):
			matchIndex[id] = len(sharedNei)/float(degree[node1]+degree[node2]-2)
		else:
			matchIndex[id] = 0
	meanMatchIndex = np.mean(matchIndex)
	return None, meanMatchIndex

def calcBetweenessCentrality(g):
	result = np.array(g.betweenness(directed=g.is_directed()))
	return result,np.average(result)

def calcBetweenessCentralityWeighted(g):
	if("weight" not in g.edge_attributes()):
		return (None,None)
	result = np.array(g.betweenness(weights="weight"))
	return result,np.average(result)

def calcBetweennessCentralization(G):
	vnum = G.vcount()
	if vnum < 3:
		return None,0
	denom = (vnum-1)*(vnum-2)
	temparr = [2*i/denom for i in G.betweenness()]
	max_temparr = max(temparr)
	return None,sum(max_temparr-i for i in temparr)/(vnum-1)

def calcRichClubCoefficient(graph, highest=True, scores=None, indices_only=False):
	Trc = richClubPercentage
	degree = np.array(g.degree())
	edges = np.array(g.get_edgelist())
	sourceDegree,targetDegree = degree[edges[:,0]],degree[edges[:,1]]
	dT = int(np.percentile(degree,Trc))
	indNodes = np.nonzero(degree>=dT)[0]
	indEdges = np.nonzero((sourceDegree>=dT)&(targetDegree>=dT))[0]
	if (indNodes.size>1):
		RC = 2.*indEdges.size/(indNodes.size*(indNodes.size-1))
	else:
		RC = 0
	return None,RC

def calcDegreeAssortativity(graph):
	return None,g.assortativity_degree(directed=graph.is_directed())

def calcDiameter(graph):
	if("weight" in g.edge_attributes()):
		return None,g.diameter(directed=graph.is_directed(),weights="weight")
	else:
		return None,g.diameter(directed=graph.is_directed())

def reindexList(names,returnDict=False):
	d = {ni: indi for indi, ni in enumerate(set(names))}
	numbers = [d[ni] for ni in names]
	if(returnDict):
		return numbers,d
	else:
		return numbers

def getNeighborhoods(g,mode="ALL"):
	if("weight" in g.edge_attributes()):
		return [[(e.target,e["weight"]) if e.target!=i else (e.source,e["weight"]) for e in g.es[g.incident(i,mode=mode)]] for i in range(g.vcount())]
	else:
		return [[(e.target,1) if e.target!=i else (e.source,1) for e in g.es[g.incident(i,mode=mode)]] for i in range(g.vcount())]

def calcModuleDegreeZScore(g,mode="ALL"):
	if("Community" in g.vertex_attributes()):
		Ci = reindexList(g.vs["Community"])
	else:
		return (None,None)
	neighs = getNeighborhoods(g,mode=mode)
	cneighs = [[(Ci[vertexID],weigth) for vertexID,weigth in neigh] for neigh in neighs]
	kappa = np.zeros(g.vcount())
	kappaSi = [[] for _ in range(max(Ci)+1)]
	
	for i in range(g.vcount()):
		kappa[i] = np.sum([weight for community,weight in cneighs[i] if community==Ci[i]])
		kappaSi[Ci[i]].append(kappa[i])

	avgKappaSi = np.zeros(max(Ci)+1)
	stdKappaSi = np.zeros(max(Ci)+1)

	for ci in range(len(kappaSi)):
		avgKappaSi[ci] = np.average(kappaSi[ci])
		stdKappaSi[ci] = np.std(kappaSi[ci])
	
	zmodule = np.zeros(g.vcount())
	for i in range(g.vcount()):
		ci = Ci[i]
		if(stdKappaSi[ci]>0):
			zmodule[i] = (kappa[i]-avgKappaSi[ci])/stdKappaSi[ci]
	return zmodule,None


def calcParticipationCoeff(g,mode="ALL"):
	if("Community" in g.vertex_attributes()):
		Ci = reindexList(g.vs["Community"])
	else:
		return (None,None)
	neighs = getNeighborhoods(g,mode=mode)
	cneighs = [[(Ci[vertexID],weigth) for vertexID,weigth in neigh] for neigh in neighs]
	
	if("weight" in g.edge_attributes()):
		degrees = np.array(g.strength(mode=mode,weights="weight"))
	else:
		degrees = np.array(g.degree(mode=mode))

	kappasi = np.zeros(g.vcount())
	for i in range(g.vcount()):
		nodeCommunities = set([community for community,weight in cneighs[i]])
		communityDegrees = {community:0 for community in nodeCommunities}
		for community,weight in cneighs[i]:
			communityDegrees[community]+=weight
		kappasi[i] = np.sum(np.power(list(communityDegrees.values()),2))
	
	result = 1.0-kappasi/np.power(degrees,2.0)
	result[degrees==0.0] = 0
	return result,None


measurements = {
	"Degree" : calcDegree,
	"InDegree" : calcInDegree,
	"OutDegree" : calcOutDegree,
	"Strength" : calcStrength,
	"InStrength" : calcInStrength,
	"OutStrength" : calcOutStrength,
	"ClusteringCoefficient" : calcClusteringCoefficient,
	"Coreness" : calcCoreness,
	"MatchIndex" : calcMatchIndex,
	"BetweenessCentrality" : calcBetweenessCentrality,
	"BetweenessCentralityWeighted" : calcBetweenessCentralityWeighted,
	"BetweennessCentralization" : calcBetweennessCentralization,
	"RichClubCoefficient" : calcRichClubCoefficient,
	"DegreeAssortativity" : calcDegreeAssortativity,
	"Diameter" : calcDiameter,
	"ModuleDegreeZScore" : calcModuleDegreeZScore,
	"ParticipationCoeff" : calcParticipationCoeff,
	"Modularity" : calcModularity,
}



def check_symmetric(a, rtol=1e-05, atol=1e-08):
	return np.allclose(a, a.T, rtol=rtol, atol=atol)

def isFloat(value):
	if(value is None):
		return False
	try:
		numericValue = float(value)
		return np.isfinite(numericValue)
	except ValueError:
		return False

def loadCSVMatrix(filename):
	return np.loadtxt(filename,delimiter=",")


configFilename = "config.json"
argCount = len(sys.argv)
if(argCount > 1):
		configFilename = sys.argv[1]

outputDirectory = "output"
csvOutputDirectory = PJ(outputDirectory, "csv")
figuresOutputDirectory = PJ(outputDirectory, "figures")

if(not os.path.exists(outputDirectory)):
		os.makedirs(outputDirectory)

if(not os.path.exists(csvOutputDirectory)):
		os.makedirs(csvOutputDirectory)

if(not os.path.exists(figuresOutputDirectory)):
		os.makedirs(figuresOutputDirectory)

with open(configFilename, "r") as fd:
		config = json.load(fd)

# "index": "data/index.json",
# "label": "data/label.json",
# "csv": "data/csv",
# "transform":"absolute", //"absolute" or "signed"
# "retain-weights":false,
# "threshold": "none"

indexFilename = config["index"]
labelFilename = config["label"]
CSVDirectory = config["csv"]

richClubPercentage = 90

if("richClubPercentage" in config):
	richClubPercentage = config["richClubPercentage"];

shallPlot = True

if("generatePlots" in config):
	shallPlot = config["generatePlots"];


with open(indexFilename, "r") as fd:
	indexData = json.load(fd)

with open(labelFilename, "r") as fd:
	labelData = json.load(fd)


for entry in indexData:
	entryFilename = entry["filename"]

	alreadySigned = ("separated-sign" in entry) and entry["separated-sign"]

	baseName,extension = os.path.splitext(entryFilename)

	#inputfile,aggregatorName,isNullModel
	filenames = [(entryFilename,baseName,False)]

	if(alreadySigned):
		filenames += [(baseName+"_negative%s"%(extension),baseName+"_negative",False)]

	if("null-models" in entry):
		nullCount = int(entry["null-models"])
		filenames += [(baseName+"-null_%d%s"%(i,extension),baseName,True) for i in range(nullCount)]
		if(alreadySigned):
			filenames += [(baseName+"_negative-null_%d%s"%(i,extension),baseName+"_negative",True) for i in range(nullCount)]


	hasCommunities = False
	if("community" in entry):
		hasCommunities = (entry["community"]==True)
	
	measurementEntries = set()
	if("measurements" in entry):
		measurementEntries.update(set(entry["measurements"]))
	
	if("properties" in entry):
		del entry["properties"];


	finalNodeMeasurements = {};
	nullmodelNodeMeasurements = {};

	finalNetworkMeasurements = {};
	nullmodelNetworkMeasurements = {};

	for filename,aggregateName,isNullModel in tqdm(filenames):
		adjacencyMatrix = np.abs(loadCSVMatrix(PJ(CSVDirectory, filename))) #taking the ABS temporarily
		directionMode=ig.ADJ_DIRECTED
		weights = adjacencyMatrix
		if(check_symmetric(adjacencyMatrix)):
			directionMode=ig.ADJ_UPPER
			weights = weights[np.triu_indices(weights.shape[0], k = 0)]
		g = ig.Graph.Adjacency((adjacencyMatrix != 0).tolist(), directionMode)
		weighted = False
		if(not ((weights==0) | (weights==1)).all()):
			g.es['weight'] = weights[weights != 0]
			weighted = True
		
		inputBaseName,_ = os.path.splitext(filename)
		communitiesFilePath = PJ(CSVDirectory,"%s_community.txt"%os.path.basename(inputBaseName))

		if(hasCommunities and os.path.exists(communitiesFilePath)):
			communities = []
			with open(communitiesFilePath, "r") as fd:
				for line in fd:
					communities.append(line.strip())
			g.vs["Community"] = communities
		
		weightsProperty = None
		if(weighted):
			weightsProperty = "weight"
		
		nodeProperties = {};
		networkProperties = {};

		for measurement,measurementFunction in measurements.items():
			nodePropData,networkPropData = measurementFunction(g)
			# print("%s: "%measurement,end=" ");
			if(nodePropData is not None):
				nodeProperties[measurement] = nodePropData
				# print("n%d"%len(nodeProp),end="   ");
			if(networkPropData is not None):
				networkProperties[measurement] = networkPropData
			# 	print(networkProp,end="   ");
			# print("\n",flush=True);

		if(not isNullModel):
			finalNodeMeasurements[aggregateName] = nodeProperties 
			finalNetworkMeasurements[aggregateName] = networkProperties
		else:
			if(aggregateName not in nullmodelNodeMeasurements):
				nullmodelNodeMeasurements[aggregateName] = {}
				nullmodelNetworkMeasurements[aggregateName] = {}
			for measurement,propData in nodeProperties.items():
				if(measurement not in nullmodelNodeMeasurements[aggregateName]):
					nullmodelNodeMeasurements[aggregateName][measurement] = []
				nullmodelNodeMeasurements[aggregateName][measurement].append(propData)
			for measurement,propData in networkProperties.items():
				if(measurement not in nullmodelNetworkMeasurements[aggregateName]):
					nullmodelNetworkMeasurements[aggregateName][measurement] = []
				nullmodelNetworkMeasurements[aggregateName][measurement].append(propData)

		
		outputBaseName,outputExtension = os.path.splitext(filename)

		if("Community" in g.vertex_attributes()):
			with open(PJ(csvOutputDirectory,"%s_community.txt"%os.path.basename(outputBaseName)), "w") as fd:
				for item in g.vs["Community"]:
					fd.write("%s\n"%str(item))
		
		for nodeProperty,nodePropData in nodeProperties.items():
			propFilename = PJ(csvOutputDirectory,"%s_prop_%s.txt"%(os.path.basename(inputBaseName),nodeProperty))
			np.savetxt(propFilename,nodePropData);

		if("properties" not in entry):
			entry["properties"] = list(nodeProperties.keys());

		with open(PJ(csvOutputDirectory,os.path.basename(filename)), "w") as fd:
			if(weighted):
				outputData = g.get_adjacency(attribute='weight').data
			else:
				outputData = g.get_adjacency().data
			np.savetxt(fd,outputData,delimiter=",")

	# print(finalNetworkMeasurements);
	# nullmodelNodeMeasurements = {};

	# print(nullmodelNetworkMeasurements);
	# nullmodelNetworkMeasurements = {};
	for aggregatorName in finalNetworkMeasurements:
		with open(PJ(csvOutputDirectory,"%s__measurements.csv"%aggregatorName), "w") as fd:
			fd.write("Measurement,Value,NullModels\n")
			for measurement,value in finalNetworkMeasurements[aggregatorName].items():
				fd.write("%s,%0.18g"%(measurement,value))
				if(aggregatorName in nullmodelNetworkMeasurements
					and measurement in nullmodelNetworkMeasurements[aggregatorName]):
					nullValues = nullmodelNetworkMeasurements[aggregatorName][measurement]
					fd.write(","+",".join(["%0.18g"%nullValue for nullValue in nullValues]))
					_,bins = np.histogram([value]+nullValues,bins=30)
					if(shallPlot):
						fig = plt.figure(figsize= (8,5))
						ax = plt.axes()
						ax.hist(nullValues,bins=bins,density=True,color="#888888")
						ax.hist([value],bins=bins,density=True,color="#cc1111")
						ax.set_xlabel(measurement);
						ax.set_ylabel("Density");
						fig.savefig(PJ(figuresOutputDirectory,"network_hist_%s_%s.pdf"%(aggregatorName,measurement)));
						plt.close(fig)
				fd.write("\n")
	
	if(shallPlot):
		for aggregatorName in finalNodeMeasurements:
			for measurement,values in finalNodeMeasurements[aggregatorName].items():
					nullValues = [];
					if(aggregatorName in nullmodelNodeMeasurements
						and measurement in nullmodelNodeMeasurements[aggregatorName]):
						nullValues = list(np.array(nullmodelNodeMeasurements[aggregatorName][measurement]).flatten())
					_,bins = np.histogram(list(values)+nullValues,bins=30)
					fig = plt.figure(figsize= (8,5))
					ax = plt.axes()
					if(nullValues):
						ax.hist(nullValues,bins=bins,density=True,color="#888888")
					ax.hist(values,bins=bins,density=True,color="#cc1111",alpha=0.75)
					ax.set_xlabel(measurement);
					ax.set_ylabel("Density");
					fig.savefig(PJ(figuresOutputDirectory,"nodes_hist_%s_%s.pdf"%(aggregatorName,measurement)));
					plt.close(fig)
		


with open(PJ(outputDirectory,"index.json"), "w") as fd:
	json.dump(indexData,fd)

with open(PJ(outputDirectory,"label.json"), "w") as fd:
	json.dump(labelData,fd)

