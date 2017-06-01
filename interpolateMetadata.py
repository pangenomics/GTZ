#!/usr/bin/env python

import os
import sys
import argparse
from collections import defaultdict
import pandas
import numpy

def parseBottleIDs(bottleID):
	
	cruiseNr = str(bottleID)[0:3]
	sampleSite = str(bottleID)[3:5]
	castNumber = str(bottleID)[5:8]
	bottleNumber = str(bottleID)[8:10]

	return(cruiseNr, sampleSite, castNumber, bottleNumber)


def parseMissingValuesFile(missingValueFileName):
	missingValueFile = open(missingValueFileName, "r")
	#HOT226  499.7   26.611  chl
	
	dictMissingValues = defaultdict(lambda: defaultdict(list))
	
	for line in missingValueFile:
		line = line.rstrip('\n')
		arrLine = line.split('\t')

		sample = arrLine[0]
		cruise = arrLine[1]
		depth = float(arrLine[2])
		density = float(arrLine[3])
		index = arrLine[4]
		
		dictMissingValues[cruise][index].append((sample,density))
		
	return(dictMissingValues)


#first make HOT-DOGS table tab delimited
def interpolateMetadata(inFile, dictMissingValues, metaDataTableFile, outputFile):
	
	dictCruise2index2BottleIDs = defaultdict(lambda: defaultdict(list))	
	
	dictBottleID2Cruise = defaultdict(str)
	dictCruise2BottleIDs = defaultdict(list)
	
	
	
	inDataFrame = pandas.read_csv(inFile, sep='\t', header=0, index_col = 0)
	df_metaDataTable = pandas.read_csv(metaDataTableFile, sep='\t', header=0, index_col = 0)
	
	bottleIDs = list(inDataFrame.index)
	for bottleID in bottleIDs:
		(cruiseNr, sampleSite, castNumber, bottleNumber) = parseBottleIDs(bottleID)
		cruise = "HOT" + str(cruiseNr)
		castNumber = int(castNumber)
		bottleNumber = int(bottleNumber)
		
		#dataForBottle = pandas.Series(inDataFrame.loc[bottleID])
		dictBottleID2Cruise[bottleID] = cruise
		dictCruise2BottleIDs[cruise].append(bottleID)
		
	for cruise in dictMissingValues:
		for index in dictMissingValues[cruise]:
			listMissingDensities = dictMissingValues[cruise][index]
			
			dictDensity2value = {}
			listDensities = []
			
			listCruiseBottleIDs = dictCruise2BottleIDs[cruise]
			for bottleID in listCruiseBottleIDs:
				value = float(inDataFrame.loc[bottleID, index])
				density = float(inDataFrame.loc[bottleID, "sigma"])
				if (value > -1):
					dictDensity2value[density] = value
					listDensities.append(density)
			
			listDensities = sorted(listDensities)
			listValues = []
			for density in listDensities:
				value = dictDensity2value[density]
				listValues.append(value)
				
			if (len(listValues) > 0):
				for (sampleName, missingDensity) in listMissingDensities:
					left_input = -9
					right_input = -9
					print missingDensity
					print listDensities
					print listValues
					listOutValues = numpy.interp([missingDensity], listDensities, listValues, left=left_input, right=right_input)				
					value = listOutValues[0]
					df_metaDataTable.set_value(sampleName, index, value)
	
	df_metaDataTable.to_csv(outputFile, sep='\t')	

		
def runInterpolateMetadata(missingValueFileName, metaDataTableFile, HOTDOGS_tableFile, outputFile):

	dictMissingValues = parseMissingValuesFile(missingValueFileName)
	interpolateMetadata(HOTDOGS_tableFile, dictMissingValues, metaDataTableFile, outputFile)
	

########################
def main(argv=None):
	if(not argv):
		argv = sys.argv[1:]
	
	parser = argparse.ArgumentParser(description='This script gets metadate for certain samples from a HOT-DOGS table', add_help = True)
	parser.add_argument('missingValueFileName', action="store", help='Tab delimited file, columns "samplename", "cruise", "depth", "density" and "index/datatype"; no header (produced by getMetadataTable_HOTDOGS_densityLink.py)')
	parser.add_argument('metaDataTable', action="store", help='Infile representing a tab delimited table, includes header (produced by getMetadataTable_HOTDOGS_densityLink.py)')
	parser.add_argument('HOTDOGS_tableFile', action="store", help='Infile representing a tab delimited table, includes header')
	parser.add_argument('outputFile', action="store", help='outfile for resulting tables')
	parser.add_argument('--version', '-v', action='version', version='%(prog)s 1.0')		
	args = parser.parse_args()

	runInterpolateMetadata( args.missingValueFileName,  args.metaDataTable,  args.HOTDOGS_tableFile,  args.outputFile)
	return 0		# success

	
	
if __name__ == '__main__':
	status = main()
	sys.exit(status)		
				
		
		
	
		

	
		
		
		