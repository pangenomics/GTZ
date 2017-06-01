#!/usr/bin/env python

import os
import sys
import argparse
from collections import defaultdict
import pandas

#HOT233_1_0025m	C01
#dictData = parseHOT_Samples_casts(infileName)
def parseHOT_Samples_casts(infileName):
	from collections import defaultdict
	#dictData = defaultdict(list)
	infile = open(infileName, "r")

	dictData = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
	dictData_samples = defaultdict(tuple)
	
	for strInLine in infile:
		strInLine = strInLine.rstrip('\n')
		arrInLine = strInLine.split('\t')

		sampleName = arrInLine[0]
		cast = int(arrInLine[1])
		bottle = int(arrInLine[2])
		
		arrSamplesName = sampleName.split('_')
		
		sampleCruise = arrSamplesName[0]
		#sampleInternalNumber = arrSamplesName[1]
		sampleDepth = arrSamplesName[-1]
		sampleDepth = sampleDepth.replace('m', '')
		sampleDepth = int(sampleDepth)
		
		dictData[sampleCruise][cast][bottle].append((sampleName, sampleDepth))
		dictData_samples[sampleName] = (sampleCruise, cast, bottle, sampleDepth)
#		dictData[samplesName]["cast"] = cast
#		dictData[samplesName]["cruise"] = sampleCruise
#		dictData[samplesName]["number"] = sampleInternalNumber
#		dictData[samplesName]["depth"] = sampleDepth

	return(dictData, dictData_samples)
	
	
def parseBottleIDs(bottleID):
	
	cruiseNr = str(bottleID)[0:3]
	sampleSite = str(bottleID)[3:5]
	castNumber = str(bottleID)[5:8]
	bottleNumber = str(bottleID)[8:10]

	return(cruiseNr, sampleSite, castNumber, bottleNumber)


def parseCruiseMixedLayerDepths(cruise2mixedLayerDepthFileName):
	dictCruise2mixLayerDepth = defaultdict(float)
	
	df_cruise2mixedLayerDepth = pandas.read_csv(cruise2mixedLayerDepthFileName, sep='\t', header=0, index_col = 0)
	
	for index, row in df_cruise2mixedLayerDepth.iterrows():
		mixLayerDepth = row.loc["mean"]

		cruiseName = "HOT" + str(index)
		dictCruise2mixLayerDepth[cruiseName] = float(mixLayerDepth)
		
	return(dictCruise2mixLayerDepth)
		
		
#first make HOT-DOGS table tab delimited
def parseHOTDOGSTable(inFile, dictData, dictData_samples, dictCruise2mixLayerDepth, outputFile, missingValueFileName):
	
	maxDiffDepth = 25
	maxDiffDepth_deep = 100
	maxDiffDepth_mix = 2
	
	maxDiffDensity = 0.1
	
	missingValueFile = open(missingValueFileName, "w")
	
	dictHOTDOGS_Data = defaultdict(lambda: defaultdict(lambda: defaultdict(pandas.Series)))
	dictHOTDOGS_Data_filtered = defaultdict(lambda: defaultdict(lambda: defaultdict(pandas.Series)))
	#dictHOTDOGS_Data_filtered_bottlenumbers = defaultdict(lambda: defaultdict(list))
	dictBottleID2Cruise = defaultdict(str)
	dictCruise2BottleIDs = defaultdict(list)
	
	#HOT-DOGS data has 2 index rows
	#inDataFrame = pandas.read_csv(inFile, sep='\t', header=[0,1], index_col = 0)
	inDataFrame = pandas.read_csv(inFile, sep='\t', header=0, index_col = 0)
	tempDataFrame = pandas.read_csv(inFile, sep='\t', header=0, index_col = 0)
	
	bottleIDs = list(inDataFrame.index)
	for bottleID in bottleIDs:
		(cruiseNr, sampleSite, castNumber, bottleNumber) = parseBottleIDs(bottleID)
		cruise = "HOT" + str(cruiseNr)
		castNumber = int(castNumber)
		bottleNumber = int(bottleNumber)
		
		#dataForBottle = pandas.Series(inDataFrame.loc[bottleID])
		dictHOTDOGS_Data[cruise][castNumber][bottleNumber] = bottleID
		dictBottleID2Cruise[bottleID] = cruise
		dictCruise2BottleIDs[cruise].append(bottleID)

		
	outputDataSeriesList = []
	listOutBottleIDs = []
	for sampleName in dictData_samples:
		dictMissingValues = defaultdict(list)
		dictMissingValues_density = defaultdict(list)
		(currentCruise, castNumber, bottleNumber, sampleDepth) = dictData_samples[sampleName]
		
		if (currentCruise in dictHOTDOGS_Data and castNumber in dictHOTDOGS_Data[currentCruise] and bottleNumber in dictHOTDOGS_Data[currentCruise][castNumber]):
			bottleID = dictHOTDOGS_Data[currentCruise][castNumber][bottleNumber] 

			dataForBottle = pandas.Series(inDataFrame.loc[bottleID])
			print(bottleID)
			listOutBottleIDs.append(bottleID)
			bottlename = dataForBottle.name
			pressure = dataForBottle.press
			density = dataForBottle.sigma
			#print(density)
			
			for (index, value) in dataForBottle.iteritems():
				
				#if the value is missing save it to recover it from other casts later
				if (int(value) < -1):
					dictMissingValues[pressure].append((bottleID, index, currentCruise, castNumber, bottleNumber, density, sampleName))			
			#print(len(dictMissingValues[pressure]))
			#dictHOTDOGS_Data_filtered[currentCruise][castNumber][bottleNumber] = dataForBottle
			#for (index, value) in enumerate(dataForBottle):
				#if the value is missing save it to recover it from other casts later
			#	if (int(value) < -1):
					#this needs to be fixed
			#		dictMissingValues[pressure].append((bottlename, index, currentCruise, castNumber, bottleNumber, density))
					#dictMissingValues_density[density].append((bottlename, index, currentCruise, castNumber,bottleNumber))		
		else:
			print("Could not find proper depth at given resolution " + str(currentCruise) + " " +  str(castNumber) + " " + str(bottleNumber)) 	


		#look for missing values
		bestMissingValues = defaultdict(lambda: defaultdict(list))
		bestMissingValues_diff = defaultdict(lambda: defaultdict(lambda: 1000000))
		bestMissingValues_diff = defaultdict(lambda: defaultdict(lambda: 1000000))
		
		for	targetDepth in dictMissingValues:
			currMissingValuesList = dictMissingValues[targetDepth]
			
			setFixed = set()
			
			for valueTuple in currMissingValuesList:
				(bottleID, index, currentCruise, castNumber, bottleNumber, targetDensity, sampleName) = valueTuple
				#print(targetDensity)
				listCruiseBottleIDs = dictCruise2BottleIDs[currentCruise]
				mixLayerDepth = dictCruise2mixLayerDepth[currentCruise]
				
				foundMatch = False
				
				for cruiseBottleID in listCruiseBottleIDs:
					#print(cruiseBottleID)
					dataForBottle = pandas.Series(inDataFrame.loc[cruiseBottleID])
					pressure = dataForBottle.press
					density = dataForBottle.sigma
					currDiffFromDepth = abs(pressure - targetDepth)
					currDiffFromDensity = abs(density - targetDensity)				
					
					
					#if in mixed layer: match on depth because temperature influence on density etc
					if ( targetDepth <= mixLayerDepth):
						if (currDiffFromDepth < maxDiffDepth_mix):
							if (float(dataForBottle.loc[index]) > -1 and (currDiffFromDepth <= bestMissingValues_diff[bottlename][index])):
								value = dataForBottle.loc[index]
								tempDataFrame.set_value(bottleID, index, value)
								bestMissingValues_diff[bottlename][index] = currDiffFromDepth
								setFixed.add(index)								
								foundMatch = True
				
					#print(dataForBottle)
					elif ((currDiffFromDepth < maxDiffDepth and targetDepth < 400) or (currDiffFromDepth < maxDiffDepth_deep and targetDepth >= 400)):
						if (currDiffFromDensity < maxDiffDensity):
							if (float(dataForBottle.loc[index]) > -1 and (currDiffFromDensity <= bestMissingValues_diff[bottleID][index])):
							#	tempDataFrame.loc[bottleID, index] = dataForBottle.loc[index]
								value = dataForBottle.loc[index]
								tempDataFrame.set_value(bottleID, index, value)
								bestMissingValues_diff[bottleID][index] = currDiffFromDensity
								setFixed.add(index)
								foundMatch = True
				
				if (not foundMatch):
					writeLine = "\t".join([str(sampleName), str(currentCruise), str(targetDepth), str(targetDensity), str(index)])
					missingValueFile.write(writeLine + "\n")
					
			print(targetDepth)			
			print(len(setFixed))
			
	for sampleName in dictData_samples:
		(sampleCruise, castNumber, bottleNumber, sampleDepth) = dictData_samples[sampleName]
		bottleID = dictHOTDOGS_Data[sampleCruise][castNumber][bottleNumber] 
		#dataForEntry = dictHOTDOGS_Data_filtered[sampleCruise][castNumber][bottleNumber]
		dataForEntry = pandas.Series(tempDataFrame.loc[bottleID])
		
		helperDict = {}
		helperDict["sample"] = sampleName
		helperDict["cruise"] = sampleCruise
		helperDict["cast"] = castNumber
		helperDict["depth"] = sampleDepth
		helperDict["bottleNumber"] = bottleNumber
		helperDict["bottleID"] = dataForEntry.name
		
		seriesHelper = pandas.Series(helperDict, index=["sample", "cruise", "cast", "depth", "bottleNumber", "bottleID"])
		newDataSeries = pandas.concat([seriesHelper, dataForEntry])
		newDataSeries.name = sampleName
		outputDataSeriesList.append(newDataSeries)	
		
	outputDataFrame = pandas.DataFrame(outputDataSeriesList)
	outputDataFrame.to_csv(outputFile, sep='\t')	

	
def runMakeMetadataTable(samplesCastFile, HOTDOGS_tableFile, cruise2mixedLayerDepthFileName, outputFile, missingValueFileName):

	(dictData, dictData_samples) = parseHOT_Samples_casts(samplesCastFile)
	dictCruise2mixLayerDepth = parseCruiseMixedLayerDepths(cruise2mixedLayerDepthFileName)
	parseHOTDOGSTable(HOTDOGS_tableFile, dictData, dictData_samples, dictCruise2mixLayerDepth, outputFile, missingValueFileName)

########################
def main(argv=None):
	if(not argv):
		argv = sys.argv[1:]
	
	parser = argparse.ArgumentParser(description='This script gets metadate for certain samples from a HOT-DOGS table', add_help = True)
	parser.add_argument('samplesCastFile', action="store", help='Tab delimited file, columns "samplenames" and "casts", no header')
	parser.add_argument('HOTDOGS_tableFile', action="store", help='Infile representing a tab delimited table, includes header')
	parser.add_argument('cruise2mixedLayerDepthFileName', action="store", help='File from HOT dogs that has the mean mixed layer depth in a column called "mean"')
	parser.add_argument('outputFile', action="store", help='outfile for resulting tables')
	parser.add_argument('missingValueFileName', action="store", help='outfile for values that could not be found and should be added fro a rerun')
	parser.add_argument('--version', '-v', action='version', version='%(prog)s 1.0')		
	args = parser.parse_args()

	runMakeMetadataTable(args.samplesCastFile, args.HOTDOGS_tableFile, args.cruise2mixedLayerDepthFileName, args.outputFile, args.missingValueFileName)
	return 0		# success

	
	
if __name__ == '__main__':
	status = main()
	sys.exit(status)		
		