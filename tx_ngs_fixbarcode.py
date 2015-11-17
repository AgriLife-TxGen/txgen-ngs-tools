#!/usr/bin/env python
from Bio import SeqIO
import csv
import gzip
import os.path
import argparse


def fixOneBarcode( keyOut, folderInput, folderOutput, SampleName, P5barcode, P7barcode, Enzyme, Lane):


	# New Data
	newLane           = Lane
	newSampleName     = SampleName.replace("_","")
	newBarcode        = P7barcode + P5barcode
	newFlowCell       = newSampleName + "-" + newBarcode
	newDNASample      = newSampleName
	newLibPlate       = newSampleName
	newRow	          = "A"
	newCol            = "1"
	newLibPrepID      = newSampleName
	newLibPlateID     = newSampleName
	newEnzyme	      = "psti-mluci"
	newBarcodeWell    = newRow + newCol
	newDNAPlate	      = ""
	newSampleDNAWell  = "" 
	newGenus	      = "zea"
	newSpecies	      = "mays"
	newPedigree	      = "" 
	newPopulation	  = ""
	newSeedLot	      = ""
	newFullSampleName = newSampleName + ":" + newFlowCell + ":" + newLane + ":" + newLibPrepID
	
	
	# Input and Output Files
	fastq_input  = folderInput + SampleName + "-" + newBarcode + "_" + Lane + "_fastq.gz"
	fastq_output = folderOutput + newSampleName + "-" + P7barcode + P5barcode + "_" + Lane + "_fastq.gz"
	if os.path.isfile(fastq_input):
		fi = gzip.open(fastq_input, 'r')
		fo = gzip.open(fastq_output, 'w')
	else:
		fastq_input  = folderInput + SampleName + "-" + newBarcode + "_" + Lane + "_fastq"
		fastq_output = folderOutput + newSampleName + "-" + P7barcode + P5barcode + "_" + Lane + "_fastq"
		if os.path.isfile(fastq_input):
			fi = open(fastq_input, 'r')
			fo = open(fastq_output, 'w')
		else:
			print("File " + fastq_input + " does not exist" )

	#print("\n\nProcessing File " + fastq_input)

	# Replacement
	toReplace  = P7barcode + P5barcode
	Replacement = newBarcode + Enzyme
	toReplaceQC = "!" * len(toReplace)
	ReplacementQC = "f" * len(Replacement)

	lenStr     = len(toReplace)
	HistDiff = [0] * (lenStr+1)
	SumOk = 0;
	SumNoOk = 0;
	posiSeq = 0
	idline = "xxx"
	#maxIters = 1000000000
	#while idline and posiSeq<maxIters:
	while idline:
		idline=fi.readline()
		seq   =fi.readline()
		spacer=fi.readline()
		quals =fi.readline()
		posiSeq = posiSeq + 1;
		numDiff = 0;
		if idline and seq:
			for Posi in range(0,lenStr):
				if(seq[Posi]!=toReplace[Posi]):
					numDiff = numDiff + 1
			HistDiff[numDiff] = HistDiff[numDiff] +1
			if numDiff<=1:
				#print("Seq = " + str(posiSeq) + " - " + seq[:lenStr] + "  OK")
				seq2 = seq.replace(seq[:lenStr], Replacement,1)
				quals2 = quals.replace(toReplaceQC, ReplacementQC,1)
				SumOk = SumOk + 1
				fo.write(idline)
				fo.write(seq2)
				fo.write(spacer)
				fo.write(quals2)
			else:
				#print("Seq = " + str(posiSeq) + " - " + seq[:lenStr] + "  Bad")
				SumNoOk = SumNoOk + 1

	keyOut.write(newFlowCell + "\t" + newLane + "\t" + newBarcode + "\t" + newDNASample + "\t"+ newLibPlate + "\t" + newRow  + "\t" + newCol  + "\t" + newLibPrepID + "\t" + newLibPlateID + "\t" + newEnzyme + "\t" + newBarcodeWell + "\t" + newDNAPlate + "\t" + newSampleDNAWell + "\t" + newGenus + "\t" + newSpecies + "\t" + newPedigree + "\t" + newPopulation + "\t" + newSeedLot + "\t" + newFullSampleName + "\n")	#if outputKeyFile!="":
	#	keyOut = open(outputKeyFile, 'a')
	#	#keyOut.write(newFlowCell + "\t" + newLane + "\t" + newBarcode + "\t" + newDNASample + "\t"+ newLibPlate + "\t" + newRow  + "\t" + newCol  + "\t" + newLibPrepID + "\t" + newLibPlateID + "\t" + newEnzyme + "\t" + newBarcodeWell + "\t" + newDNAPlate + "\t" + newSampleDNAWell + "\t" + newGenus + "\t" + newSpecies + "\t" + newPedigree + "\t" + newPopulation + "\t" + newSeedLot + "\t" + newFullSampleName + "\n")
	#	keyOut.close
	#	#print(newFlowCell + "\t" + newLane + "\t" + newBarcode + "\t" + newDNASample + "\t"+ newLibPlate + "\t" + newRow  + "\t" + newCol  + "\t" + newLibPrepID + "\t" + newLibPlateID + "\t" + newEnzyme + "\t" + newBarcodeWell + "\t" + newDNAPlate + "\t" + newSampleDNAWell + "\t" + newGenus + "\t" + newSpecies + "\t" + newPedigree + "\t" + newPopulation + "\t" + newSeedLot + "\t" + newFullSampleName + "\n")
				
	print("Input        = " + fastq_input)
	print("Output       = " + fastq_output)
	#print("Output Key   = " + outputKeyFile)
	print("To replace   = " + toReplace)
	print("Replacement  = " + Replacement)
	print("To replace Q = " + toReplaceQC)
	print("Replacement Q= " + ReplacementQC)
	print("Replaced     = " + str(SumOk))	
	print("Not Replaced = " + str(SumNoOk))	
	for Posi in range(0,lenStr+1):
		if HistDiff[Posi]>0:
			print("Differences = " + str(Posi) + " -> " + str(HistDiff[Posi]))	
	print("Flow Cell    = " + newFlowCell)
	print("Lane         = " + newLane)
	return


def fixBarcodes( keyfile, folderInput, folderOutput, lengthP7, lengthP5, Enzyme):

	print "Fixing Barcodes";
		
	#folderInput = "./unzipped/"
	#folderOutput = "./processedgz/"
	folderInput = folderInput + "/"
	folderOutput = folderOutput + "/"
	#Enzyme     = "TGCAG"

	keyFileName = os.path.basename(keyfile)
	inputKeyFile = keyfile
	outputKeyFile = folderOutput + "/" + keyFileName
	#inputKeyFile = folderInput + "keyfile.txt"
	#outputKeyFile = folderOutput + "KeyfileAuto.txt"

	
	FlowCellNames = [x[0] for x in csv.reader(open(inputKeyFile,'r'),delimiter='\t')]
	Lanes         = [x[1] for x in csv.reader(open(inputKeyFile,'r'),delimiter='\t')]
	Barcodes      = [x[2] for x in csv.reader(open(inputKeyFile,'r'),delimiter='\t')]
	SampleNames   = [x[3] for x in csv.reader(open(inputKeyFile,'r'),delimiter='\t')]

	# Writes header for new Keyfile
	keyOut = open(outputKeyFile, 'w')
	keyOut.write("Flowcell	Lane	Barcode	DNASample	LibraryPlate	Row	Col	LibraryPrepID	LibraryPlateID	Enzyme	BarcodeWell	DNA_Plate	SampleDNA_Well	Genus	Species	Pedigree	Population	SeedLot	FullSampleName\n")

	print("Num of Files = " + str(len(FlowCellNames)))

	# Loops trough Files
	for Posi in range(1,len(FlowCellNames)):
		SampleName = SampleNames[Posi]
		Barcode = Barcodes[Posi]
		Lane = Lanes[Posi]
		#print(len(Barcode))
		#print(lengthP7+lengthP5)
		if len(Barcode)>=(lengthP7+lengthP5):
			P7barcode = Barcode[0:lengthP7]
			P5barcode = Barcode[lengthP7:(lengthP7+lengthP5)]
			print("\n\nSample = " + SampleName + ", Barcode = " + Barcode + ", P7 = " + P7barcode + ", P5 = " + P5barcode + ", Lane = " + Lane)
			fixOneBarcode( keyOut, folderInput, folderOutput, SampleName, P5barcode, P7barcode, Enzyme, Lane)

	keyOut.close

	
	
# example 
# ./fixbarcodegz.py -i ./zipped/ -o ./processedgz/ -k ./zipped/keyfile.txt -e TGCAG -lp7 6 -lp5 8
	
	
parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", help="Output Folder",  required=True, metavar="Output Folder")
parser.add_argument("-i", "--input", help="Input  Folder",  required=True, metavar="Input Folder")
parser.add_argument("-k", "--keyfile", help="Keyfile containing the files to process",  required=True, metavar="Keyfile Name")
parser.add_argument("-lp7", "--lengthp7", help="length of p7 barcode",  required=True, metavar="Length of P7", type=int)
parser.add_argument("-lp5", "--lengthp5", help="length of p5 barcode",  required=True, metavar="Length of P5", type=int)
parser.add_argument("-e", "--enzyme", help="Enzyme Seq.",  required=True, metavar="Enzyme Sequence")
args = parser.parse_args()
print("\n" + os.path.basename(__file__))
print("Input Folder . = " + args.input)
print("output  Folder = " + args.output)
print("Key File ..... = " + args.keyfile)
print("Length P7 .... = " + str(args.lengthp7))
print("length P5 .... = " + str(args.lengthp5))
print("Enzyme Seq. .. = " + args.enzyme)
fixBarcodes( args.keyfile, args.input, args.output, args.lengthp7, args.lengthp5, args.enzyme)
