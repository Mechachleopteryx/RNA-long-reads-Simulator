#!/usr/bin/env python3 
# -*- coding: utf-8 -*-

# ############################################################################
#	Pipeline for the simulation of transcriptomics data sets of ONT long reads
#	authors: Camille Marchet & Leandro Ishi
# 	contact: camille.marchet@irisa.fr
# ############################################################################

import os
import re
import sys
import time
import shlex, subprocess
import struct
import shutil
import os.path
import tempfile
import argparse
from subprocess import Popen, PIPE, STDOUT

# ############################################################################
#									Utils functions
# ############################################################################

# get the platform
def getPlatform():
	if sys.platform == "linux" or sys.platform == "linux2":
		return "linux"
	elif sys.platform == "darwin":
		return "OSX"
	else:
		print("[ERROR] BWISE is not compatible with Windows.")
		sys.exit(1);


# get the timestamp as string
def getTimestamp():
	return "[" + time.strftime("%H:%M:%S") + " " + time.strftime("%d/%m/%Y") + "] "



# check if reads files are present
def checkReadFiles(readfiles):
	if readfiles is None:
		return True
	allFilesAreOK = True
	#~ for file in readfiles:
	if not os.path.isfile(readfiles):
		print("[ERROR] File \""+file+"\" does not exist.")
		allFilesAreOK = False
	if not allFilesAreOK:
		dieToFatalError("One or more read files do not exist.")

		
# check if files written are present
def checkWrittenFiles(files):
	allFilesAreOK = True
	if not os.path.isfile(files):
		print("[ERROR] There was a problem writing \"" + files + "\".")
		allFilesAreOK = False
	if not allFilesAreOK:
		dieToFatalError("One or more files could not be written.")



# to return if an error makes the run impossible
def dieToFatalError (msg):
  print("[FATAL ERROR] " + msg)
  print("Try `Bwise --help` for more information")
  sys.exit(1);


# launch subprocess
def subprocessLauncher(cmd, argstdout=None, argstderr=None,	 argstdin=None):
	args = shlex.split(cmd)
	p = subprocess.Popen(args, stdin = argstdin, stdout = argstdout, stderr = argstderr).communicate()
	return p

# print time in hh:mm:ss
def printTime(msg, seconds):
	m, s = divmod(seconds, 60)
	h, m = divmod(m, 60)
	return msg + " %d:%02d:%02d" % (h, m, s)

# print a warning message
def printWarningMsg(msg):
	print("[Warning] " + msg)



def main():

	wholeT = time.time()
	print("\n*** This is RNACreator by Leandro & Camille - We simulate your RNA long reads! ***\n")
	SIMULATOR_MAIN = os.path.dirname(os.path.realpath(__file__))
	SIMULATOR_PATH = SIMULATOR_MAIN + "/simulator" 
	ERROR_PROFILE_BUILDER_PATH = SIMULATOR_MAIN + "/errorProfileBuilder"
	SEQ_GTF_EXTRACTOR_PATH = SIMULATOR_MAIN + "/extractSequencesFromGTF"
	
	# ========================================================================
	#						 Manage command line arguments
	# ========================================================================
	parser = argparse.ArgumentParser(description='RNACreator - Simulation of whole transcriptome data sets of ONT long reads')

	# ------------------------------------------------------------------------
	#							 Define allowed options
	# ------------------------------------------------------------------------
	parser.add_argument('-g',	action="store",	dest="gtfFilePath",	type=str,	default = "",	help="Path to the GTF file")
	parser.add_argument('-b',	action="store",	dest="bamFilePath",	type=str,	default = "",	help="Path to the BAM file - used only to produced the error profile")
	parser.add_argument('-i',	action="store",	dest="baiFilePath",	type=str,	default = "",	help="Path to the BAI file - used only to produced the error profile")
	parser.add_argument('-r',	action="store",	dest="genomeRefPath",	type=str,	default = "",	help="Path to the reference genome file")
	parser.add_argument('-c', action="store", dest="coverage",	type=str,	default = "1",	help="An integer that represents the desired coverage (default=1)")
	parser.add_argument('-o', action="store", dest="outputDirPath",	type=str,	default = ".",	help="Path to the output directory (default: .)")
	parser.add_argument('--version', action='version', version='%(prog)s 0.0.1')

	# ------------------------------------------------------------------------
	#				Parse and interpret command line arguments
	# ------------------------------------------------------------------------
	options = parser.parse_args()
		
	# ------------------------------------------------------------------------
	#				  Print command line
	# ------------------------------------------------------------------------
	print("The command line was: " + ' '.join(sys.argv))

	# ------------------------------------------------------------------------
	#				  Misc parameters
	# ------------------------------------------------------------------------
	gtfFilePath			= options.gtfFilePath
	bamFilePath			= options.bamFilePath
	baiFilePath			= options.baiFilePath
	genomeRefPath		= options.genomeRefPath
	coverage			= options.coverage
	outputDirPath 		= options.outputDirPath
	
	# ------------------------------------------------------------------------
	#				Check if mandatory arguments are missing
	# ------------------------------------------------------------------------
	if gtfFilePath == "" or bamFilePath == "" or baiFilePath == "" or genomeRefPath == "":
		dieToFatalError("There are missing arguments.\n Mandatory arguments: -b -g -i -r\nUsage: ./RNACreator CHROMOSOME_FILES_DIRECTORY_PATH -g PATH_TO_GTF -b PATH_TO_BAM -i PATH_TO_BAI -r PATH_TO_REFERENCE_GENOME -c COVERAGE -o OUTPUT_DIRECTORY") 


	# ------------------------------------------------------------------------
	#				Create output dir and log files
	# ------------------------------------------------------------------------
	
	try:
		if not os.path.exists(outputDirPath):
			os.mkdir(outputDirPath)
		else:
			printWarningMsg(outputDirPath+ " directory already exists, we will use it.")
		#log files
		OUT_LOG_FILES = outputDirPath + "/logs"
		if not os.path.exists(OUT_LOG_FILES):
			os.mkdir(OUT_LOG_FILES)
		# intermediate output files
		OUT_RESULTS_FILES = outputDirPath + "/files"
		if not os.path.exists(OUT_RESULTS_FILES):
			os.mkdir(OUT_RESULTS_FILES)
		outName = outputDirPath.split("/")[-1]
		outputDirPath = os.path.dirname(os.path.realpath(outputDirPath)) + "/" + outName
		print("Results will be stored in: ", outputDirPath)
	except:
		print("Could not write in output directory :", sys.exc_info()[0])
		dieToFatalError('')
	#TODO print cmd line and path to the different input files in log

	# ========================================================================
	#						Get error profile
	# ========================================================================
	try:
		print("Getting error profile...")
		errorProfileCmd = "python " + ERROR_PROFILE_BUILDER_PATH + "/errorProfileBuilder.py -b " + bamFilePath + " -i " + baiFilePath + " -r " + genomeRefPath
		print getTimestamp() + "Running " + errorProfileCmd
		subprocessLauncher(errorProfileCmd)
		cmdMv = "mv error_profile.basic " + OUT_RESULTS_FILES
		subprocess.check_output(['bash','-c', cmdMv])
		checkWrittenFiles(OUT_RESULTS_FILES + "/error_profile.basic")
	except SystemExit:	# happens when checkWrittenFiles() returns an error
		sys.exit(1);
	except KeyboardInterrupt:
		sys.exit(1);
	except:
		print("Unexpected error during error profile computation:", sys.exc_info()[0])
		dieToFatalError('')

	#TODO check if keys in gtf are the same than in the fasta files

	# ========================================================================
	#						Get expression levels
	# ========================================================================
	print("Getting expression levels...")
	
	try:
		# ------------------------------------------------------------------------
		#				  Write a .par file for Flux Simulator
		# ------------------------------------------------------------------------
		parameterFile = open(OUT_RESULTS_FILES + "/file_for_expression.par", 'w')
		parameterFile.write("REF_FILE_NAME\t" + gtfFilePath + "\nPRO_FILE\t"  + OUT_RESULTS_FILES + "/expression.pro")
		checkWrittenFiles(OUT_RESULTS_FILES + "/file_for_expression.par")
		# ------------------------------------------------------------------------
		#				  Launch Flux Simulator for expression generation
		# ------------------------------------------------------------------------
		expressionLevelsCmd = "flux-simulator -x -p " + + OUT_RESULTS_FILES + "/file_for_expression.par"
		checkWrittenFiles( OUT_RESULTS_FILES + "/expression.pro")
	except SystemExit:
		sys.exit(1);
	except KeyboardInterrupt:
		sys.exit(1);
	except:
		print("Unexpected error during expression levels computation:", sys.exc_info()[0])
		dieToFatalError('')

	# ========================================================================
	#						Get reference transcripts
	# ========================================================================
	try:
		print("Getting reference transcripts...")
		referenceTranscriptsCmd = SEQ_GTF_EXTRACTOR_PATH + "gffread/gffread -g " + genomeRefPath + " -w " + OUT_RESULTS_FILES + "/transcripts.fa " + gtfFilePath
		subprocessLauncher(referenceTranscriptsCmd)
		checkWrittenFiles(OUT_RESULTS_FILES + "/transcripts.fa")
	except SystemExit:
		sys.exit(1);
	except KeyboardInterrupt:
		sys.exit(1);
	except:
		print("Unexpected error during reference transcripts computation:", sys.exc_info()[0])
		dieToFatalError('')
		
	# ========================================================================
	#						Simulate long reads
	# ========================================================================
	try:
		print("Simulating reads...")
		simulationCmd = SIMULATOR_PATH + "./theReadCreator -c " + coverage + " -e " + OUT_RESULTS_FILES + "/error_profile.basic -t " + OUT_RESULTS_FILES + "/transcripts.fa -p "+ OUT_RESULTS_FILES + "/expression.pro"
		subprocessLauncher(referenceTranscriptsCmd)
		cmdMv = "mv simulatedReads.fa simulatedPerfectSequences.fa " + outputDirPath
		subprocess.check_output(['bash','-c', cmdMv])
		checkWrittenFiles(outputDirPath + "/simulatedReads.fa")
		checkWrittenFiles(outputDirPath + "/simulatedPerfectSequences.fa")
		print(printTime("\nThe end !\nSimulation took: ", time.time() - wholeT))
	except SystemExit:
		sys.exit(1);
	except KeyboardInterrupt:
		sys.exit(1);
	except:
		print("Unexpected error during read simulation process:", sys.exc_info()[0])
		dieToFatalError('')


#TODO rm other files

if __name__ == '__main__':
	main()
