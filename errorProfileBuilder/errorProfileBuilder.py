#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pprint as pp
import argparse, os, sys

"""
Builds an error profile given alignment files and reference
Uses AlignQC's modules
"""

class AlignQCParser:
  """
  Parses files .raw and .stats
  """
  def __init__(self, prefix, output):
    """
    Initializes the parser
    :param prefix: prefix of the .raw and .stats files
    :return: nothing
    """
    self.__prefix = prefix
    self.__output = output

  def computeBasicErrorProfile(self):
    """
    Computes the basic error profile
    :return: nothing
    """
    #process .stats
    with open(self.__prefix+".error.stats") as errorStatsFile:
      errorStatsVector = [line.rstrip().split() for line in errorStatsFile]
      errorStats={}
      for (key, value) in errorStatsVector:
        errorStats[key] = int(value)

    #process .raw
    with open(self.__prefix+".error.raw") as errorRawFile:
      errorRawVector = [line.rstrip().split() for line in errorRawFile]
      #skips the header
      errorRawVector=errorRawVector[1:]

    with open(self.__output+".basic", 'w') as basicErrorProfileFile:
      sys.stdout = basicErrorProfileFile
      print("mismatches {:.06f}".format(float(errorStats["MISMATCHES"])/errorStats["ALIGNMENT_BASES"]))
      print("non-homopolymer_ins {:.06f}".format(float(errorStats["COMPLETE_INSERTION"])/errorStats["ALIGNMENT_BASES"]))
      print("homopolymer_ins {:.06f}".format(float(errorStats["HOMOPOLYMER_INSERTION"])/errorStats["ALIGNMENT_BASES"]))
      print("non-homopolymer_del {:.06f}".format(float(errorStats["COMPLETE_DELETION"])/errorStats["ALIGNMENT_BASES"]))
      print("homopolymer_del {:.06f}".format(float(errorStats["HOMOPOLYMER_DELETION"])/errorStats["ALIGNMENT_BASES"]))
      sys.stdout = sys.__stdout__

    with open(self.__output+".detailed", 'w') as detailedErrorProfileFile:
      sys.stdout = detailedErrorProfileFile
      print("mismatches {:.06f}".format(float(errorStats["MISMATCHES"])/errorStats["ALIGNMENT_BASES"]))
      for (targetBase, queryBase, count, total) in errorRawVector:
        if (targetBase!="-" and queryBase!="-" and targetBase!=queryBase):
          print("{}->{} {:0.06f}".format(targetBase, queryBase, float(count)/errorStats["MISMATCHES"]))

      print("ins {:.06f}".format(float(errorStats["ANY_INSERTION"])/errorStats["ALIGNMENT_BASES"]))
      for (targetBase, queryBase, count, total) in errorRawVector:
        if (targetBase=="-" and queryBase!="-"):
          print("{} {:0.06f}".format(queryBase, float(count)/errorStats["ANY_INSERTION"]))

      print("del {:.06f}".format(float(errorStats["ANY_DELETION"])/errorStats["ALIGNMENT_BASES"]))
      for (targetBase, queryBase, count, total) in errorRawVector:
        if (targetBase!="-" and queryBase=="-"):
          print("{} {:0.06f}".format(targetBase, float(count)/errorStats["ANY_DELETION"]))

      sys.stdout = sys.__stdout__


class ErrorProfileBuilder:
  def main(self):
    #read the args
    args = self.__doInputs()

    #run alignqc to get the error profile
    #TODO: this is not well-coded
    print("Invoking {wd}/AlignQC/alignqc/bam_to_alignment_error_plot.py to compute the error profile!".format(wd=os.path.dirname(os.path.realpath(__file__))))
    commandLine = "python {wd}/AlignQC/alignqc/bam_to_alignment_error_plot.py --max_length {maxBases} -random  -r {ref} --input_index {index} -o {bamFileOnly}.error.pdf \
    --output_stats {bamFileOnly}.error.stats --output_raw {bamFileOnly}.error.raw {bam}".format(wd=os.path.dirname(os.path.realpath(__file__)), index=args.sorted_bam_index,\
    ref=args.reference, maxBases=args.max_bases, bam=args.sorted_bam, bamFileOnly=os.path.basename(args.sorted_bam))
    print(commandLine)
    os.system(commandLine)
    print("Done!")

    #parse the error profile and create the files we want
    alignQCParser = AlignQCParser(os.path.basename(args.sorted_bam), args.output)
    alignQCParser.computeBasicErrorProfile()

  def __doInputs(self):
    # Setup command line inputs
    parser=argparse.ArgumentParser(description="Builds an error profile given alignment files and reference",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    requiredArgs = parser.add_argument_group('Required arguments')
    requiredArgs.add_argument('-b', '--sorted_bam',help="sorted BAMFILE", required=True)
    requiredArgs.add_argument('-i', '--sorted_bam_index',help="sorted BAMFILE index", required=True)
    requiredArgs.add_argument('-r','--reference',help="Fasta reference file", required=True)
    parser.add_argument('-o','--output', default="error_profile", help="OUTPUTFILE for error profile")
    parser.add_argument('--max_bases',type=int,default=100000,help="Maximum number of bases to build the profile from")

    args = parser.parse_args()
    return args




if __name__ == "__main__":
  #testing
  errorProfileBuilder = ErrorProfileBuilder()
  errorProfileBuilder.main()
