'''
Created on 5 Feb 2011

@author: Simon Bull
'''

import sys
import os
import shutil
import argparse

import checkfastaformat
import performBLAST
import adjlistcreation
import Leafcull


def main():
    """Runs the protein culling.
    
    @param args: The command line arguments.
    @type args: A Python list.
    
    """
    
    #===========================================================================
    # Parse the user's input.
    #===========================================================================
    parser = argparse.ArgumentParser(description=('Generate a non-redundant dataset of sequences from a FASTA file of input sequences. ' +
                                                  'Please see the README for more information on how to use this program.'),
                                     epilog=('This program is designed to cull a dataset of protein sequences so that no ' +
                                             'two sequences have a sequence identity greater than the specified threshold ' +
                                             'percentage. The method used is the Leaf heuristic, which is described in PAPER. ' +
                                             'A server to perform the culling can be found at http://www.bioinf.manchester.ac.uk/leaf/.')
                                     )
    parser.add_argument('inputFile', help='The location of the input FASTA file. (Required type: %(type)s).', type=str)
    parser.add_argument('-p', '--percent', help='The maximum percent sequence identity between sequences 5 <= maxPercent < 100 must be true. (Required type: %(type)s, default value: %(default)s).',
                        metavar="maxPercent", type=float, default=20, required=False)
    parser.add_argument('-m', '--minLen', help='The maximum sequence length permissible. A negative value means not to use a minimum sequence length. Must not be greater than the maximum sequence length. (Required type: %(type)s, default value: Not Used).',
                        metavar="minLength", type=int, required=False, default=-1)
    parser.add_argument('-a', '--maxLen', help='The minimum sequence length permissible A negative value means not to use a maximum sequence length. Must not be less than the minimum sequence length. (Required type: %(type)s, default value: Not Used).',
                        metavar="maxLength", type=int, required=False, default=-1)
    parser.add_argument('-s', '--seg', help='Whether or not to run SEG prior to BLASTing. (Default value: Not used).',
                        action='store_true', default=False, required=False)
    parser.add_argument('-c', '--cores', help='The number of processor cores to use for BLASTing. (Required type: %(type)s, default value: %(default)s).',
                        metavar="cores", type=int, default=2, required=False)
    parser.add_argument('-o', '--output', help='The name of the output directory to create in the current working directory. (Required type: %(type)s, default value: %(default)s).',
                        metavar="outputFolder", type=str, default='CullResults', required=False)
    parser.add_argument('-v', '--verbose', help='Whether status updates should be displayed. (Default value: No status updates).',
                        action='store_true', default=False, required=False)
    args = parser.parse_args()
    
    inputFile = args.inputFile
    sequenceIdentity = args.percent
    minLength = args.minLen
    maxLength = args.maxLen
    SEG = args.seg
    cores = args.cores
    cullOperationID = args.output
    verboseOutput = args.verbose

    #===========================================================================
    # Validate the user's input.
    #===========================================================================
    toExit = False
    if not os.path.isfile(inputFile):
        print 'The location supplied for the file of input sequences is not a valid file location.'
        toExit = True

    if sequenceIdentity < 5 or sequenceIdentity >= 100:
        print 'The maximum allowable percentage sequence similarity must be no less than 5, and less than 100.'
        toExit = True
    
    if minLength < 0:
        minLength = -1
    if maxLength < 0:
        maxLength = -1
    
    if minLength > maxLength:
        print 'The minimum sequence length must be less than the maximum sequence length.'
        toExit = True
    
    if toExit:
        sys.exit()

    #===========================================================================
    # Perform the culling.
    #===========================================================================

    # Create to directory to store the output in.
    if verboseOutput:
        print 'Creating the output directory.'
    cwd = os.getcwd()
    outputLocation = cwd + '\\' + cullOperationID
    try:
        if os.path.isdir(outputLocation):
            shutil.rmtree(outputLocation)
        elif os.path.exists(outputLocation):
            os.remove(outputLocation)
        os.mkdir(outputLocation)
    except:
        print 'The output directory could not be created. Please check the location specified in  the input parameters.'
        print 'If you did not specify a location then consider changing the default output location (the variable cullOperationID)'
        sys.exit()
    
    # Ensure that the FASTA file input is appropriately formatted.
    if verboseOutput:
        print 'Validating the input file.'
    fileToBLAST = outputLocation + '\\InputCopy.fasta'
    inputFileToLoad = open(inputFile, 'r')
    inputFile = inputFileToLoad.read()
    inputFileToLoad.close()
    errorCode, message = checkfastaformat.main(inputFile, minLength, maxLength)
    if errorCode != 0:
        print message
        sys.exit()
    writeOut = open(fileToBLAST, 'w')
    writeOut.write(message)
    writeOut.close()
    
    # Perform the BLASTing.
    if verboseOutput:
        print 'Starting the BLASTing.'
    similarities = performBLAST.main(fileToBLAST, fileToBLAST, cullOperationID + '\\BLASTOutput', SEG, cores, verboseOutput=verboseOutput)
    
    # Create the sparsematrix of the protein similarity graph.
    if verboseOutput:
        print 'Creating the adjacency matrix'
    adjList, protNames = adjlistcreation.user_seq_main(similarities, sequenceIdentity)
    
    # Choose which proteins to remove from the similarity graph.
    if verboseOutput:
        print 'Performing the culling.'
    if protNames == []:
        # This is True if there are no similarities greater than the given percentage sequence identity. If there are no
        # proteins that are too similar, then there is no need to cull any proteins from the network.
        proteinsToCull = []
    else:
        proteinsToCull, proteinsToKeep = Leafcull.main(adjList, protNames)
    
    if verboseOutput:
        print 'Writing out the results.'
    
    # Write out the proteins that were removed.
    writeOutRem = open(outputLocation + '\\Removed.txt', 'w')
    for i in proteinsToCull:
        writeOutRem.write(i + '\n')
    writeOutRem.close()
    
    # Write out a FASTA file of the proteins kept.
    writeOutKeepFasta = open(outputLocation + '\\KeptFasta.fasta', 'w')
    writeOutKeepList = open(outputLocation + '\\KeptList.txt', 'w')
    writeOutKeepList.write('IDs\tLength\n')
    readFasta = open(fileToBLAST, 'r')
    recording = False
    uniqueProteins = []  # Used to ensure no duplicates get through.
    for line in readFasta:
        if line[0] == '>' and not line[1:-1] in proteinsToCull and not line in uniqueProteins:
            # If the line starts a new protein definition, and that protein is one of the ones to keep.
            recording = True
            uniqueProteins.append(line)
            writeOutKeepFasta.write(line)
            writeOutKeepList.write(line[1:-1])
        elif line[0] == '>':
            # If the line start a new protein definition, but the protein is not one of the ones to keep.
            recording = False
        else:
            # Otherwise the line is a protein sequence.
            if recording:
                # If we are currently working on a protein that is being kept.
                writeOutKeepFasta.write(line)
                writeOutKeepList.write('\t' + str(len(line[:-1])) + '\n')
    readFasta.close()
    writeOutKeepList.close()
    writeOutKeepFasta.close()

if __name__ == '__main__':
    main()