'''
Created on 5 Feb 2011

@author: Simon Bull
'''

import argparse
import os
import shutil
import sys

import Leafcull


def main(args):
    """Runs the protein culling.

    :param args:    The command line arguments.
    :type args:     list

    """

    #===========================================================================
    # Parse the user's input.
    #===========================================================================
    parser = argparse.ArgumentParser(description=('Generate a non-redundant dataset of chains from the PDB. ' +
                                                  'Please see the README for more information on how to use this program.'),
                                     epilog=('This program is designed to cull a dataset of protein sequences so that no ' +
                                             'two sequences have a sequence identity greater than the specified threshold ' +
                                             'percentage. The method used is the Leaf heuristic, which is described in a paper located at ' +
                                             'http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0055484.' +
                                             'A server to perform the culling can be found at http://leaf-protein-culling.appspot.com/.')
                                     )
    parser.add_argument('datalocation', help='The location of the directory that contains the processed PDB data.')
    parser.add_argument('-i', '--inputFile', help='The location of the file containing the identifiers of the chains to cull. (Required type: %(type)s, default value: Not used).',
                        metavar="inputFile", type=str, default='', required=False)
    parser.add_argument('-p', '--percent', help='The maximum percent sequence identity between sequences. 5 <= maxPercent < 100 must be true. (Required type: %(type)s, default value: %(default)s).',
                        metavar="maxPercent", type=float, default=20, required=False)
    parser.add_argument('-q', '--minRes', help='The minimum resolution a chain/entry can have. 0 <= minResolution <= 100 must be true. Must not be greater than the maximum resolution. (Required type: %(type)s, default value: %(default)s).',
                        metavar="minResolution", type=float, default=0.0, required=False)
    parser.add_argument('-r', '--maxRes', help='The maximum resolution a chain/entry can have. 0 <= maxResolution <= 100 must be true. Must not be less than the minimum resolution. (Required type: %(type)s, default value: %(default)s).',
                        metavar="maxResolution", type=float, default=3.0, required=False)
    parser.add_argument('-l', '--rval', help='The maximum R value a chain/entry can have. 0 <= maxRValue <= 1 must be true. (Required type: %(type)s, default value: %(default)s).',
                        metavar="maxRValue", type=float, default=0.5, required=False)
    parser.add_argument('-m', '--minLen', help='The maximum sequence length permissible. A negative value means not to use a minimum sequence length. Must not be greater than the maximum sequence length. (Required type: %(type)s, default value: Not Used).',
                        metavar="minLength", type=int, required=False, default=-1)
    parser.add_argument('-a', '--maxLen', help='The minimum sequence length permissible A negative value means not to use a maximum sequence length. Must not be less than the minimum sequence length. (Required type: %(type)s, default value: Not Used).',
                        metavar="maxLength", type=int, required=False, default=-1)
    parser.add_argument('-x', '--usenonxray', help='Whether non-X-ray chains/entries should be included. (Default value: Not used).',
                        action='store_true', default=False, required=False)
    parser.add_argument('-u', '--usealpha', help='Whether alpha carbon only chains/entries should be included. (Default value: Not used).',
                        action='store_true', default=False, required=False)
    parser.add_argument('-o', '--output', help='The name of the output directory to create in the current working directory. (Required type: %(type)s, default value: a directory called %(default)s in the current working directory).',
                        metavar="outputFolder", type=str, default='PDBCullResults', required=False)
    parser.add_argument('-v', '--verbose', help='Whether status updates should be displayed. (Default value: No status updates).',
                        action='store_true', default=False, required=False)
    args = parser.parse_args()

    dataDirectory = args.datalocation
    fileUserInputChains = args.inputFile
    sequenceIdentity = args.percent
    minResolution = args.minRes
    maxResolution = args.maxRes
    maxRValue = args.rval
    minLength = args.minLen
    maxLength = args.maxLen
    skipNonXray = not args.usenonxray
    skipAlphaCarbon = not args.usealpha
    cullOperationID = args.output
    verboseOutput = args.verbose

    #===========================================================================
    # Validate the user's input.
    #===========================================================================
    toExit = False
    if not os.path.isdir(dataDirectory):
        print('The data directory supplied is not a valid directory.')
        toExit = True

    cullSubset = fileUserInputChains != ''  # Whether the user specified a subset of the chains should be culled.
    if cullSubset:
        if not os.path.isfile(fileUserInputChains):
            print('The file of chains to cull is not a valid file location.')
            toExit = True

    if sequenceIdentity < 5 or sequenceIdentity >= 100:
        print('The maximum allowable percentage sequence similarity must be no less than 5, and less than 100.')
        toExit = True

    if minResolution < 0 or minResolution > 100:
        print('The valid range for the minimum resolution is 0 - 100.')
        toExit = True

    if maxResolution < 0 or maxResolution > 100:
        print('The valid range for the maximum resolution is 0 - 100.')
        toExit = True

    if minResolution > maxResolution:
        print('The minimum resolution must be less than or equal to the maximum resolution.')
        toExit = True

    if maxRValue < 0 or maxRValue > 1:
        print('The valid range for the maximum R value is 0 - 1.')
        toExit = True

    if minLength < 0:
        minLength = -1
    if maxLength < 0:
        maxLength = -1

    if minLength > maxLength:
        print('The minimum sequence length must be less than the maximum sequence length.')
        toExit = True

    if toExit:
        sys.exit()

    # Create the directory to store the output in.
    cwd = os.getcwd()
    if cullOperationID == 'PDBCullResults':
        outputLocation = cwd + '/' + cullOperationID
    else:
        outputLocation = cullOperationID
    try:
        if os.path.isdir(outputLocation):
            shutil.rmtree(outputLocation)
        elif os.path.exists(outputLocation):
            os.remove(outputLocation)
        os.mkdir(outputLocation)
    except:
        print('The output directory could not be created. Please check the location specified in  the input parameters.')
        print('If you did not specify a location then consider changing the default output location (the variable cullOperationID)')
        sys.exit()

    #===========================================================================
    # PDB data files.
    #===========================================================================
    similarityData = dataDirectory + '/' + 'Similarity.tsv'
    chainData = dataDirectory + '/' + 'Chains.tsv'

    #===========================================================================
    # Extract the chains meeting the user's criteria.
    #===========================================================================
    if verboseOutput:
        print('Now extracting the chains to cull.')

    chainsOfInterest = {}  # The chains that meet the user specified quality criteria.
    readChainData = open(chainData, 'r')
    readChainData.readline()  # Strip the header.
    for i in readChainData:
        # Parse the data file containing all the chains in the PDB, and record only those that meet the quality criteria.
        chunks = (i.strip()).split('\t')
        chain = chunks[0]
        resolution = float(chunks[1])
        rValue = float(chunks[2])
        sequenceLength = int(chunks[3])
        xRayNotUsed = chunks[4] == 'yes'
        alphaCarbonOnly = chunks[5] == 'yes'
        representativeGroup = chunks[6]
        invalid = ((xRayNotUsed and skipNonXray) or
                   (resolution < minResolution) or
                   (resolution > maxResolution) or
                   (rValue > maxRValue) or
                   (alphaCarbonOnly and skipAlphaCarbon) or
                   (minLength != -1 and sequenceLength < minLength) or
                   (maxLength != -1 and sequenceLength > maxLength)
                   )
        if not invalid:
            chainsOfInterest[chain] = representativeGroup
    readChainData.close()

    allValidChains = {}  # The chains that will be submitted for culling.
    if cullSubset:
        readUserChains = open(fileUserInputChains, 'r')
        for line in readUserChains:
            chain = line.strip()
            if chain in chainsOfInterest:
                allValidChains[chain] = chainsOfInterest[chain]
        readUserChains.close()
    else:
        allValidChains = chainsOfInterest

    if verboseOutput:
        print('Number of chains meeting the specified criteria: {0}.'.format(len(allValidChains)))

    #===========================================================================
    # Determine the representative chains to submit for culling.
    #===========================================================================
    representativeGroupings = {}
    for i in allValidChains:
        representativeGroupings[allValidChains[i]] = i

    if verboseOutput:
        print('Number of representative chains: {0}.'.format(len(representativeGroupings)))

    #===========================================================================
    # Create the adjacency matrix.
    #===========================================================================
    if verboseOutput:
        print('Now creating the adjacency matrix')

    readSimilarities = open(similarityData, 'r')
    readSimilarities.readline()  # Strip the header.
    adjList = {}
    for line in readSimilarities:
        chunks = (line.strip()).split('\t')
        representativeGroupA = chunks[0]
        representativeGroupB = chunks[1]
        similarity = float(chunks[2])
        addSimilarity = ((representativeGroupA in representativeGroupings) and
                         (representativeGroupB in representativeGroupings) and
                         (similarity >= sequenceIdentity)
                        )
        if addSimilarity:
            # The sequences are of interest and too similar.
            chainA = representativeGroupings[representativeGroupA]
            chainB = representativeGroupings[representativeGroupB]
            if chainA in adjList:
                adjList[chainA].add(chainB)
            else:
                adjList[chainA] = set([chainB])
            if chainB in adjList:
                adjList[chainB].add(chainA)
            else:
                adjList[chainB] = set([chainA])
    readSimilarities.close()

    if verboseOutput:
        print('Number of similarity relationships: {0}.'.format(int(sum([len(adjList[i]) for i in adjList]) / 2)))

    #===========================================================================
    # Run the culling.
    #===========================================================================
    if verboseOutput:
        print('Now performing the culling.')
    proteinsToCull = Leafcull.main(adjList)
    proteinsToKeep = [i for i in representativeGroupings.values() if not i in proteinsToCull]

    if verboseOutput:
        print('{0} protins kept and {1} removed.'.format(len(proteinsToKeep), len(proteinsToCull)))

    #===========================================================================
    # Save the results.
    #===========================================================================
    if verboseOutput:
        print('Now saving the results.')

    # Write out the non-redundant chains.
    writeOutKeepList = open(outputLocation + '/KeptList.txt', 'w')
    writeOutKeepList.write('\n'.join(proteinsToKeep))
    writeOutKeepList.close()

    # Write out the redundant chains.
    writeOutKeepList = open(outputLocation + '/CulledList.txt', 'w')
    writeOutKeepList.write('\n'.join(proteinsToCull))
    writeOutKeepList.close()

    if verboseOutput:
        print('Results saved.')

if __name__ == '__main__':
    main(sys.argv[1:])