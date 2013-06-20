'''
Created on 5 Feb 2011

@author: Simon Bull
'''

import os
import shutil
import sys
import getopt
import argparse
import time

import adjlistcreation
import Leafcull
import checkPDBinput


def main(args):
    """Runs the protein culling.

    @param args: The command line arguments.
    @type args: A Python list.

    """

    #===========================================================================
    # Parse the user's input.
    #===========================================================================
    parser = argparse.ArgumentParser(description=('Generate a non-redundant dataset of chains or entries from the PDB. ' +
                                                  'Please see the README for more information on how to use this program.'),
                                     epilog=('This program is designed to cull a dataset of protein sequences so that no ' +
                                             'two sequences have a sequence identity greater than the specified threshold ' +
                                             'percentage. The method used is the Leaf heuristic, which is described in a paper located at ' +
											 'http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0055484. ' +
                                             'A server to perform the culling can be found at http://www.bioinf.manchester.ac.uk/leaf/.')
                                     )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-i', '--inputFile', help='The location of the input file. (Required type: %(type)s).',
                       metavar="inputFile", type=str, default='')
    group.add_argument('-w', '--whole', help='Cull chains/entries from the entire PDB. (Default value: %(default)s).',
                       action='store_true', default=False)
    group.add_argument('-n', '--organism', help="Cull chains/entries from a specified organism. Replace spaces with underscores ('_'). (Required type: %(type)s).",
                       metavar="organismName", type=str, default='')

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
    parser.add_argument('-c', '--cullbyentry', help='Whether the culling should be by entry. (Default value: Do not cull by entry).',
                        action='store_true', default=False, required=False)
    parser.add_argument('-e', '--intraentry', help='Whether culling should be performed within entries. Only useable when culling by entry is being used. (Default value: No intra entry culling).',
                        action='store_true', default=False, required=False)
    parser.add_argument('-s', '--intrapercent', help='The maximum percent sequence identity between chains in the same entry 5 <= maxIntraEntryPercent < 100 must be true. Only usable when culling by entry with intra entry culling is being used. (Required type: %(type)s, default value: %(default)s).',
                        metavar="maxIntraEntryPercent", type=float, default=20, required=False)
    parser.add_argument('-d', '--dataloc', help='The location of the directory that contains the processed PDB data (required type: %(type)s, default value: in the directory this script is being called from).',
                        metavar="dataFileDir", type=str, default='', required=False)
    parser.add_argument('-o', '--output', help='The name of the output directory to create in the current working directory. (Required type: %(type)s, default value: a directory called %(default)s in the current working directory).',
                        metavar="outputFolder", type=str, default='PDBCullResults', required=False)
    parser.add_argument('-v', '--verbose', help='Whether status updates should be displayed. (Default value: No status updates).',
                        action='store_true', default=False, required=False)
    args = parser.parse_args()

    wholePDB = args.whole
    cullByOrganism = ' '.join(args.organism.split('_'))
    userInput = args.inputFile
    sequenceIdentity = args.percent
    minResolution = args.minRes
    maxResolution = args.maxRes
    maxRValue = args.rval
    minLength = args.minLen
    maxLength = args.maxLen
    skipNonXray = not args.usenonxray
    skipAlphaCarbon = not args.usealpha
    cullByChain = not args.cullbyentry
    performIntraEntryCulling = args.intraentry
    intraEntrySequenceIdentity = args.intrapercent
    dataFileDir = args.dataloc
    cullOperationID = args.output
    verboseOutput = args.verbose

    #===========================================================================
    # Validate the user's input.
    #===========================================================================
    toExit = False
    if not wholePDB:
        if userInput != '' and not os.path.isfile(userInput):
            print 'The input location supplied is not a valid file location.'
            toExit = True

    if dataFileDir == '':
        # Load the default location that contains the data directory.
        srcLocation = os.path.abspath(__file__)
        srcLocation = '\\'.join(srcLocation.split('\\')[:-1])
        dataFileDir = srcLocation + '\\PDBData'
    elif not os.path.isdir(dataFileDir):
        # Check that the user submitted a directory as the data directory.
        print 'The location supplied for the parsed PDB data is not a valid directory location.'
        toExit = True

    if sequenceIdentity < 5 or sequenceIdentity >= 100:
        print 'The maximum allowable percentage sequence similarity must be no less than 5, and less than 100.'
        toExit = True

    if minResolution < 0 or minResolution > 100:
        print 'The valid range for the minimum resolution is 0 - 100.'
        toExit = True

    if maxResolution < 0 or maxResolution > 100:
        print 'The valid range for the maximum resolution is 0 - 100.'
        toExit = True

    if minResolution > maxResolution:
        print 'The minimum resolution must be less than or equal to the maximum resolution.'
        toExit = True

    if maxRValue < 0 or maxRValue > 1:
        print 'The valid range for the maximum R value is 0 - 1.'
        toExit = True

    if not cullByChain and performIntraEntryCulling:
        intraEntrySequenceIdentity = float(intraEntrySequenceIdentity)
        if intraEntrySequenceIdentity < 5 or intraEntrySequenceIdentity >= 100:
            print 'The maximum allowable intr entry percentage sequence similarity must be no less than 5, and less than 100.'
            toExit = True
    elif cullByChain and performIntraEntryCulling:
        print 'WARNING: Culling by entry is not enabled, but intra entry culling is selected. If you choose to continue the culling will NOT be done by entry. Continue (y for yes, anything else for no)?'
        cont = raw_input('--> ')
        if cont.upper() != 'Y':
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

    # Create the directory to store the output in.
    cwd = os.getcwd()
    if cullOperationID == 'PDBCullResults':
        outputLocation = cwd + '\\' + cullOperationID
    else:
        outputLocation = cullOperationID
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

    #===========================================================================
    # Extract the PDB data.
    #===========================================================================
    PDBEntriesData = dataFileDir + '\\' + 'AllPDBEntries.txt'
    chainTypeData = dataFileDir + '\\' + 'ChainType.txt'
    representativeData = dataFileDir + '\\' + 'Representative.txt'
    similarityData = dataFileDir + '\\' + 'Similarity.txt'
    proteinData = dataFileDir + '\\' + 'ProteinInformation.txt'

    PDBEntriesList = []
    readPDBEntriesData = open(PDBEntriesData, 'r')
    for i in readPDBEntriesData:
        PDBEntriesList.append(i.strip())
    readPDBEntriesData.close()

    chainTypeDict = {}
    readChainTypeData = open(chainTypeData, 'r')
    for i in readChainTypeData:
        chunks = (i.strip()).split('\t')
        chainTypeDict[chunks[0]] = chunks[1]
    readChainTypeData.close()

    proteinDict = {}
    readProteinData = open(proteinData, 'r')
    for i in readProteinData:
        chunks = (i.strip()).split('\t')
        chain = chunks[0]
        entry = chunks[1]
        organism = chunks[10]
        proteinDict[chain] = {'entry' : entry, 'organism' : organism}
    readProteinData.close()

    #===========================================================================
    # Process the PDB data.
    #===========================================================================
    if verboseOutput:
        print 'Now parsing the processed PDB data.'
    if userInput != '':
        # Check that the user's input contains no invalid chains/entries. Only necessary if the whole PDB is not being culled.
        userInputFile = open(userInput, 'r')
        userInput = userInputFile.read()
        userInputFile.close()
        userInput = userInput.split('\n')
        inputList = [i.strip() for i in userInput]
        if cullByChain:
            # If the user has selected to cull 'by chain'.
            processedInputList = set([])
            for i in inputList:
                # Check the user input to see if any of the 'chains' supplied by the user might actually be entries.
                # This is fine, but it needs to be checked for. Any entries supplied will be converted into the chains
                # that correspond to the entry supplied (i.e. entry E will be converted to chains EA, EB, EC, etc.).
                if len(i) == 4:
                    # Possibly an entry.
                    chainsFromEntry = [j for j in proteinDict.keys() if proteinDict[j]['entry'] == i]
                    if chainsFromEntry == []:
                        # If there are not chains corresponding to the 4 character input, then it is likely that the 4
                        # character input is invalid. In this case simply keep the 4 character input, and wait for the
                        # input checking script to raise an error.
                        processedInputList.add(i)
                    else:
                        # Chains were found, and therefore the 4 character input was in fact an entry.
                        for j in chainsFromEntry:
                            processedInputList.add(j)
                else:
                    # If i could not be an entry then just add it to the list of chains to be checked.
                    processedInputList.add(i)
            processedInputList = list(processedInputList)
            retCode, retVal = checkPDBinput.main(processedInputList, allChains=chainTypeDict, checkType='chain')
        else:
            # If the user has selected to cull 'by entry'.
            allProtEntries = set([chainTypeDict[i]['chain'][:4] for i in chainTypeDict.keys() if chainTypeDict[i]['chainType'] == 'Protein'])
            retCode, retVal = checkPDBinput.main(inputList, allEntries=PDBEntriesList, allProtEntries=allProtEntries,
                                                 checkType='entry')
        if retCode != 0:
            # An error was found in the input.
            print retVal
            sys.exit()
    elif cullByOrganism != '':
        if cullByChain:
            # Collect all chains belonging to the organism specified.
            retVal = [i for i in proteinDict.keys() if proteinDict[i]['organism']]
        else:
            # Collect all entries belonging to the organism specified.
            retVal = list(set([proteinDict[i]['entry'] for i in proteinDict.keys() if proteinDict[i]['organism']]))
        if len(retVal) < 2:
            # Not enough chains/entries of the given organism type.
            print 'There are less than 2 chains/entries in the PDB from the organism you entered. This is possibly due ',
            print 'to a misspelling of the organism name.'
            sys.exit()
        retVal = '\n'.join(retVal)

    if verboseOutput:
        startTime = time.time()
    if not cullByChain:
        # If the method of culling is 'by entry', record the entries and convert the entries to their corresponding chains.
        if not wholePDB:
            userInput = retVal.split('\n')
            entriesUsed = userInput  # All redundant and non-redundant entries.
            potentialChains = []  # All redudant and non-redundant chains.
            chainsToCull = set([])  # The chains that will be used in the culling.
            readProteinData = open(proteinData, 'r')
            for i in readProteinData:
                # Parse the data file containing all the chains in the PDB, and record only those which are members of an entry in the user's input.
                chunks = (i.strip()).split('\t')
                chain = chunks[0]
                entry = chunks[1]
                if entry in userInput:
                    potentialChains.append(chain)
                    experimentType = chunks[2]
                    resolution = float(chunks[3])
                    rValueObs = float(chunks[4])
                    alphaCarbonOnly = False if chunks[6] == '0' else True
                    sequence = chunks[11]
                    invalid = ((experimentType != 'XRAY' and skipNonXray) or
                               (resolution < minResolution) or
                               (resolution > maxResolution) or
                               (rValueObs > maxRValue) or
                               (alphaCarbonOnly and skipAlphaCarbon) or
                               (minLength != -1 and len(sequence) < minLength) or
                               (maxLength != -1 and len(sequence) > maxLength)
                               )
                    if not invalid:
                        chainsToCull.add(chain)
        else:
            entriesUsed = set([])  # All redundant and non-redundant entries.
            potentialChains = []  # All redudant and non-redundant chains.
            chainsToCull = set([])  # The chains that will be used in the culling.
            readProteinData = open(proteinData, 'r')
            for i in readProteinData:
                # Parse the data file containing all the chains in the PDB, and record only those which are members of an entry in the user's input.
                chunks = (i.strip()).split('\t')
                chain = chunks[0]
                entry = chunks[1]
                potentialChains.append(chain)
                entriesUsed.add(entry)
                experimentType = chunks[2]
                resolution = float(chunks[3])
                rValueObs = float(chunks[4])
                alphaCarbonOnly = False if chunks[6] == '0' else True
                sequence = chunks[11]
                invalid = ((experimentType != 'XRAY' and skipNonXray) or
                           (resolution < minResolution) or
                           (resolution > maxResolution) or
                           (rValueObs > maxRValue) or
                           (alphaCarbonOnly and skipAlphaCarbon) or
                           (minLength != -1 and len(sequence) < minLength) or
                           (maxLength != -1 and len(sequence) > maxLength)
                           )
                if not invalid:
                    chainsToCull.add(chain)
            entriesUsed = list(entriesUsed)
        entriesToCull = set([i[:4] for i in chainsToCull])
    else:
        # If the method of culling is 'by chain', record the chains input by the user.
        if not wholePDB:
            userInput = retVal.split('\n')
            potentialChains = []  # All redudant and non-redundant chains.
            chainsToCull = set([])  # The chains that will be used in the culling.
            readProteinData = open(proteinData, 'r')
            for i in readProteinData:
                # Parse the data file containing all the chains in the PDB, and record only those which are members in the user's input.
                chunks = (i.strip()).split('\t')
                chain = chunks[0]
                if chain in userInput:
                    potentialChains.append(chain)
                    experimentType = chunks[2]
                    resolution = float(chunks[3])
                    rValueObs = float(chunks[4])
                    alphaCarbonOnly = False if chunks[6] == '0' else True
                    sequence = chunks[11]
                    invalid = ((experimentType != 'XRAY' and skipNonXray) or
                               (resolution < minResolution) or
                               (resolution > maxResolution) or
                               (rValueObs > maxRValue) or
                               (alphaCarbonOnly and skipAlphaCarbon) or
                               (minLength != -1 and len(sequence) < minLength) or
                               (maxLength != -1 and len(sequence) > maxLength)
                               )
                    if not invalid:
                        chainsToCull.add(chain)
            readProteinData.close()
        else:
            potentialChains = []  # All redudant and non-redundant chains.
            chainsToCull = set([])  # The chains that will be used in the culling.
            readProteinData = open(proteinData, 'r')
            for i in readProteinData:
                # Parse the data file containing all the chains in the PDB, and record only those which are members in the user's input.
                chunks = (i.strip()).split('\t')
                chain = chunks[0]
                potentialChains.append(chain)
                experimentType = chunks[2]
                resolution = float(chunks[3])
                rValueObs = float(chunks[4])
                alphaCarbonOnly = False if chunks[6] == '0' else True
                sequence = chunks[11]
                invalid = ((experimentType != 'XRAY' and skipNonXray) or
                           (resolution < minResolution) or
                           (resolution > maxResolution) or
                           (rValueObs > maxRValue) or
                           (alphaCarbonOnly and skipAlphaCarbon) or
                           (minLength != -1 and len(sequence) < minLength) or
                           (maxLength != -1 and len(sequence) > maxLength)
                           )
                if not invalid:
                    chainsToCull.add(chain)
            readProteinData.close()

    if verboseOutput:
        print 'Potential chains: ', len(potentialChains), 'Valid chains: ', len(chainsToCull)

    # Determine representative chain information.
    # representatives records the non-representative to representative chain mapping for the non-representative
    # chains in the set of chains to cull.
    representatives = {}
    readRepresentativeData = open(representativeData, 'r')
    for i in readRepresentativeData:
        chunks = (i.strip()).split('\t')
        nonreprChain = chunks[0]
        reprChain = chunks[1]
        if nonreprChain in chainsToCull:
            representatives[nonreprChain] = reprChain
    readRepresentativeData.close()
    # representativeChains records the set of representative chains that cover all the chains in the set of chains to cull.
    # This means that if a chain in chainsToCull is a representative itself then it is in representativeChains, and if
    # a chain in chainsToCull is not a representative, then its representative chain is in representativeChains.
    representativeChains = set([i if not representatives.has_key(i) else representatives[i] for i in chainsToCull])
    # representativesReverse records for each chain that represents at least one chain in chainsToCull, a set of the
    # non-representative chains in chainsToCull that it represents.
    # For example, if chainsToCull == [a, b, c, d], and a and b are non-representative chains represented by chain q,
    # then representativesReverse[q] = set([a, b]).
    representativesReverse = {}
    for i in representatives.keys():
        reprChain = representatives[i]
        if representativesReverse.has_key(reprChain):
            representativesReverse[reprChain].add(i)
        else:
            representativesReverse[reprChain] = set([i])
    representativesReverseKeys = representativesReverse.keys()

    if verboseOutput:
        print 'Now beginning the culling. Time elapsed: ', time.time() - startTime

    if not cullByChain:
        # Perform the 'by entry' culling.
        # Determine the redundant user input entries.
        removedInput = cull_main(similarityData, sequenceIdentity, representativeChains, 'entry', representativesReverse, verboseOutput, startTime)
        # Determine the non-redundant user input entries.
        keptInput = set([i[:4] for i in entriesToCull if i[:4] not in removedInput])
        if performIntraEntryCulling and intraEntrySequenceIdentity < 100:
            # Perform intra-entry culling.
            if verboseOutput:
                print 'Now performing intra-entry culling. Time elapsed: ', time.time() - startTime
            entryToChain = {}  # Records all the chains within each non-redudnant user input entry.
            chainsOfInterest = set([])  # Records all the chains that are in a non-redundant user input entry.
            readProteinData = open(proteinData, 'r')
            for i in readProteinData:
                # Determine all the chains for each non-redundant user input entry.
                chunks = (i.strip()).split('\t')
                chain = chunks[0]
                entry = chunks[1]
                sequence = chunks[11]
                invalid = (minLength != -1 and len(sequence) < minLength) or (maxLength != -1 and len(sequence) > maxLength)
                if entry in keptInput and not invalid:
                    chainsOfInterest.add(chain)
                    if entryToChain.has_key(entry):
                        entryToChain[entry].append(chain)
                    else:
                        entryToChain[entry] = [chain]
            readProteinData.close()

            representatives = {}  # Records the non-representative to representative chain mapping for all chains in non-redundant user input entries.
            readRepresentativeData = open(representativeData, 'r')
            for i in readRepresentativeData:
                # Determine the representative for each chain that is in a non-redundant user input entry (for all chains in chainsOfInterest).
                chunks = (i.strip()).split('\t')
                nonreprChain = chunks[0]
                reprChain = chunks[1]
                if nonreprChain in chainsOfInterest:
                    representatives[nonreprChain] = reprChain
            readRepresentativeData.close()
            representativeChains = set([i if not representatives.has_key(i) else representatives[i] for i in chainsOfInterest])  # Records all representative chains for the set of chains that are in the non-redundant user input entries.
            representativesReverse = {}  # Records all the non-representative chains represented by the representative chains in representatives.
            for i in representatives.keys():
                reprChain = representatives[i]
                if representativesReverse.has_key(reprChain):
                    representativesReverse[reprChain].add(i)
                else:
                    representativesReverse[reprChain] = set([i])

            entryToRepChain = dict([(i, set([])) for i in entryToChain.keys()])  # Maps entries to their representative chains.
            for i in chainsOfInterest:
                entry = i[:4]
                if representatives.has_key(i):
                    entryToRepChain[entry].add(representatives[i])
                else:
                    entryToRepChain[entry].add(i)

            keptInputChains = set([])

            for i in keptInput:
                if len(entryToRepChain[i]) == 1:
                    # If the entry's chains are all representated by one chain, then all the chains are identical. A random chain from the entry should be kept.
                    keptInputChains.add(entryToChain[i][0])
                    del entryToRepChain[i]

            adjList, namesList = adjlistcreation.intra_entry_main(similarityData, intraEntrySequenceIdentity, representativeChains, entryToRepChain)

            for i in range(len(adjList)):
                # Perform the intra-entry culling for each entry that needs it.
                chainsToCull = Leafcull.main(adjList[i], namesList[i])
                keptReprChains = [j for j in namesList[i] if not j in chainsToCull]
                for j in keptReprChains:
                    # Calculate the kept input chains.
                    if representativesReverse.has_key(i):
                        # If the representative chain that was kept has non-representative chains in the input, then select one of them.
                        keptInputChains.add(iter(representativesReverse[j]).next())
                    else:
                        # If the representative chain that was kept has no non-representative chains in the input, then the representative chain was in the input. Keep it.
                        keptInputChains.add(j)
        else:
            # If intra-entry cuilling is not being used.
            keptInputChains = set([i for i in chainsToCull if i[:4] in keptInput])
    else:
        # Perform the 'by chain' culling.
        # Calculate the redundant representative chains.
        removedReprChains = set(cull_main(similarityData, sequenceIdentity, representativeChains, 'chain', {}, verboseOutput, startTime))
        # Determine the non-redundant representative chains.
        keptReprChains = [i for i in representativeChains if i not in removedReprChains]
        keptInputChains = set([])
        # Determine the non-redundant user input chains.
        for i in keptReprChains:
            # Calculate the kept input chains.
            if representativesReverse.has_key(i):
                # If the representative chain that was kept has non-representative chains in the input, then select one of them.
                keptInputChains.add(iter(representativesReverse[i]).next())
            else:
                # If the representative chain that was kept has no non-representative chains in the input, then the representative chain was in the input. Keep it.
                keptInputChains.add(i)
        # Determine the redundant user input chains.
        removedInput = sorted([i for i in potentialChains if i not in keptInputChains])

    if verboseOutput:
        print 'Now saving results. Time elapsed: ', time.time() - startTime

    # Write out the redundant chains/entries.
    writeOutRem = open(outputLocation + '\\Removed.txt', 'w')
    writeOutRem.write('\n'.join(removedInput))
    writeOutRem.close()

    # Determine the non-redundant chain/entry data to be output.
    keptInputOutput = 'IDs\tlength\tExptl.\tresolution\tR-factor\tFreeRvalue\n'
    fastaOutput = ''
    entryStats = {}
    readProteinData = open(proteinData, 'r')
    for i in readProteinData:
        chunks = (i.strip()).split('\t')
        chain = chunks[0]
        entry = chunks[1]
        experimentType = chunks[2]
        resolution = chunks[3]
        rValueObs = chunks[4]
        rValueFree = chunks[5]
        alphaCarbon = 'no' if chunks[6] == '0' else 'yes'
        description = chunks[7]
        dbName = chunks[8]
        dbCode = chunks[9]
        organism = chunks[10]
        sequence = chunks[11]
        if not cullByChain:
            if entry in keptInput:
                entryStats[entry] = {'len' : str(len(sequence)), 'expt' : experimentType, 'res' : resolution, 'rval' : rValueObs, 'freeRval' : rValueFree}
            if chain in keptInputChains:
                fastaOutput += '>' + '\t'.join([chain, str(len(sequence)), experimentType, resolution, rValueObs, rValueFree, alphaCarbon, description, '<' + dbName + ' ' + dbCode + '>', '[' + organism + ']']) + '\n' + sequence + '\n'
        else:
            if chain in keptInputChains:
                keptInputOutput += '\t'.join([chain, str(len(sequence)), experimentType, resolution, rValueObs, rValueFree]) + '\n'
                fastaOutput += '>' + '\t'.join([chain, str(len(sequence)), experimentType, resolution, rValueObs, rValueFree, alphaCarbon, description, '<' + dbName + ' ' + dbCode + '>', '[' + organism + ']']) + '\n' + sequence + '\n'
    readProteinData.close()

    if not cullByChain:
        keptInputOutput += '\n'.join(['\t'.join([i, entryStats[i]['len'], entryStats[i]['expt'], entryStats[i]['res'], entryStats[i]['rval'], entryStats[i]['freeRval']]) for i in sorted(entryStats.keys())])

    # Write out the non-redundant chain/entry statistics.
    writeOutKeepList = open(outputLocation + '\\KeptList.txt', 'w')
    writeOutKeepList.write(keptInputOutput)
    writeOutKeepList.close()

    # Write out the non-redundant FASTA file of chains.
    writeOutKeepFasta = open(outputLocation + '\\KeptFasta.fasta', 'w')
    writeOutKeepFasta.write(fastaOutput)
    writeOutKeepFasta.close()

    if verboseOutput:
        print 'Results saved. Total time elapsed: ', time.time() - startTime

def cull_main(similarities, thresholdPercentage, representativeChains, adjType, representativesReverse={}, verboseOutput=False, startTime=0):
    """Perform the PDB redundancy removal.

    @param similarities: A record of the percentage sequence identity between the chains/entries up for culling.
    @type similarities : file name (string)
    @param thresholdPercentage: The maximum permissible percentage sequence identity that any two chains/entries may possess.
    @type thresholdPercentage : float
    @param representativeChains: The names of the chains/entries that will compose the protein similarity graph.
    @type representativeChains:  set
    @param adjType: 'chain' or 'entry' indicating the type of culling to be performed.
    @type adjType:  string
    @param representativesReverse: A mapping of representative chains to the chains that they represent.
    @type representativesReverse:  dictionary
    @param verboseOutput: Whether status updates should be printed out to the user.
    @type verboseOutput:  boolean
    @param startTime: The time when the culling began. Used to output elapsed time.
    @type startTime:  integer
    return @type: list
    return @use:  The redundant proteins to be removed from teh dataset.

    """

    # Create the sparsematrix of the protein similarity graph.
    if verboseOutput:
        print 'Creating the adjacency matrix. Time elapsed: ', time.time() - startTime
    if adjType == 'chain':
        adjacent, proteinNames = adjlistcreation.pdb_chain_main(similarities, thresholdPercentage, representativeChains)
    elif adjType == 'entry':
        adjacent, proteinNames = adjlistcreation.pdb_entry_main(similarities, thresholdPercentage, representativeChains, representativesReverse)

    # Choose which proteins to remove from the similarity graph.
    if proteinNames == []:
        if verboseOutput:
            print 'No similarities found. Culling not needed. Time elapsed: ', time.time() - startTime
        # This is True if there are no similarities greater than the given percentage sequence identity. If there are no
        # chains that are too similar, then there is no need to cull any chains from the network.
        proteinsToCull = []
    else:
        if verboseOutput:
            print 'Performing the culling. Time elapsed: ', time.time() - startTime
        # Choose which chains to remove from the similarity graph.
        proteinsToCull, proteinsToKeep = Leafcull.main(adjacent, proteinNames)

    if verboseOutput:
        print 'Culling finished. Time elapsed: ', time.time() - startTime

    return proteinsToCull

if __name__ == '__main__':
    main(sys.argv[1:])