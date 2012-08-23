'''
Created on 5 Feb 2011

@author: Simon Bull
'''

import sparsematrix

def pdb_chain_main(similarities, cutoffPercent=20, representativeChains=set([])):
    """Create a sparse matrix representation of a protein similarity graph.
    
    Returns a sparsematrix representing the protein similarity graph, and a list of the names of the protein in the sparsematrix.
    The proteins in the list are ordered such that the name of the node in the graph identified as node i is at position i
    in the list of names.
    
    @param similarities: A record of the percentage sequence identity between the chains up for culling.
    @type similarities : string (file name)
    @param cutoffPercent: A percentage similarity > this parameter is deemed to be too similar.
    @type cutoffPercent :  float
    @param representativeChains: The names of the chains that will compose the protein similarity graph.
    @type representativeChains:  set
    return @type: SparseMatrix, list
    return @use:  The SparseMatrix representation of the protein similarity graph, a list of the names of the proteins in the graph (element i of this list is the name of the protein represented by node i in the graph)
    
    """

    proteinNames = set([])  # Store the names of all the proteins found to be too similar to another protein
    similarProteins = set([])  # Store the pairs that are too similar

    readSimilarityData = open(similarities, 'r')
    for i in readSimilarityData:
        chunks = (i.strip()).split('\t')
        chainA = chunks[1]
        chainB = chunks[3]
        similarity = float(chunks[5])
        matchLength = int(chunks[6])
        if chainA in representativeChains and chainB in representativeChains and similarity > cutoffPercent and chainA != chainB:
            proteinNames.add(chainA)
            proteinNames.add(chainB)
            similarProteins.add(tuple(sorted([chainA, chainB])))
    readSimilarityData.close()
        
    proteinNames = list(proteinNames)
    proteinNames.sort()
    similarProteins = list(similarProteins)
    indexDict = dict((proteinNames[x], x) for x in range(len(proteinNames)))
    
    # Create the sparse matrix
    adjacent = sparsematrix.SparseMatrix(len(proteinNames))
    xValues = [indexDict[x] for (x,y) in similarProteins]
    yValues = [indexDict[y] for (x,y) in similarProteins]
    adjacent.addlist(xValues, yValues)
    adjacent.addlist(yValues, xValues)
    
    return adjacent, proteinNames

def pdb_entry_main(similarities, cutoffPercent=20, representativeChains=set([]), representativesReverse=set([])):
    """Create a sparse matrix representation of a protein similarity graph.
    
    Returns a sparsematrix representing the protein similarity graph, and a list of the names of the protein in the sparsematrix.
    The proteins in the list are ordered such that the name of the node in the graph identified as node i is at position i
    in the list of names.
    
    @param similarities: A record of the percentage sequence identity between the chains up for culling.
    @type similarities : string (file name)
    @param cutoffPercent: A percentage similarity > this parameter is deemed to be too similar.
    @type cutoffPercent :  float
    @param representativeChains: The names of the chains that will compose the protein similarity graph.
    @type representativeChains:  set
    @param representativesReverse: A mapping of representative chains to the chains that they represent.
    @type representativesReverse:  dictionary
    return @type: SparseMatrix, list
    return @use:  The SparseMatrix representation of the protein similarity graph, a list of the names of the proteins in the graph (element i of this list is the name of the protein represented by node i in the graph)
    
    """

    entryNames = set([])  # Store the names of all the proteins found to be too similar to another protein
    similarEntries = set([])  # Store the pairs that are too similar

    # Create a representaive chain to entry mapping.
    repChainToEntry = {}
    for i in representativeChains:
        if representativesReverse.has_key(i):
            repChainToEntry[i] = (iter(representativesReverse[i]).next())[:4]  # Only use one entry, as all entries will have the same edgs in the similarity graph.
        else:
            repChainToEntry[i] = i[:4]  # Only use one entry, as all entries will have the same edgs in the similarity graph.

    readSimilarityData = open(similarities, 'r')
    for i in readSimilarityData:
        chunks = (i.strip()).split('\t')
        chainA = chunks[1]
        chainB = chunks[3]
        similarity = float(chunks[5])
        matchLength = int(chunks[6])
        if chainA in representativeChains and chainB in representativeChains and similarity > cutoffPercent and chainA != chainB:
            entryA = repChainToEntry[chainA]
            entryNames |= set([entryA])
            entryB = repChainToEntry[chainB]
            entryNames |= set([entryB])
            if entryA != entryB:
                similarEntries.add(tuple(sorted([entryA, entryB])))
    readSimilarityData.close()

    entryNames = list(entryNames)
    entryNames.sort()
    similarEntries = list(similarEntries)
    indexDict = dict((entryNames[x], x) for x in range(len(entryNames)))
    
    # Create the sparse matrix
    adjacent = sparsematrix.SparseMatrix(len(entryNames))
    xValues = [indexDict[x] for (x,y) in similarEntries]
    yValues = [indexDict[y] for (x,y) in similarEntries]
    adjacent.addlist(xValues, yValues)
    adjacent.addlist(yValues, xValues)
    
    return adjacent, entryNames

def intra_entry_main(similarities, cutoffPercent, representativeChains, entryToRepChain):
    """Create a sparse matrix representation of a protein similarity graph.
    
    Returns a sparsematrix representing the protein similarity graph, and a list of the names of the protein in the sparsematrix.
    The proteins in the list are ordered such that the name of the node in the graph identified as node i is at position i
    in the list of names.
    
    @param similarities: A record of the percentage sequence identity between the chains up for culling.
    @type similarities : string (file name)
    @param cutoffPercent: A percentage similarity > this parameter is deemed to be too similar.
    @type cutoffPercent :  float
    @param representativeChains: The names of the chains that will compose the protein similarity graph.
    @type representativeChains:  set
    @param entryToRepChain: A mapping of entries to their representative chains.
    @type entryToRepChain:  dictionary
    return @type: SparseMatrix, list
    return @use:  The SparseMatrix representation of the protein similarity graph, a list of the names of the proteins in the graph (element i of this list is the name of the protein represented by node i in the graph)
    
    """

    similarProteins = {}  # Store the pairs that are too similar
    
    readSimilarityData = open(similarities, 'r')
    for i in readSimilarityData:
        chunks = (i.strip()).split('\t')
        chainA = chunks[1]
        chainB = chunks[3]
        similarity = float(chunks[5])
        matchLength = int(chunks[6])
        if chainA in representativeChains and chainB in representativeChains and similarity > cutoffPercent and chainA != chainB:
            similarProteins[tuple(sorted([chainA, chainB]))] = None
    readSimilarityData.close()

    adjList = []
    namesList = []
    for entry in entryToRepChain:
        sortedRepChains = sorted(entryToRepChain[entry])
        chainNames = set([])
        adjacencies = set([])
        for i in range(len(sortedRepChains)):
            for j in sortedRepChains[1+1:]:
                key = tuple([sortedRepChains[i], j])
                if similarProteins.has_key(key):
                    chainNames.add(sortedRepChains[i])
                    chainNames.add(j)
                    adjacencies.add(key)

        if chainNames == set([]):
            # If there are no similarities between the representative chains.
            continue

        chainNames = list(chainNames)
        chainNames.sort()
        adjacencies = list(adjacencies)
        indexDict = dict((chainNames[x], x) for x in range(len(chainNames)))
        
        # Create the sparse matrix
        adjacent = sparsematrix.SparseMatrix(len(chainNames))
        xValues = [indexDict[x] for (x,y) in adjacencies]
        yValues = [indexDict[y] for (x,y) in adjacencies]
        adjacent.addlist(xValues, yValues)
        adjacent.addlist(yValues, xValues)

        adjList.append(adjacent)
        namesList.append(chainNames)

    return adjList, namesList

def user_seq_main(similarities, cutoffPercent=20):
    """Create a sparse matrix representation of a protein similarity graph.
    
    Returns a sparsematrix representing the protein similarity graph, and a list of the names of the protein in the sparsematrix.
    The proteins in the list are ordered such that the name of the node in the graph identified as node i is at position i
    in the list of names.
    
    @param similarities: A dictionary containing the information about the simialrity between the proteins.
    @type similarities : dictionary
    @param cutoffPercent: A percentage similarity > this parameter is deemed to be too similar.
    @type cutoffPercent :  float
    return @type: SparseMatrix, list
    return @use:  The SparseMatrix representation of the protein similarity graph, a list of the names of the proteins in the graph (element i of this list is the name of the protein represented by node i in the graph)
    
    """

    proteinNames = set([])  # Store the names of all the proteins found to be too similar to another protein
    similarProteins = set([])  # Store the pairs that are too similar

    for i in similarities.keys():
        if i[0] == i[1] or similarities[i]['Identity'] <= cutoffPercent:
            continue
        proteinNames.add(i[0])
        proteinNames.add(i[1])
        similarProteins.add(i)
        
    proteinNames = list(proteinNames)
    proteinNames.sort()
    similarProteins = list(similarProteins)
    indexDict = dict((proteinNames[x], x) for x in range(len(proteinNames)))
    
    # Create the sparse matrix
    adjacent = sparsematrix.SparseMatrix(len(proteinNames))
    xValues = [indexDict[x] for (x,y) in similarProteins]
    yValues = [indexDict[y] for (x,y) in similarProteins]
    adjacent.addlist(xValues, yValues)
    adjacent.addlist(yValues, xValues)
    
    return adjacent, proteinNames

        
