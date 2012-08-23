'''
Created on 28 Mar 2011

@author: Simon Bull
'''

def pruneGraph(adjList, IDs):
    """The method by which the Leaf algorithm determines which nodes to remove from the graph.
    
    @param adjList: An adjacency list representation of the protein similarity graph.
    @type adjList : dictionary
    @param IDs: A list of the numerical indices of the nodes in the graph.
    @type IDs : list
    
    """
    
    removeList = []
    nodesInGraph = adjList.keys()

    # Determine number of neighbours for each node.
    neighbours = {}
    for i in adjList.keys():
        numNeighbours = len(adjList[i])
        if neighbours.has_key(numNeighbours):
            neighbours[numNeighbours] |= set([i])
        else:
            neighbours[numNeighbours] = set([i])
    # Fill in the blank keys.
    neighbourKeys = neighbours.keys()
    for i in set(range(max(neighbourKeys))) - set(neighbourKeys):
        if neighbours.has_key(i):
            continue
        else:
            neighbours[i] = set([])
        
    
    while True:

        # Determine the maximum number of neighbours.
        maxNeighbours = max(neighbours.keys())
        while neighbours[maxNeighbours] == set([]):
            del neighbours[maxNeighbours]
            maxNeighbours = max(neighbours.keys())
        
        # If there are no nodes with neighbours then exit.
        if maxNeighbours == 0:
            return removeList

        nClique = 1
        while nClique <= maxNeighbours:
            # Get the nodes with nClique neighbours.
            if neighbours.has_key(nClique):
                nodesOfInterest = neighbours[nClique]
                # For every node of interest see if the neighbours of the node are all neighbours of each other (i.e. a clique).
                for i in nodesOfInterest:
                    neighboursOfInterest = set(adjList[i])
                    while len(neighboursOfInterest) > 1:
                        toCheck = neighboursOfInterest.pop()
                        if neighboursOfInterest <= set(adjList[toCheck]):
                            continue
                        else:
                            break
                    else:
                        toRemove = [j for j in adjList[i]]
                        neighbours[nClique] -= set([i])
                        neighbours[0] |= set([i])
                        for j in toRemove:
                            removeList.append(j)
                            # Update the list of neighbours for each node that toRemove is adjacent to.
                            for k in adjList[j]:
                                numNeighbours = len(adjList[k])
                                neighbours[numNeighbours] -= set([k])
                                neighbours[numNeighbours - 1] |= set([k])
                                adjList[k].remove(j)
                            # Update the adjacency list to reflect the removal of to remove.
                            numNeighbours = len(adjList[j])
                            neighbours[numNeighbours] -= set([j])
                            neighbours[0] |= set([j])
                            adjList[j] = []
                        nClique = 1
                        break
                else:
                    # No clique found.
                    nClique += 1
            else:
                # No nodes with the desired number of neighbours.
                nClique += 1

        ########################################
        # Perform the NeighbourCull operation. #
        ########################################
        # Re-calculate this, as it may have changed since it was last calculated.
        maxNeighbours = max(neighbours.keys())
        while neighbours[maxNeighbours] == set([]):
            del neighbours[maxNeighbours]
            maxNeighbours = max(neighbours.keys())

        # If there are no nodes with neighbours then exit.
        if maxNeighbours == 0:
            return removeList

        # Get the IDs of the nodes with the max number of neighbours.
        nodesWithMaxNeighbours = list(neighbours[maxNeighbours])
        # If there is more than one node with the maximum number of neighbours determine which node to remove.
        if len(nodesWithMaxNeighbours) != 1:
            # Determine the number of neighbours for each node.
            extendedNeighbourhood = [adjList[x] + [x] for x in nodesWithMaxNeighbours]
            extendedNeighbourhood = [set().union(*[adjList[i] for i in a]) for a in extendedNeighbourhood]
            # Determine the size of each extended neighbourhood, and which nodes have the min size.
            sizes = [len(x) for x in extendedNeighbourhood]
            minSize = min(sizes)
            toRemove = nodesWithMaxNeighbours[sizes.index(minSize)]
        else:
            toRemove = nodesWithMaxNeighbours[0]
      
        removeList.append(toRemove)
        # Update the list of neighbours for each node that toRemove is adjacent to.
        for i in adjList[toRemove]:
            numNeighbours = len(adjList[i])
            neighbours[numNeighbours] -= set([i])
            neighbours[numNeighbours - 1] |= set([i])
            adjList[i].remove(toRemove)
        # Update the adjacency list to reflect the removal of to remove.
        adjList[toRemove] = []
        neighbours[maxNeighbours] -= set([toRemove])
        neighbours[0] |= set([toRemove])

def main(adj, names):
    """Use the Leaf heuristic method to calculate an approximation to the maximum independent set.
    
    Returns a list of the proteins to keep and a list of the proteins to cull. The list of proteins to keep only contains the
    names of the proteins in the protein similarity graph that should be kept. If there are any proteins that were not
    included in adj (for example proteins with no neighbours), then these will NOT be included in the list of proteins to keep.
    See the README for a more in depth description of this.
    
    @param adj: A sparsematrix representation of the protein similarity graph
    @type adj : sparsematrix
    @param names: A list of the names of the proteins in adj. Ordered such that the name of the protein represented by node i
                  in adj is located at names[i].
    @type names : list
    return @type: list, list
    return @use:  the proteins from the graph to be removed, the proteins from the graph to be kept (NOT THE SAME AS THE MIS APPROXIMATION)

    """

    rem = pruneGraph(adj.adjList(), range(len(names)))
    proteinsToCull = [names[x] for x in rem]
    proteinsToKeep = [names[x] for x in range(len(names)) if x not in rem]

    return proteinsToCull, proteinsToKeep
