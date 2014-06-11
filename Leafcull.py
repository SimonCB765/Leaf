'''
Created on 28 Mar 2011

@author: Simon Bull
'''

import collections

def main(adjList):
    """The method by which the Leaf algorithm determines which nodes to remove from the dataset.

    :param adjList: An adjacency list representation of the protein similarity graph.
    :type adjList:  dictionary
    :returns :      The proteins that should be removed from the dataset.
    :type :         list

    """

    removeList = []
    if not adjList:
        # If the graph supplied is empty (i.e. no redundancy is present)
        return removeList

    # Determine number of neighbours for each node.
    neighbours = collections.defaultdict(set)
    for i in adjList:
        numNeighbours = len(adjList[i])
        neighbours[numNeighbours].add(i)
    # Fill in the blank keys.
    for i in set(range(max(neighbours))) - set(neighbours):
        neighbours[i]  # Defaultdict automatically creates an empty set here.

    while True:
        # Determine the maximum number of neighbours.
        maxNeighbours = max(neighbours)
        while maxNeighbours > 0 and (not neighbours[maxNeighbours]):
            del neighbours[maxNeighbours]
            maxNeighbours -= 1

        # If there are no nodes with neighbours then exit.
        if maxNeighbours == 0:
            return removeList

        nClique = 1
        while nClique <= maxNeighbours:
            nodesOfInterest = neighbours[nClique]  # Get the nodes with nClique neighbours.
            # For every node of interest see if the neighbours of the node are all neighbours of each other (i.e. a clique).
            for i in nodesOfInterest:
                neighboursOfInterest = adjList[i]
                if neighboursOfInterest.intersection(*[adjList[j].union([j]) for j in neighboursOfInterest]) == neighboursOfInterest:
                    # i's neighbours are all connected to one another, and therefore i participates in a clique with all of its neighbours where
                    # it is connected only to nodes in the clique.
                    toRemove = set(adjList[i])  # Make a duplicate to prevent Set changed size during iteration errors.
                    removeList.extend(toRemove)  # Mark all the removed nodes as removed.
                    neighbours[nClique].remove(i)  # Node i no longer has nClique neighbours.
                    neighbours[0].add(i)  # Node i now has no neighbours.
                    adjList[i] = set([])  # Update the adjacency list to reflect the fact that i has no neighbours.
                    neighboursOfRemoved = set([k for j in toRemove for k in adjList[j]])
                    for j in neighboursOfRemoved:
                        neighboursRemoved = adjList[j].intersection(toRemove)
                        numNeighbours = len(adjList[j])
                        neighbours[numNeighbours].remove(j)
                        neighbours[numNeighbours - len(neighboursRemoved)].add(j)
                        adjList[j] -= neighboursRemoved
                    # Update the adjacency list to reflect the removal of the nodes in toRemove.
                    for j in toRemove:
                        numNeighbours = len(adjList[j])
                        neighbours[numNeighbours].remove(j)
                        neighbours[0].add(j)
                        adjList[j] = set([])
                    nClique = 1
                    break
            else:
                # No clique found.
                nClique += 1

        ########################################
        # Perform the NeighbourCull operation. #
        ########################################
        # Re-calculate this, as it may have changed since it was last calculated.
        maxNeighbours = max(neighbours)
        while maxNeighbours > 0 and (not neighbours[maxNeighbours]):
            del neighbours[maxNeighbours]
            maxNeighbours -= 1

        # If there are no nodes with neighbours then exit.
        if maxNeighbours == 0:
            return removeList

        # Get the IDs of the nodes with the max number of neighbours.
        nodesWithMaxNeighbours = list(neighbours[maxNeighbours])
        # If there is more than one node with the maximum number of neighbours determine which node to remove.
        if len(nodesWithMaxNeighbours) != 1:
            # Determine the number of neighbours for each node.
            extendedNeighbourhood = [adjList[x].union([x]) for x in nodesWithMaxNeighbours]
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
            neighbours[numNeighbours].remove(i)
            neighbours[numNeighbours - 1].add(i)
            adjList[i].remove(toRemove)
        # Update the adjacency list to reflect the removal of to remove.
        adjList[toRemove] = set([])
        neighbours[maxNeighbours].remove(toRemove)
        neighbours[0].add(toRemove)