'''
Created on 7 Mar 2011

@author: Simon Bull
'''

from scipy import zeros

class SparseMatrix():
    """Class for storing and manipulating sparse matrices as a dictionary of lists.
    
    As the class is designed to represent adjacency matrices it is intended to work with square matrices, but
    can be used with non-square matrices (if there are more columns than rows this becomes difficult).
    Additionally, the class is designed to represent graphs with unweighted edges, and therefore only records
    the presence or absence of an edge.

    """

    def __init__(self, arg):
        """Initialise the sparse matrix.
        
        @param arg: arg must be either an integer or a dictionary of lists. If an integer is provided then arg is assumed
                    to be the number of rows in the matrix. If a dictionary of lists is provided the formatting of the
                    dictionary and lists must be the same as a dictionary created from scratch using the class, otherwise
                    subsequent methods will fail.
        @type arg:  integer or dictionary of lists
        
        """
        
        if isinstance(arg, int):
            self.dict = dict((x,[]) for x in range(arg))
        else:
            self.dict = arg

    def display(self):
        """Display the matrix."""
        
        for i in self.dict:
            print i, '\t', self.dict[i]

    def add(self, xDim, yDim):
        """Adds entry [xDim, yDim] where xDim and yDim are indices.
        
        @param xDim: Index of the row at which to add the entry.
        @type xDim:  integer
        @param yDim: Index of the column at which to add the entry.
        @type yDim:  integer
        
        """
        
        self.dict[xDim].append(yDim)

    def addlist(self, xDims, yDims):
        """Add multiple edges at once.
        
        The two lists must be of the same length. xDims[i] and yDims[i] correspond to the x and y indices for the
        ith edge to add.
        
        @param xDims: The rows at which to add entries.
        @type xDims:  list of integers
        @param yDims: The columns at which to add entries.
        @type yDims:  list of integers
        
        """
        
        for i in range(len(xDims)):
            self.dict[xDims[i]].append(yDims[i])

    def remove(self, xDim, yDim):
        """Removes entry [xDim, yDim] where xDim and yDim are indices.
        
        @param xDim: Index of the row at which to remove the entry.
        @type xDim:  integer
        @param yDim: Index of the column at which to remove the entry.
        @type yDim:  integer
        
        """
        
        self.dict[xDim].remove(yDim)

    def removeList(self, xDims, yDims):
        """Remove multiple edges at once.
        
        The two lists must be of the same length. xDims[i] and yDims[i] correspond to the x and y indices for the
        ith edge to remove.
        
        @param xDims: The rows at which to remove entries.
        @type xDims:  list of integers
        @param yDims: The columns at which to remove entries.
        @type yDims:  list of integers
        
        """
        
        for i in range(len(xDims)):
            self.dict[xDims[i]].remove(yDims[i])

    def take(self, dimList, column = True):
        """Returns a sparse matrix with only the rows or columns specified.
        
        @param dimList: The dimensions which should be returned.
        @type dimList:  list
        @param column: True indicates that dim~List corresponds to the columns to return, False means return rows.
        @type column:  boolean
        return @type: SparseMatrix
        return @use:  A SparseMatrix that contains only the rows/columns specified.
        
        """
        
        if not column:
            result = {}
            for x in self.dict.keys():
                y = [y for y in self.dict[x] if y in dimList]
                result[x] = y
        else:
            result = dict((x,self.dict[x]) for x in dimList)
        return SparseMatrix(result)

    def takesquare(self, dimensions):
        """Returns a sparse matrix with only those dimensions specified.
        
        Returns the same result as two calls to take where column is True in one call and False in the other,
        and dimList is kept the same for both calls. E.g.:
        rows = sparse.take(dimList, True)
        subset = rows.take(dimList, False)
        
        @param dimensions: The rows and columns to return.
        @type dimensions:  list
        return @type: SparseMatrix
        return @use:  A SparseMatrix containing the subset of rows and columns specified.
        
        """
        
        result = {}
        for y in dimensions:
            yAdd = [i for i in self.dict[y] if i in dimensions]
            result[y] = yAdd
        return SparseMatrix(result)

    def todense(self):
        """Returns a dense 2D scipy array of the sparse matrix.

        return @type: A 2 dimensional SciPy array, a list
        return @use:  An adjacency matrix of the SparseMatrix, the indices in the matrix

        """
        
        xVals = self.dict.keys()
        xValsIndices = dict((xVals[i], i) for i in range(len(xVals)))
        result = zeros((len(xVals), len(xVals)))
        for i in xVals:
            for j in self.dict[i]:
                result[xValsIndices[i],xValsIndices[j]] = 1
        return result, xVals

    def connectedcomponents(self):
        """Return a list of the connected components of the sparse matrix.
        
        Performs a breadth first search to determine the components.

        return @type: list
        return @use:  Each element contains the integer indices of the nodes of the graph that are in the component.
                      One element of the list per connected component.
        
        """
        
        subgraphs = []
        inSub = []
        toCheck = range(len(self.dict.keys()))
        while toCheck != []:
            subgraph = []
            start = toCheck.pop(0)
            subgraph.append(start)
            inSub.append(start)
            check = [x for x  in self.dict[start] if x in toCheck]
            while check != []:
                current = check.pop(0)
                toCheck.remove(current)
                subgraph.append(current)
                inSub.append(current)
                nonzero = self.dict[current]
                for i in nonzero:
                    if i in toCheck and i not in check:
                        check.append(i)
            subgraphs.append(subgraph)

        return subgraphs

    def adjList(self):
        """Returns an adjacency list with the indices altered to be in the range 0..len(self.dict).

        return @type: dictionary
        return @use:  a copy of the SparseMatrix with the indices of the matrix altered to be within the range 0..n where n is the number of nodes.

        """
        
        keys = self.dict.keys()
        keys = sorted(keys)
        indices = dict((keys[i], i) for i in range(len(keys)))
        result = {}
        for i in keys:
            tempList = [indices[x] for x in self.dict[i]]
            result[indices[i]] = tempList
        return result