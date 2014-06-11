'''
Created on 3 Feb 2011

@author: Simon Bull
'''

def main(PSIoutput, minAlignLength, maxEValue):
    """Extracts the relevant information from the PSI-BLAST output.

    :param PSIoutput:       The location of the file containing the PSI-BLAST results.
    :type PSIoutput:        string
    :param minAlignLength:  The minimum permissible length for the BLAST sequence alignments.
    :type minAlignLength:   integer
    :param maxEValue:       The maximum permissible value which the BLAST EValue can take.
    :type maxEValue:        float
    :returns :              A record of the similarities between pairs of proteins.
    :type :                 dictionary

    """

    currentQuery = ''
    hitsFound = {}
    similaritiesFound = {}

    BLASTOutput = open(PSIoutput, 'r')
    for line in BLASTOutput:
        chunks = line.split()
        if len(chunks) == 0:
            # If the line is a blank line, then ignore it.
            continue
        elif chunks[0] == '#' and chunks[1] == 'Query:':
            # The end of a round has been reached.
            nextQuery = chunks[2]
            if currentQuery != nextQuery:
                # A new query has been found. Write out the hits from the last query.
                for hit in hitsFound:
                    if hit != currentQuery:
                        # Only record the similarity if the query and hit are not the same
                        pair = tuple(sorted([currentQuery, hit]))
                        if not (pair in similaritiesFound and similaritiesFound[pair] >= hitsFound[hit]):
                            # If the pair exists and the recorded similarity is less than the newly found one or the pair does not exist,
                            # then record the new value for the similarity.
                            similaritiesFound[pair] = hitsFound[hit]
            currentQuery = nextQuery
            hitsFound = {}
        elif chunks[0] == currentQuery:
            # An alignment is recorded on the line if the line starts with the query protein.
            hit = chunks[1]
            alignLength = int(chunks[3])
            evalue = float(chunks[4])
            if alignLength >= minAlignLength and evalue <= maxEValue:
                # Only record the hit if the alignment length is long enough and the evalue is large enough.
                hitsFound[hit] = float(chunks[2])
    BLASTOutput.close()

    return similaritiesFound