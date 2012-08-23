'''
Created on 3 Feb 2011

@author: Simon Bull
'''

def main(PSIoutput, processedLoc):
    """Extracts the relevant information from the PSI-BLAST output.
    
    @param PSIoutput: The location of the file containing the PSI-BLAST results.
    @type PSIoutput:  string
    @param processedLoc: The location to write the results of the processing to. If this location already exists
                         then the results will be appended to the existing file. If not the file will be created.
    @type processedLoc:  string
    
    """
    
    queryProtein = ''  # Records the protein used to query the database
    previousRoundOutput = ''  # Record the information to output as a string
    currentRoundOutput = ''  # Record the information to output as a string
    previousRoundHits = {}  # Record the hits from the previous round to prevent drift in the PSSM
    currentRoundHits = {}  # Record the hit proteins from the current round to prevent drift
    driftFound = False  # Used to exit and write out the results if drift is found
    
    BLASTOutput = open(PSIoutput, 'r')
    
    for line in BLASTOutput:
        
        chunks = line.split()

        if len(chunks) == 0:
            # If the line is a blank line, then ignore it.
            continue
        elif chunks[0] == '#' and chunks[1] == 'Query:':
            # If the current line is a comment line with the query protein information on it.
            for key in previousRoundHits:
                # Compare currentRoundHits with previousRoundHits to check for drift.
                if not currentRoundHits.has_key(key):
                    # If a hit is recorded in previousRoundHits but not in currentRoundHits then drift is deemed
                    # to have occurred.
                    currentRoundOutput = previousRoundOutput
                    driftFound = True
            if driftFound:
                break
            queryProtein = chunks[2]
            previousRoundHits = currentRoundHits
            currentRoundHits = {}
            previousRoundOutput = currentRoundOutput
            currentRoundOutput = ''
        elif chunks[0] == queryProtein:
            # An alignment is recorded on the line if the line starts with the query protein.
            hit = chunks[1]
            currentRoundHits[hit] = True
            currentRoundOutput += line

    BLASTOutput.close()
    
    outputResults = open(processedLoc, 'a')
    outputResults.write(currentRoundOutput)
    outputResults.close()