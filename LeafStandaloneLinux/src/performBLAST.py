'''
Created on 7 Feb 2011

@author: Simon Bull
'''

import os
import subprocess
import shutil
import time

import processPSIoutput

def main(inputFile, databaseFile, blastOperationID, SEG=False, cores=2, minAlignLength=20, maxEValue=1.0, verboseOutput=False):
    """Perform the BLASTing of the proteins in an input file (inputFile) against those in another file (databaseFile).
    
    Returns a dictionary of the similarities between proteins, as determined by BLAST. The dictionary is indexed by a
    alphanumerically ordered tuple (index[0] < index[1]), and the entry for each index is a dictionary recording
    {'Identity' : identity, 'Length' : alignLength, 'EValue' : evalue}.
    
    The BLAST version used must be the C++ version. If a different version is being used (i.e. the old C version), then
    parameters like the number of cores can not be used.
    
    @param inputFile: The location of a FASTA format file of the proteins to BLAST against the proteins in databaseFile.
    @type inputFile : string
    @param databaseFile: The location of a FASTA format file of the proteins from which the BLASTable database should be generated.
    @type databaseFile : string
    @param blastOperationID: The name for the directory where the results of the BLASTing -> parsing will be stored.
    @type blastOperationID : string
    @param SEG: Whether or not SEG should be run on the input proteins.
    @type SEG : boolean
    @param cores: The number of CPU cores on which BLAST will be run.
    @type cores : integer
    @param minAlignLength: The minimum permissible length for the BLAST sequence alignments.
    @type minAlignLength : integer
    @param maxEValue: The maximum permissible value which the BLAST EValue can take.
    @type maxEValue : float
    @param verboseOutput: Whether status updates of the BLASTing should be printed out to the user.
    @type verboseOutput:  boolean
    return @type: dictionary
    return @use:  A record of the similarities between the proteins
    
    """
    
    # Get the location of the BLAST executables.
    srcLocation = os.path.abspath(__file__)
    srcLocation = '/'.join(srcLocation.split('/')[:-1])
    BLASTExecutables = srcLocation + '/BLASTExecutables'
    cwd = os.getcwd()
    outputLocation = cwd + '/' + blastOperationID
    if os.path.exists(outputLocation):
        shutil.rmtree(outputLocation)
    os.mkdir(outputLocation)
    
    # Make a BLASTable database from the database file.
    if verboseOutput:
        print 'Creating the BLASTable database.'
    databaseDir = outputLocation + '/TempDatabase'
    os.mkdir(databaseDir)
    os.mkdir(databaseDir + '/TempDB')
    makeDBArgs = BLASTExecutables + '/makeblastdb -in ' + databaseFile + ' -out ' + databaseDir + '/TempDB -dbtype prot'
    subprocess.call(makeDBArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    
    # Loop through the input file and create a FASTA format file for each individual protein.
    if verboseOutput:
        print 'Generating a FASTA file of each sequence.'
    proteinDir = outputLocation + '/TempProteins'
    os.mkdir(proteinDir)
    fullFASTA = open(inputFile, 'r')
    protCount = 0
    for line in fullFASTA:
        if line[0] == '>':
            # If the line starts a new protein definition.
            if protCount == 0:
                # If this is the first protein definition found.
                proteinWrite = open(proteinDir + '\Prot' + str(protCount) + '.fasta', 'w')
                proteinWrite.write(line)
            else:
                # If this is not the first protein definition found.
                proteinWrite.close()
                proteinWrite = open(proteinDir + '\Prot' + str(protCount) + '.fasta', 'w')
                proteinWrite.write(line)
            protCount += 1
        else:
            # Otherwise the line is a protein sequence.
            proteinWrite.write(line)
    
    proteinWrite.close()
    fullFASTA.close()
    
    # BLAST each of the individual protein FASTA files just made against the database generated from databaseFile.
    if verboseOutput:
        print 'Starting to BLAST each file.'
        fileCount = 1
    processedBLAST = outputLocation + '/Processed.txt'
    proteinFiles = os.listdir(proteinDir)
    for file in proteinFiles:
        if verboseOutput:
            if fileCount % 100 == 0:
                print 'Currently BLASTing file ', fileCount, ' out of ', len(proteinFiles), '...'
            fileCount += 1
        sequence_BLAST(processedBLAST, proteinDir + '/' + file, databaseDir + '/TempDB', BLASTExecutables + '/psiblast',
                       SEG, cores)
    
    # Parse the processed BLAST output, and record the similarities between the different proteins.
    if verboseOutput:
        print 'Now parsing the processed BLAST output.'
    similarities = {}
    readProcessedBLAST = open(processedBLAST, 'r')
    for line in readProcessedBLAST:
        chunks = line.split('\t')
        key = tuple(sorted([chunks[0], chunks[1]]))
        identity = float(chunks[2])
        alignLength = int(chunks[3])
        if alignLength <= minAlignLength:
            # If the alignment length is too short, then ignore the alignment.
            continue
        evalue = float(chunks[4])
        if evalue >= maxEValue:
            # If the EValue is too great, then ignore the alignment.
            continue
        if similarities.has_key(key):
            oldSimilarity = similarities[key]['Identity']
            if identity > oldSimilarity:
                similarities[key] = {'Identity' : identity, 'Length' : alignLength, 'EValue' : evalue}
        else:
            similarities[key] = {'Identity' : identity, 'Length' : alignLength, 'EValue' : evalue}
    readProcessedBLAST.close()

    # Remove the temporary directory used for manipulating and processing the BLAST output.
    try:
        shutil.rmtree(outputLocation)
    except:
        time.sleep(60)
        shutil.rmtree(outputLocation)
    
    return similarities


def sequence_BLAST(processedBLAST, inputFile, database, BLASTLoc, SEG, cores):
    """Will perform the process of BLAST -> PROCESS OUTPUT on inputFile.
    
    @param processedBLAST: The location at which to write the output of the processing of the BLAST output.
    @type processedBLAST: string
    @param inputFile: The FASTA file which needs to be submitted to PSI-BLAST.
    @type inputFile: string
    @param database: The database to BLAST the inputFile protein against
    @param database: string
    @param BLASTLoc: The location of the PSI-BLAST executable.
    @type BLASTLoc: string
    @param SEG: Set to True to use SEG to mask low complexity regions of the query.
    @type SEG: boolean
    @param cores: The number of threads to create to run BLAST with.
    @type cores: character
    
    """ 

    # Setup the parameters for the BLASTing.
    outputLoc = inputFile.split('.')[0] + '.tmp'    
    query = ' -query ' + inputFile
    out = ' -out ' + outputLoc
    evalue = ' -evalue 1'
    inclusionEThresh = ' -inclusion_ethresh 0.0001'
    numIterations = ' -num_iterations 3'
    gapTrigger = ' -gap_trigger 18'
    numDescriptions = ' -num_descriptions 10000'
    numAlignments = ' -num_alignments 10000'
    dbsize = ' -dbsize 0'
    db = ' -db ' + database
    outputFormat = ' -outfmt "7 qseqid sseqid pident length evalue"'
    if SEG:
        seg = ' -seg yes'
    else:
        seg = ' -seg no'
    numThreads = ' -num_threads ' + str(cores)
    argsPSI = (query + out + evalue + inclusionEThresh + numIterations + gapTrigger + numDescriptions +
               numAlignments + dbsize + db + outputFormat + seg + numThreads
               )
    # Perform the BLASTing.
    subprocess.call(BLASTLoc + argsPSI, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # Process the BLAST output.
    processPSIoutput.main(outputLoc, processedBLAST)