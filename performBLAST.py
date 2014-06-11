'''
Created on 7 Feb 2011

@author: Simon Bull
'''

import os
import subprocess
import shutil
import time

import processPSIoutput

def main(inputFile, blastOperationID, cores=2, minAlignLength=20, maxEValue=1.0, verboseOutput=False):
    """Perform the BLASTing of the proteins in an input file (inputFile) against those in another file (databaseFile).

    Returns a dictionary of the similarities between proteins, as determined by BLAST. The dictionary is indexed by a
    alphanumerically ordered tuple (index[0] < index[1]), and the entry for each index is the percentage sequence similarity.

    The BLAST version used must be the C++ version. If a different version is being used (i.e. the old C version), then
    parameters like the number of cores can not be used.

    :param inputFile:           The location of a FASTA format file of the proteins to BLAST against each other.
    :type inputFile:            string
    :param blastOperationID:    The name for the directory where the results of the BLASTing -> parsing will be stored.
    :type blastOperationID:     string
    :param cores:               The number of CPU cores on which BLAST will be run.
    :type cores:                integer
    :param minAlignLength:      The minimum permissible length for the BLAST sequence alignments.
    :type minAlignLength:       integer
    :param maxEValue:           The maximum permissible value which the BLAST EValue can take.
    :type maxEValue:            float
    :param verboseOutput:       Whether status updates of the BLASTing should be printed out to the user.
    :type verboseOutput:        boolean
    :returns :                  A record of the similarities between the proteins
    :type :                     dictionary

    """

    # Get the location of the BLAST executables.
    srcLocation = os.path.dirname(os.path.realpath(__file__))
    BLASTExecutables = os.path.join(srcLocation, 'BLASTExecutables')
    cwd = os.getcwd()
    outputLocation = blastOperationID
    if os.path.exists(outputLocation):
        shutil.rmtree(outputLocation)
    os.mkdir(outputLocation)

    # Generate BLAST database.
    databaseDir = outputLocation + '/TempDatabase'
    os.mkdir(databaseDir)
    makeDBArgs = [BLASTExecutables + '/makeblastdb', '-in', inputFile, '-out', databaseDir + '/TempDB', '-dbtype', 'prot']
    subprocess.call(makeDBArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    # Perform BLAST.
    if verboseOutput:
        print('Now BLASTing.')
    resultsBLAST = outputLocation + '/ResultsBLAST.txt'
    sequence_BLAST(resultsBLAST, inputFile, databaseDir + '/TempDB', BLASTExecutables + '/psiblast', cores)

    # Determine the similarities between proteins.
    if verboseOutput:
        print('Now determining similarities.')
    similarities = processPSIoutput.main(resultsBLAST, minAlignLength, maxEValue)

    # Remove the temporary directory used for manipulating and processing the BLAST output.
    try:
        shutil.rmtree(outputLocation)
    except:
        time.sleep(60)
        shutil.rmtree(outputLocation)

    return similarities


def sequence_BLAST(resultsBLAST, inputFile, database, BLASTLoc, cores):
    """Will perform the process of BLAST -> PROCESS OUTPUT on inputFile.

    :param resultsBLAST:    The location at which to write the output of the processing of the BLAST output.
    :type resultsBLAST:     string
    :param inputFile:       The FASTA file which needs to be submitted to PSI-BLAST.
    :type inputFile:        string
    :param database:        The database to BLAST the inputFile protein against
    :type database:         string
    :param BLASTLoc:        The location of the PSI-BLAST executable.
    :type BLASTLoc:         string
    :param cores:           The number of threads to create to run BLAST with.
    :type cores:            character

    """

    # Setup the parameters for the BLASTing.
    argsPSI = []
    argsPSI.append(BLASTLoc)
    argsPSI.append('-query')
    argsPSI.append(inputFile)
    argsPSI.append('-out')
    argsPSI.append(resultsBLAST)
    argsPSI.append('-evalue')
    argsPSI.append('1')
    argsPSI.append('-num_iterations')
    argsPSI.append('3')
    argsPSI.append('-gap_trigger')
    argsPSI.append('18')
    argsPSI.append('-num_descriptions')
    argsPSI.append('10000')
    argsPSI.append('-num_alignments')
    argsPSI.append('10000')
    argsPSI.append('-dbsize')
    argsPSI.append('0')
    argsPSI.append('-db')
    argsPSI.append(database)
    argsPSI.append('-outfmt')
    argsPSI.append('7 qseqid sseqid pident length evalue')
    argsPSI.append('-num_threads')
    argsPSI.append(str(cores))

    # Perform the BLASTing.
    subprocess.call(argsPSI, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
