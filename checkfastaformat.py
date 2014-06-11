'''
Created on 2 Feb 2011

@author: Simon Bull
'''

def main(fileToCheck, minLength=-1, maxLength=-1):
    """Determines whether fileToCheck is an appropriately formatted FASTA file.

    An appropriately formatted FASTA file is returned in correctFormatting.

    FASTA files are accepted if they have the format:
    >PID1
    letters
    >PID2
    letters

    Where PID1 and PID2 can be anything, and letters are a (possibly multiline) sequence of alphabetic letters.
    The letters can be upper or lower case, and each letter is interpreted as one amino acid.
    The correctly formatted FASTA file is returned with the sequence only going over one line, and all letters
    in upper case.

    If a protein (i.e. a FASTA information line) appears in the file more than one time, then the final appearance is
    taken to be the correct one. Prior appearances are discarded.

    :param fileToCheck: A string containing the contents of a potential FASTA format file.
    :type fileToCheck:  string
    :param minLength:   The minimum length that a protein sequence in the FASTA file is permitted to be.
    :type minLength:    integer
    :param maxLength:   The maximum length that a protein sequence in the FASTA file is permitted to be.
    :type maxLength:    integer
    :returns :          The numerical code for the error and the error message or the valid FASTA format contents.
    :type :             integer and a unicode string

    """

    # Initialise variables.
    lineCount = 1  # The number of the line being examined. Used for displaying error messages.
    protDescription = True  # Whether or not we are currently expecting a line starting with >.
    firstLine = True  # Whether or not we are currently examining the first line of the file.
    proteinsInFile = {}  # A dictionary indexed by the protein description line of the FASTA file.
                         # The value of each entry is the correctly formatted protein sequence corresponding to the index.

    # Strip off all excess whitespace, and split the string into the individual lines of the file.
    checking = fileToCheck.rstrip()
    checking = checking.lstrip()
    checking = checking.split('\n')
    for line in checking:
        line = line.rstrip()
        if firstLine:
            # True if we have just started parsing the file string, and haven;t yet examined any lines.
            if line[0] == '>':
                currentProt = line  # Record the description line of the protein which is about to have its sequence inspected.
                currentSeq = ''  # Initialise the sequence of the protein.
                protDescription = False  # We are now expecting a protein sequence, not a protein description.
                firstLine = False
            else:
                # The first line of the file MUST be a protein description line (i.e. start with '>'). If the line was not
                # the beginning of a protein record, terminate the program.
                errorMessage = "Expected line " + str(lineCount) + " to start with a >, but instead got: " + line
                return 1, errorMessage
        elif protDescription:
            # This is true only if a line beginning with a '>' is expected.
            if line[0] == '>':
                # Expected a protein description line, and found a protein description line. This means that the entire sequence
                # of the currentProt protein has been found (i.e. we have finished inspecting the sequence of a protein, and
                # have found the protein to be valid). Now determine if the length of the sequence is within the user
                # specified bounds.
                if minLength == -1:
                    if maxLength == -1:
                        # If there are no restrictions on the protein sequence length, then record the protein and its sequence.
                        proteinsInFile[currentProt] = currentSeq
                    elif len(currentSeq) <= maxLength:
                        # If there is no minimum length restriction, and the protein sequence is not longer than the maximum
                        # sequence length permitted, then record the protein and its sequence.
                        proteinsInFile[currentProt] = currentSeq
                elif len(currentSeq) >= minLength:
                    if maxLength == -1:
                        # If there is no maximum length restriction, and the protein sequence is not shorter than the minimum
                        # sequence length permitted, then record the protein and its sequence.
                        proteinsInFile[currentProt] = currentSeq
                    elif len(currentSeq) <= maxLength:
                        # If the protein sequence is not shorter than the minimum sequence length permitted and not longer
                        # than the maximum length permitted, then record the protein and its sequence.
                        proteinsInFile[currentProt] = currentSeq
                currentProt = line  # Record the description line of the protein which is about to have its sequence inspected.
                currentSeq = ''  # Initialise the sequence of the protein.
                protDescription = False  # We are now expecting a protein sequence, not a protein description.
            else:
                # If the line does not begin with a '>', and it is expected to, it is possible that the amino acid sequence
                # is split over multiple lines.
                if line.isalpha():
                    # If every character on the line is a letter, then the line contains a valid portion of the sequence.
                    # Add the uppercase version of the sequence portion to the sequence currently being recorded.
                    currentSeq += line.upper()
                else:
                    # If the line did not contain only letters, terminate the program.
                    errorMessage = "Expected line " + str(lineCount) + " to start with a >, but instead got: " + line
                    return 1, errorMessage
        else:
            # If an amino acid sequence is expected.
            if line.isalpha():
                # If the line is all alphabetic characters, write the line out and indicate that we are expecting a
                # protein description line next (i.e. one beginning with a '>').
                currentSeq += line.upper()
                protDescription = True
            else:
                # If the line did not contain only letters, terminate the program.
                errorMessage = "Expected line " + str(lineCount) + " to contain only letters, but instead got: " + line
                return 2, errorMessage

        lineCount += 1

    # Catch the final protein from the file, and determine whether it should be recorded.
    if minLength == -1:
        if maxLength == -1:
            proteinsInFile[currentProt] = currentSeq
        elif len(currentSeq) <= maxLength:
            proteinsInFile[currentProt] = currentSeq
    elif len(currentSeq) >= minLength:
        if maxLength == -1:
            proteinsInFile[currentProt] = currentSeq
        elif len(currentSeq) <= maxLength:
            proteinsInFile[currentProt] = currentSeq

    if len(proteinsInFile.keys()) < 2:
        # There are too few protein sequences entered
        errorMessage = ("Not enough unique protein sequences have been entered." +
                        " This is possibly caused by not enough sequences of the required minimum and maximum length being provided."
                        )
        return 3, errorMessage
    elif protDescription:
        # Return an indication that the FASTA file is correctly formatted.
        outputString = ''
        for i in proteinsInFile.keys():
            outputString += i + '\n' + proteinsInFile[i] + '\n'
        return 0, outputString[:-1]
    else:
        # The file did not end with a protein sequence.
        errorMessage = "Reached the end of the file, but no protein sequence found for the final protein."
        return 3, errorMessage
