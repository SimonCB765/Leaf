'''
Created on 5 Feb 2011

@author: Simon Bull
'''

def main(listToCheck, allChains={}, allEntries=[], allProtEntries=[], checkType='chain'):
    """Checks whether the chains/entries provided as input are valid protein chains or entries that contain a protein chain.
    
    @param listToCheck: The list of chains/entries to check.
    @type listToCheck : list
    @param allChains: A dictionary mapping each chain in the PDB to its type (e.g. Protein, DNA, RNA, etc.). Only needed when checkType == 'chain'.
    @type allChains : dictionary
    @param allEntries: A list of all the entries in the PDB. Only needed when checkType == 'entry'.
    @type allEntries : list
    @param allProtEntries: A list of all the entries in the PDB that contain valid proteins. Only needed when checkType == 'entry'.
    @type allProtEntries : list
    @param checkType: 'chain' if chains are being checked, 'entry' if entries are being checked
    @type checkType : string
    return @type: integer, unicode string
    return @use:  the numerical code for the error, the error message or the valid chains/entries
    
    """

    listToCheck = set([i.strip().upper() for i in listToCheck])  # Remove any extra whitespace, and convert all chains/entries to uppercase.
    listToCheck -= set([''])  # Remove the record of any blank lines.

    if checkType == 'chain':
        # If the validity of chains is being checked.
        nonExistantChains = [i for i in listToCheck if i not in allChains.keys()]  # Chains in the input that can't be found.
        nonProtChains = [i for i in listToCheck if not i in nonExistantChains and allChains[i] != 'Protein']  # Chains that can be found, but aren't proteins.
        if nonExistantChains != [] or nonProtChains != []:
            # If there are any chains that don't exist or aren't proteins, then exit with an error message.
            if nonExistantChains != []:
                errorMessage = u'The following chains were not found in the database:\n' + u', '.join(nonExistantChains) + u'\n'
            if nonProtChains != []:
                errorMessage = u'The following chains are not known to be protein chains:\n' + u', '.join(nonProtChains) + u'\n'
            return 1, errorMessage
        elif len(listToCheck) < 2:
            # If there aren't enough chains in the input.
            errorMessage = u'The input file does must contain at least 2 chains.\n'
            return 1, errorMessage
        else:
            # All the chains are valid protein chains, and there are enough of them.
            return 0, '\n'.join(listToCheck)
    elif checkType == 'entry':
        # If the validity of entries is being checked.
        nonExistantEntries = [i for i in listToCheck if not i in allEntries]  # Entries that can;t be found.
        entriesWithNoProts = [i for i in listToCheck if not i in allProtEntries and not i in nonExistantEntries]  # Entries that can be found, but have no protein chains.

        if nonExistantEntries != [] or entriesWithNoProts != []:
            # If there are any input entries that either don;t exist or don't have any protein chains, then exit with
            # an error message.
            if nonExistantEntries != []:
                errorMessage = u'The following entries were not found in the database:\n' + u', '.join(nonExistantEntries) + u'\n'
            if entriesWithNoProts != []:
                errorMessage += (u'The following entries were found in the database, but have no protein chains recorded for them:\n' +
                                 u', '.join(entriesWithNoProts) + u'\n')
            return 1, errorMessage
        elif len(listToCheck) < 2:
            # If there aren't enough entries in the input.
            errorMessage = u'The input file does must contain at least 2 entries.\n'
            return 1, errorMessage
        else:
            # All the entries are valid protein chain containing entries, and there are enough of them.
            return 0, '\n'.join(listToCheck)
    else:
        return 2, u'The choice for parameter checkType was not one of \'chain\' or \'entry\'.'