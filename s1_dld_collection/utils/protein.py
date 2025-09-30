from Bio.SeqUtils import IUPACData

def protein_letters_3to1(protein_3_letter: list):
    '''
    Convert 3 letters to 1 letter representation.
    
    params:
        protein_3_letter: list
    return:
        protein_1_letter: list
        
    '''
    protein_1_letter = []
    protein_letters_3to1 = {k.upper():v for k, v in IUPACData.protein_letters_3to1.items()}
    try:
        protein_1_letter = [protein_letters_3to1[x.upper()] for x in protein_3_letter]
    except:
        # with ubiquitous amino acids
        for letter in protein_3_letter:
            if letter in protein_letters_3to1:
                protein_1_letter.append(protein_letters_3to1[letter])
            else:
                protein_1_letter.append('X')
    
    return protein_1_letter