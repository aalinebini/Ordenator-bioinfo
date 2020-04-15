import sys
import getopt

from Bio import SeqIO
import pandas as pd
import re

class Ordenator():
    """A class to ordenate DNA sequences by their primers
    
    Returns:
        [CSV] -- A CSV with the ordenated sequneces
    """

    def __init__(self):
        self.primers = None
        self.sequences = None

    def read_sequence(self, sequence_path):
        """this method loads a DNA sequence with a given input fastq file
        
        Arguments:
            sequence_path {str} -- the path of the fastq file
        """

        description = list()
        sequence = list()

        for record in SeqIO.parse(sequence_path, "fastq"):
    
            description.append(record.description)
            sequence.append(record.seq)

        self.sequences = pd.DataFrame({'Description': description, 'Sequences': sequence})
        self.sequences['Sequences'] = self.sequences.Sequences.apply(str)

    def read_primers(self, primers_path):
        """This method reads a primers list from a CSV or TXT file
        
        Arguments:
            primers_path {str} -- the path of the primers list file
        """

        self.primers = pd.read_csv(primers_path)

    def catalog(self, sequence):
        """A function to map the primer of a given sequence. 
           The values os the primers are between 0 and 29,
           where -1 represents a primer mismatch
        
        Returns:
            int -- an int that representes the primer type
        """

        for i, primer in enumerate(self.primers.Primers):
            if re.match('^%s.*' % primer, sequence):
                return i
            
        return -1

    def ordenate_sequence(self):
        """A method responsable for create a new column on the sequences Dataframe,
           and then ordenate that column based on the primers numbers
        """
        
        self.sequences['primer_number'] = self.sequences.Sequences.apply(self.catalog)
        self.sequences.sort_values(by="primer_number", inplace=True)

    def save_to_csv(self, path_to_save):
        """A method to save the ordenated sequence into a CSV file
        
        Arguments:
            path_to_save {str} -- the path to save the ordenated sequence
        """
        
        self.sequences.to_csv(path_to_save, index=False)

if __name__ == "__main__":
    
    try:
        OPTS, ARGS = getopt.getopt(sys.argv[1:], 'f:p:s:h', ['fastq_path', 
                                                             'primer_path', 
                                                             'save_path', 
                                                             'help'])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(1)

    ORDENATOR = Ordenator()
    SEQUENCE_PATH = None
    PRIMER_PATH = None
    SAVE_PATH = None

    for opt, arg in OPTS:

        if opt in ('-h', '--help'):
            print('''
            This program sorts the sequences according to the primers. 
            ordenator.py -f --fastq_path -p --primer_path -s --save_path -h --help
           ''')
            sys.exit(2)
        
        elif opt in ('-f', '--fastq_path'):
            if re.match('.+fastq$|.+txt$', arg):
                SEQUENCE_PATH = arg
            else:
                print("-f isn't a fastq or txt file")
                sys.exit(3)

        elif opt in ('-p', '--primer_path'):
            if re.match('.+csv$', arg):
                PRIMER_PATH = arg
            else:
                print("-p isn't a csv file")
                sys.exit(4)

        elif opt in ('-s', '--save_path'):
            if re.match('.+csv$', arg):
                SAVE_PATH = arg
            else:
                print("-s isn't a csv file")
                sys.exit(5)

    if SEQUENCE_PATH and PRIMER_PATH and SAVE_PATH:

        ORDENATOR.read_sequence(SEQUENCE_PATH)
        ORDENATOR.read_primers(PRIMER_PATH)
        ORDENATOR.ordenate_sequence()
        ORDENATOR.save_to_csv(SAVE_PATH)

    else:
        print('missing arguments -f, -p or -s')
        sys.exit(6)