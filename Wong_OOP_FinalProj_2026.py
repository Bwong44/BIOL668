## RENAME this file YourLastName_OOP_FinalProject_2026.py

##Assignment: Add to the constructor and methods of a parent class and child classes
##            which inherit the base class properties. NOTE: You are not allowed
##            to import any specialized libraries for this project (e.g., no Biopython)
##            The idea is for you to write these methods from scratch.

## Begin with the parent Seq class and the child DNA class we created in lecture below.
## 

### Seq Class
#
#  Constructor:
#  (1) Use the string functions upper and strip to clean up self.sequence.
#  (2) Add a variable self.kmers to the constructor and make it equal to an empty list.

#  Methods:
#  (1) Add a method called make_kmers that makes overlapping kmers of a given length from self.sequence
#      appends these to self.kmers. Default kmer parameter=3.
#  (2) Add a method called fasta that returns a fasta formatted string like this:
#      >species gene
#      AGATTGATAGATAGATAT


### DNA Class: INHERITS Seq class
#   
#  Constructor:
#  Use re.sub to change any non nucleotide characters in self.sequence into an 'N'.
#      re.sub('[^ATGCU]','N',sequence) will change any character that is not a
#      capital A, T, G, C or U into an N. (Seq already uppercases and strips.)

#  Methods:
#  (1) Add a method called print_info that is like print_record, but adds geneid and an
#      empty space to the beginning of the string.
#  (2) Add a method called reverse_complement that returns the reverse complement of
#      self.sequence
#  (3) Add a method called six_frames that returns all 6 frames of self.sequence
#      This include the 3 forward frames, and the 3 reverse complement frames

### RNA Class:  INHERITS DNA class
#  
#  Construtor:
#  Use the super() function (see DNA Class example).
#  (1) Automatically change all Ts to Us in self.sequence. 
#  (2) Add self.codons equals to an empty list

#  Methods:
#  (1) Add make_codons which breaks the self.sequence into 3 letter codons
#      and appends these codons to self.codons unless they are less than 3 letters long.
#  (2) Add translate which uses the Global Variable standard_code below to
#      translate the codons in self.codons and returns a protein sequence.

### Protein Class: INHERITS Seq class
#
#  Construtor:
#  Use the super() function (see DNA Class example).
#  Use re.sub to change any non LETTER characters in self.sequence into an 'X'.

#  Methods:
#  The next 2 methods use a kyte_doolittle and the aa_mol_weights dictionaries.
#  (2) Add total_hydro, which return the sum of the total hydrophobicity of a self.sequence
#  (3) Add mol_weight, which returns the total molecular weight of the protein
#      sequence assigned to the protein object. 



import re

standard_code = {
     "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S",
     "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y",
     "UAA": "*", "UAG": "*", "UGA": "*", "UGU": "C", "UGC": "C",
     "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
     "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAU": "H",
     "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R",
     "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I", "AUA": "I",
     "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
     "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S",
     "AGC": "S", "AGA": "R", "AGG": "R", "GUU": "V", "GUC": "V",
     "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A",
     "GCG": "A", "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
     "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

kyte_doolittle={'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,'H':-3.2,'I':4.5,'K':-3.9,'L':3.8,
                'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,'R':-4.5,'S':-0.8,'T':-0.7,'V':4.2,'W':-0.9,'X':0,'Y':-1.3}

aa_mol_weights={'A':89.09,'C':121.15,'D':133.1,'E':147.13,'F':165.19,
                'G':75.07,'H':155.16,'I':131.17,'K':146.19,'L':131.17,
                'M':149.21,'N':132.12,'P':115.13,'Q':146.15,'R':174.2,
                'S':105.09,'T':119.12,'V':117.15,'W':204.23,'X':0,'Y':181.19}



class Seq:

    def __init__(self,sequence,gene,species):
        """Initialize the sequence class object with the given parameters.

        >>> s=Seq(" AGTt ","tmp","m")
        >>> s.sequence
        'AGTT'
        >>> s.gene
        'tmp'
        >>> s.species
        'm'
        """

        self.sequence=sequence.upper().strip()
        self.gene=gene
        self.species=species    
        self.kmers=[]

    def __str__(self): #Added by Brandon Wong
        """Overload str function to return the object's sequence.
        
        >>> s=Seq("AGT","tmp","m")
        >>> print(s)
        AGT
        """

        return self.sequence
    
    def __len__(self): #Adds len overload so we can access the seq length of the object - Brandon Wong
        """Len overload to return the length of the sequence from the Seq class.

        >>> s=Seq("AGTAGC","tmp","m")
        >>> len(s)
        6
        """

        return len(self.sequence)

    def __eq__(self, other): # adds eq overload to compare two seq objects - Troy Haynes
        """
        Comparing two seq objects by sequence.

        >>> Seq("ATGC", "gene1", "human") == Seq("ATGC", "gene2", "dog")
        True
        >>> Seq("ATGC", "gene1", "human") == Seq("AAAA", "gene2", "dog")
        False
        """
        return self.sequence == other.sequence

    def print_record(self):
        """Prints the record in the format 'species gene: sequence'.

        >>> s=Seq("AGT","test_gene","test_species")
        >>> s.print_record()
        test_species test_gene: AGT
        """
        
        print(self.species + " " + self.gene + ": " + self.sequence)

    def make_kmers(self, k=3): 
        """Makes overlapping kmers from a given sequence and appends them to self.kmer list.

        >>> s=Seq("AGTAGC","tmp","m")
        >>> s.make_kmers()
        >>> s.kmers
        ['AGT', 'GTA', 'TAG', 'AGC']
        """

        self.kmers=[]
        for i in range(len(self.sequence)):
            if len(self.sequence[i:i+k]) == k:
                self.kmers.append(self.sequence[i:i+k])
            
    def fasta(self):
        """Returns a fasta formatted string of the sequence.

        >>> s=Seq("AGT","test_gene","test_species")
        >>> print(s.fasta())
        >test_species test_gene
        AGT
        """

        fasta_output = ">" + self.species + " " + self.gene + "\n" + self.sequence
        return fasta_output

    def count_base(self, base): # counts the number of times a base is present in a sequence
        """
        Counts the occurrence of a specific base in sequence.

        >>> seq = Seq("AGAGTTT", "gene1", "human")
        >>> seq.count_base("T")
        3
        >>> seq.count_base("A")
        2
        """
        return self.sequence.count(base)

class DNA(Seq):

    def __init__(self,sequence,gene,species,geneid,**kwargs):
        """Initialize the DNA class object, inheriting from Seq class, adding geneid.

        >>> d=DNA("AGTXAGC","tmp","m","geneid_test")
        >>> d.sequence
        'AGTNAGC'
        >>> d.geneid
        'geneid_test'
        """

        super().__init__(sequence,gene,species) #Due to super the sequnece will already be upper and stripped
        self.sequence=re.sub("[^ATGCU]","N",self.sequence)
        self.geneid=geneid

    def analysis(self):
        """Returns length of G and C in the DNA sequence.

        >>> d=DNA("AGTXAGC","tmp","m","geneid_test")
        >>> d.analysis()
        3
        """

        gc=len(re.findall('G',self.sequence) + re.findall('C',self.sequence))
        return gc
    
    def print_info(self):
        """Inherits print_record and adds geneid to the beginning of the string.

        >>> d=DNA("AGTXAGC","tmp","m","geneid_test")
        >>> d.print_info()
        geneid_test m tmp: AGTNAGC
        """

        print(self.geneid + " ", end = "")
        super().print_record()

    def reverse_complement(self):
        """
        Returns reverse complement of sequence.

        >>> seq = DNA("ATGC","gene1","human",1)
        >>> seq.reverse_complement()
        'GCAT'
        """
        reverse=self.sequence[::-1]
        reverse=reverse.replace("A","t").replace("T","a").replace("G","c").replace("C","g").upper()
        return reverse

    def six_frames(self):
        """
        Returns the six reading frames of the sequence.

        >>> seq = DNA("ATGCCG","gene1","human",1)
        >>> seq.six_frames()
        ['ATGCCG', 'TGCCG', 'GCCG', 'CGGCAT', 'GGCAT', 'GCAT']
        """
        frames_list = []
        frames_list.append(self.sequence) #1st frame
        frames_list.append(self.sequence[1:]) #2nd frame by shifting 1
        frames_list.append(self.sequence[2:]) #3rd frame by shifting 2
        reverse_comp = self.reverse_complement()
        frames_list.append(reverse_comp) #4th frame is reverse complement
        frames_list.append(reverse_comp[1:]) #5th frame by shifting 1
        frames_list.append(reverse_comp[2:]) #6th frame by shifting 2
        return frames_list

    # calculates gc content of sequence - Troy Haynes
    # if percent True -> return percentage
    # if percent False -> return proportion
    def gc_content(self, percent=False):
        """
        Calculates GC content.

        >>> seq = DNA("ATGC","gene1","human",8417)
        >>> seq.gc_content()
        0.5
        >>> seq.gc_content(True)
        50.0
        """
        gc = self.sequence.count('G') + self.sequence.count('C')

        if percent == True:
            return (gc / len(self.sequence)) * 100
        else:
            return gc / len(self.sequence)

class RNA(DNA):
    
    def __init__(self,sequence,gene,species,geneid,**kwargs):
        """Initialize the RNA class object, inheriting from DNA class, changing Ts to Us and adding codons list.
        
        >>> r=RNA("AGTXAGC","tmp","m","geneid_test")
        >>> r.sequence
        'AGUNAGC'
        >>> r.geneid
        'geneid_test'
        """

        super().__init__(sequence,gene,species,geneid)
        self.sequence=self.sequence.replace("T","U")
        self.codons=[]
        
    def make_codons(self):
        """Returns the codon list of the sequence

        >>> r=RNA("AGTXAGC","tmp","m","geneid_test")
        >>> r.make_codons()
        ['AGU', 'NAG']
        """
        for i in range(0,len(self.sequence),3):
            if len(self.sequence[i:i+3]) == 3:
                self.codons.append(self.sequence[i:i+3])
        return self.codons
 
    def translate(self):
        """Translates the codons in self.codons into a protein sequence using the standard_code dictionary.
        >>> r=RNA("AGTXAGC","tmp","m","geneid_test")
        >>> r.make_codons()
        ['AGU', 'NAG']
        >>> r.translate()
        'SX'
        """

        protein_seq=""
        for codon in self.codons:
            if codon in standard_code:
                protein_seq += standard_code[codon]
            else:
                protein_seq += "X" 
        return protein_seq

class Protein(Seq):

    def __init__(self,sequence,gene,species,geneid,**kwargs):
        super().__init__(sequence,gene,species)
        self.sequence=re.sub("[^A-Z]","X",self.sequence)
        self.geneid=geneid

    def total_hydro(self):
        total_hydro_score=0
        for aa in self.sequence:
            if aa in kyte_doolittle:
                total_hydro_score += kyte_doolittle[aa]
        return total_hydro_score

    def mol_weight(self):
        mol_weight_score=0
        for aa in self.sequence:
            if aa in  aa_mol_weights:
                mol_weight_score += aa_mol_weights[aa]
        return mol_weight_score
    
    def amino_acid_composition(self): #Added aa comp dictionary - Brandon Wong
        """Creates a dictionary of the amino acids and their count in the sequence.
        >>> r=Protein("AGTXAGC","tmp","m","geneid_test")
        >>> r.amino_acid_composition()
        {'A': 2, 'G': 2, 'T': 1, 'X': 1, 'C': 1}
        """
        aa_composition={}
        for aa in self.sequence:
            if aa in aa_composition:
                aa_composition[aa] += 1
            else:
                aa_composition[aa] = 1
        return aa_composition

if __name__ == "__main__":
    import doctest
    doctest.testmod()