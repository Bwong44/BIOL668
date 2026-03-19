### Part I
# (1) Create a base class call Seq. Seq should have a constructor in initialize
#     and instance of Seq with sequence, gene name and species.
# (2) Create 2 methods for Seq.
#        1. print_record returns the sequence
#        2. overload the str function to return species, gene name : sequence
# 
# Test with the function calls below

### Part II
# (1) Create a class DNA that inherits Seq. In the constructor, use the super()
#        function to make the contructor. Add a variable gene_id and
#        add **kwargs also in the constructor
# (2) Add a method called analysis that returns then number of Gs and Cs
#
# Test with the function calls below


class Seq:
    def __init__(self, sequence, gene_name, species):
        self.sequence = sequence
        self.gene_name = gene_name
        self.species = species
    
    def print_record(self):
        return self.sequence

    def __str__(self): #Overload str function to return itself
        return f"{self.sequence}, {self.gene_name} : {self.species}"

#class DNA(Seq):


if __name__ == "__main__":
    myseq=Seq("ATATAG","my_gene","H.sapiens")
    print(myseq.print_record())
    print(myseq.gene_name)
    print(myseq.species)

#d=DNA("GATCTC","my_dna","D.terebrans","AX5667.2")
#d.print_record()
#print(d)

#d.source="Mexico"
#print(d.source)
#print(d.analysis())
