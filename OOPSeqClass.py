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


#class Seq:


#class DNA(Seq):


# myseq=Seq("ATATAG","my_gene","H.sapiens")
# print(myseq.print_record())
# print(myseq.gene)
# print(myseq.species)

#d=DNA("GATCTC","my_dna","D.terebrans","AX5667.2")
#d.print_record()
#print(d)

#d.source="Mexico"
#print(d.source)
#print(d.analysis())
