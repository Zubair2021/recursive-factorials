#! /usr/bin/python3

# Importing modules
from Bio import SeqIO
from Bio.Seq import translate
import sys
import re

sys.setrecursionlimit(1500)


# Defining recursive factorial formula

def get_factorial_recur(n):
    if n == 1:
        return 1
    else:
        return n * get_factorial_recur(n-1)


# Creating a class containing important attributes

class Rec:
     def __init__(self, ID, Date, Phenotype, Sequence):

         self.ID = ID
         self.Date = Date
         self.Phenotype = Phenotype
         self.Sequence = Sequence


# Creating an empty dictionary

record_dict = {"ID":[],"Date":[], "Phenotype":[], "Sequence":[]};

# Creating a for loop to parse information from fasta file and append to dictionary
def main():

  print('Opening' + fasta)
  with open(fasta, 'r') as in_stream:
        print('opening' + results)
        with open(results, 'w') as out_stream:
          for seq_record in SeqIO.parse(fasta, "fasta"):
            string = seq_record.id
            sample = string.split('_')[0]    # to split the name and date written together with _
            date = string.split('_')[1]
            sequence = seq_record.seq
            AA = translate(sequence)
            if AA[3] == 'R':
                          Phenotype = "Orange"
                          record_dict["Phenotype"].append(Phenotype)
            if AA[3] == 'S':
                          Phenotype = "Blue"
                          record_dict["Phenotype"].append(Phenotype)
            record_dict["ID"].append(sample)
            record_dict["Date"].append(date)
            record_dict["Sequence"].append(sequence)
            print(record_dict)
          print("Ending with of outstream")

          k = 0
          for i in record_dict['Phenotype']:
              if i == 'Orange':
                 k +=1
          print(k)

          p = float(sys.argv[2])
          print(type(p))
          q = 1-p
		##  Using the formula
          n = len([1 for line in open(sys.argv[1]) if line.startswith(">")])

          bern_formula = (get_factorial_recur(n)/(get_factorial_recur(n-k)*get_factorial_recur(k)))*((p**k)*(q**(n-k)))

          out_stream.write("Results\n\np (the frequency of \"orange\" in the population) = " + str(p))
          out_stream.write("\nn (the number of sampled individuals) = " + str(n))
          out_stream.write("\nk (the number of \"orange\" individuals in the sample set) = " + str(k))
          out_stream.write("\n\nProbability of collecting 32 individuals with 5 being \"orange\" (given a population frequency of 0.3) = " + str(bern_formula) + "\n")

        print("Done!")
        print(sys.argv[3] + ' is closed?', out_stream.closed)
  print(sys.argv[1] + ' is closed?', in_stream.closed)


if __name__ == '__main__':
    fasta = (sys.argv[1])
    p = float(sys.argv[2])
    results = (sys.argv[3])
    main()
