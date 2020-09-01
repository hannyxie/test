f = open('IRGSP-1.0_genome.fasta','r')
genomes = f.read().split('>')
for genome in genomes[1:]:
        print(genome.count('A')+genome.count('T')+genome.count('C')+genome.count('G')+genome.count('a')+genome.count('t')+genome.count('c')+genome.count('g')+genome.count('N'))