import sys

outlier_list = sys.argv[1]

annotated_vcf = sys.argv[2]

#mutation_type_list = []
#mutation_strength_list = []
concatenated_list = []
genome_position_list = []
with open(outlier_list, 'r') as outlierlist:

    for line2 in outlierlist: 
        if line2.startswith('#'): #ignores all the header lines

            continue

        else:

            line2 = line2.strip()

            items2 = line2.split('\t')
            # print(items2)
            contig = items2[0] #the #CHROM column in the VCF
            # print(contig)
            position = items2[1] #the POS column in the VCF

            concatenated = contig + "." + position #creates a string that contains both the chromosome and the position information
            
            genome_position = items2[2]
            #formatfield = items2[7]

            #things = formatfield.split(';')
            #ann = things[int(len(things)-1)]
            #print(ann)
            #details = ann.split('|')
            #print(details)
            #mutation_type = details[1]
            #print(mutation_type)
            #mutation_strength = details[2]
            # print(mutation_strength)
            #mutation_type_list.append(mutation_type)
            #mutation_strength_list.append(mutation_strength)
            concatenated_list.append(concatenated)
            genome_position_list.append(genome_position)
            
concatenated_and_genome_position = zip(concatenated_list, genome_position_list)
concatenated_dictionary = dict(concatenated_and_genome_position)

#concatenated_and_mutationstrength = zip(concatenated_list, mutation_strength_list)
#mutationstrength_dictionary = dict(concatenated_and_mutationstrength)

with open(annotated_vcf, 'r') as vcf:

    for line in vcf:

        if line.startswith('##'):
            continue
            
        elif line.startswith('#CHROM'):
            print(line)    
            
        else:

            line = line.strip()
            items = line.split('\t')
            contig = items[0] #the #CHROM column in the VCF
            # print(contig)
            position = items[1] #the POS column in the VCF

            concatenated = contig + "." + position #creates a string that contains both the chromosome and the position information
            if concatenated in concatenated_and_genome_position:
                print(line)
           
           
           
           