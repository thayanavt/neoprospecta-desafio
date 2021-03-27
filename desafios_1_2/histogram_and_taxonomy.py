import curation_and_analysis_db as cad
import os

path_file = os.path.dirname(__file__)
email = "thayana.vt@gmail.com"

path_database = f"{path_file}/fas/banco_16S.fasta"
ids_and_seqs = cad.get_ids_and_seq(path_database)                                           #chamando a funcao get_ids_and_seq

### Histogram ###
path_histogram = f"{path_file}/histogram.png"
lenght_seqs = [len(seq) for seq in ids_and_seqs.values()]                                   #cria uma lista apenas com os elementos do dicionario
cad.histogram_length_seqs(lenght_seqs, path_histogram, 'turquoise')                         #chamando a funcao histogram_length_seqs                    

### Pie ###
path_pie = f"{path_file}/pie.png"
genus_species = cad.get_taxonomy(ids_and_seqs.keys(), email)                        
porcent_1 = int(sum(genus_species.values()) / 100)                                          #calculando o valor 1%
colors = ['rosybrown', 'firebrick', 'darksalmon', 'lightcoral', 'palevioletred',
        'mediumpurple', 'plum', 'coral', 'violet', 'indianred','sienna', 'sandybrown', 
        'tan', 'darkkhaki', 'lightgreen', 'olivedrab', 'olive', 'palegreen', 'darkgreen', 
        'seagreen', 'mediumseagreen', 'lightseagreen', 'paleturquoise','darkcyan', 
        'slategray', 'darkorange', 'navajowhite', 'darkgoldenrod',  'cadetblue', 'skyblue']
cad.pie_more_abundant(genus_species, path_pie, porcent_1, colors)                           #chamando a funcao pie_more_abundant
