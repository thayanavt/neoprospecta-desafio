import curation_and_analysis_db as cad
import os

path_file = os.path.dirname(__file__)
email = "thayana.vt@gmail.com"

path_database = f"{path_file}/fas/banco.fasta" 
path_16S = f"{path_file}/fas/banco_16S.fasta" 
ids_and_seqs = cad.get_ids_and_seq(path_database)                                       #chamando a funcao get_ids_and_seq
ids_16S = cad.get_ids_16S(ids_and_seqs.keys(), email)                                   #chamando a funcao get_ids_16S
cad.write_db_only_16S(ids_and_seqs, ids_16S, path_16S)                                  #chamando a funcao write_db_only_16S
    
