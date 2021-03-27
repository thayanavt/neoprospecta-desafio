"""Modulo para realizacao de curadoria e analise de sequencias

Requer instalacao de biopython, matplotlib e numpy
"""

from Bio import SeqIO, Entrez
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

def get_ids_and_seq(path_database): 
    """Extrai os ids e sequencias de um arquivo fasta
    para um dicionario
    
    Parameters: 
    path_database           -- caminho para o banco de dados

    Returns: dict
    """
    ids_and_seqs = {}    
    for seq_record in SeqIO.parse(path_database, "fasta"): 
        ids_and_seqs[seq_record.id] = seq_record.seq

    #o dicionario criado possui: chaves = ids | elementos = sequencias
    return ids_and_seqs

def get_ids_16S(ids, email):
    """Acessa o Entrez (sistema de banco de dados de biologia molecular)
    e seleciona apenas os ids das sequências que seguirem as regras: 
        - Ser do gene 16S
        - Pertencer a um organismo procarioto
        
    Parameters: 
    ids                     -- um id ou mais 
    email                   -- email para acessar o Entrez
        
    Returns: list
    """
    Entrez.email = email
    ids_16S = []

    handle = Entrez.efetch(db='nucleotide', rettype='gb', retmode='text', id=ids)       #acessa o Entrez
    for seq_informations in SeqIO.parse(handle, 'gb'):                                  #analisa as informações obtidas no acesso ao Entrez
        for feature in seq_informations.features:        
            if feature.type == 'source':                                                #verifica se o gene tem como origem uma organela
                organelle = feature.qualifiers.get('organelle')
            if (feature.type=='rRNA' and                                                #verifica se o gene eh 16S e que não tem como origem uma organela
            '16S ribosomal RNA' in 
            feature.qualifiers.get('product') and 
            organelle==None):
                ids_16S.append(seq_informations.id)                                     #se sim, salva esse id na lista ids_16S
    return ids_16S 

def write_db_only_16S(ids_and_seqs, ids_16S, path_database_analyzed): 
    """Cria um arquivo fasta selecionando ids e sequências específicas

    Parameters: 
    ids_and_seqs            -- dicionário com os ids e sequências totais
    ids_16S                 -- lista com os ids correspondentes ao gene 16S
    path_database_analyzed  -- caminho para salvar o novo banco de dados 
                               (se quiser apenas substituir o já existente, 
                                coloque o mesmo caminho)
    """
    lines_db_analyzed = []
    for id in ids_and_seqs:
        if id in ids_16S:                                                              
            lines_db_analyzed.append(f'>{id}\n{ids_and_seqs[id]}\n')                   #se sim, adiciona um id e sequencia no formato fasta na lista lines_db_analyzed
    
    with open(path_database_analyzed, 'w') as db_16S: 
        db_16S.writelines(lines_db_analyzed)                                           #escreve o arquivo com os itens da lista lines_db_analyzed

def histogram_length_seqs(lenght_seqs, path_histogram, color): 
    """Cria um histograma com a distribuição do comprimento das sequencias
    
    Parameters: 
    lenght_seqs             -- lista com os comprimentos das sequências
    path_histogram          -- caminho para salvar o histograma
    color                   -- cor do histograma
    """
    plt.hist(lenght_seqs, rwidth=0.9, color=color)   
    plt.title('Comprimento das sequências')
    plt.xlabel('Quantidade de bases')
    plt.ylabel('Frequência')
    plt.savefig(path_histogram)

def get_taxonomy(ids, email): 
    """Extrai o gênero e a espécie de um organismo, a partir do id, 
    acessando o Entrez (sistema de banco de dados de biologia molecular).
    
    Parameters: 
    ids                     -- um id ou mais
    email                   -- email para acessar o Entrez
    """
    Entrez.email = email
    genus_species = {}

    handle = Entrez.efetch(db='nucleotide', rettype='gb', retmode='text', id=ids)       #acessa o Entrez    
    for seq_informations in SeqIO.parse(handle, 'gb'):                                  #analisa as informacoes obtidas no acesso ao Entrez
        taxonomy_list = seq_informations.annotations['taxonomy'][-2:]                   #seleciona apenas o genero e a especie
        taxonomy_string = ' '.join([str(taxonomy) for taxonomy in taxonomy_list])       #converte a lista taxonomy_list em string

        if taxonomy_string not in genus_species:                                        
            genus_species[taxonomy_string] = 1

        else:
            genus_species[taxonomy_string] = genus_species.get(taxonomy_string) + 1     #seleciona os valores da chave em questao e soma 1

    #o dicionario criado possui: chaves = genero especie | elementos = frequencia 
    return genus_species 

def pie_more_abundant(taxonomy, path_pie, min_value, colors): 
    """Cria um gráfico de pizza com a distribuição.
    
    Parameters: 
    taxonomy           -- dicionario com taxonomia na chave e frequencia no elemento
    path_pie           -- caminho para salvar o grafico de pizza
    min_value          -- o menor valor de frequencia valido
    colors             -- cores para cada id
    """
    labels = []
    sizes = []
    
    for name in taxonomy:                                                                
        if taxonomy[name] >= min_value:
            labels.append(name)                                                         #adiciona a chave do dicionario taxonomy na lista labels
            sizes.append(taxonomy[name])                                                #adiciona o elemento do dicionario taxonomy na lista sizes
    
    plt.figure(figsize = (35, 18))
    legend = []

    for color, name in zip(colors, labels):                                             #correlaciona os nomes com as cores
        legend.append(mpatches.Patch(color = color, label = name))                      #salva a correlação na lista legend
    
    plt.legend(handles=legend, bbox_to_anchor=(1.05, 1), fontsize='xx-large')
    plt.pie(sizes, shadow=True, autopct='%1.1f%%', textprops={'fontsize': 17},
            colors = colors, pctdistance=0.9, radius=1.2)
    plt.title('Distribuição das taxonomias mais frequentes', fontsize=30)
    plt.savefig(path_pie)

