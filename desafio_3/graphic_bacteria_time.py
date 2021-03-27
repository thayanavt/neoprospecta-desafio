"""Script para plotar um gráfico de barras, mostrando a contagem absoluta de 50 bacterias 
mais abundantes,, agrupando-as por tempo. 

Requer instalacao de pandas, matplotlib e numpy
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import os

path = os.path.dirname(__file__)

def insert_sum_column(df):
    """Insere uma coluna com a soma de todas as linhas
    
    Parameters: 
    df              -- dataframe 
    """
    df['SUM'] = df.sum(axis=1)

#salvando os arquivos iniciais como dataframes
df_bacterias = pd.read_csv(f"/{path}/tables/otu_table_tax_amostras.tsv", sep='\t')
df_group_time = pd.read_csv(f"/{path}/tables/metadata.csv", sep='\t')

df_group_time = df_group_time.groupby('time')['#group'].apply(list)                         #agrupa de acordo com a coluna time

insert_sum_column(df_bacterias)                                                             #chama a funcao insert_sum_column 

df_50bacterias = df_bacterias.nlargest(50, 'SUM', keep='first')                             #seleciona as 50 bacterias mais abundantes
df_50bacterias = df_50bacterias.rename(columns = lambda i : str(i)[:-5])                    #remove os 5 ultimos caracteres das colunas do df_50bacterias

df_early = pd.DataFrame(df_50bacterias[df_group_time['Early']])                             #cria um dataframe apenas com grupos relacionados ao Early
insert_sum_column(df_early)                                                                 #chama a funcao insert_sum_column
df_late = pd.DataFrame(df_50bacterias[df_group_time['Late']])                               #cria um dataframe  apenas com grupos relacionados ao Late
insert_sum_column(df_late)                                                                  #chama a funcao insert_sum_column

sum_early_late = {}

for counter, bacteria_taxonomy in enumerate(df_50bacterias.iloc[:,0]):                                
    sum_early_late[bacteria_taxonomy] = [df_early.iloc[counter, -1], df_late.iloc[counter, -1]] 

#o dicionario criado possui: chaves = taxonomia | elementos = soma de early e soma de late
    
colors = cm.gist_ncar(np.linspace(0, 1, len(sum_early_late)))                               #selecionando cores da paleta gist_ncar
legend = []

for color, bacteria_taxonomy in zip(colors, sum_early_late.keys()):                         #correlacionando a taxonomia com a cor
    legend.append(mpatches.Patch(color = color, label = bacteria_taxonomy))                 #salvando a correlação na lista legend

plt.figure(figsize=(22, 25))
bottom = [0,0]

for bacteria_taxonomy, color in zip(sum_early_late, colors):    
    plt.bar(['Early', 'Late'], sum_early_late[bacteria_taxonomy],                           #criando o grafico de barras
            color=color, bottom=bottom) 
    bottom = list(map(sum, zip(bottom, sum_early_late[bacteria_taxonomy])))                 #aumentando o degrau na barra para o proximo loop
    
plt.legend(handles=legend[::-1], loc=0, bbox_to_anchor=(1.05, 1), fontsize='xx-large')      #legenda eh colocada ao contrario para as cores se igualarem com a barra 
plt.yticks(fontsize=40)
plt.xticks(fontsize=40)
plt.title('50 bactérias mais abundantes agrupadas por tempo de desmame', fontsize=35)
plt.savefig(f'/{path}/bacterias.png', bbox_inches='tight')

