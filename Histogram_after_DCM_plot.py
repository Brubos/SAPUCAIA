#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 14:02:54 2023

@author: bruno.souza
"""
from optlnls.math import get_fwhm
import numpy as np
import matplotlib.pyplot as plt

# Nome do arquivo a ser lido
nome_arquivo = "Histogram_Intensity(m).txt"

# Leitura do arquivo e criação dos vetores
dados = np.loadtxt(nome_arquivo)
energia = dados[:, 0]
intensidade = dados[:, 1]

# Normalização da intensidade
max_intensidade = np.max(intensidade)
intensidade_normalizada = intensidade / max_intensidade

# Encontrando o valor médio de energia
energia_media = np.mean(energia)

# Largura à meia altura (FWHM)
valores = get_fwhm(energia, intensidade_normalizada)

# Criação do plot
plt.figure(facecolor='white')
plt.plot(energia, intensidade_normalizada, label="Dados simulados",color='#95cf52',alpha=0.5)
plt.xlabel("Energia [eV]")
plt.ylabel("Intensidade Normalizada")
plt.title("Histograma de Energia - Espelho")
plt.grid(True, alpha=0.25)  # Grid branco

# Adicionar linhas verticais
plt.axvline(x=valores[1], linestyle='--', alpha=0.4,linewidth=0.5) 
plt.axvline(x=valores[2], linestyle='--', alpha=0.4,linewidth=0.5)
plt.axvline(x=energia_media, linestyle='--', linewidth=1, alpha=0.5,color='gray')

offset = 0.08  # muda o tamanho da seta
# Adicionar linha horizontal verde com setas nas extremidades usando a cor #95cf52
plt.arrow(valores[1], valores[3], (valores[2] - valores[1]) - offset, 0, head_width=0.02, head_length=0.05, linestyle='-', linewidth=0.5,color='black')
plt.arrow(valores[2], valores[3], (valores[1] - valores[2]) + offset, 0, head_width=0.02, head_length=0.05, linestyle='-', linewidth=0.5, label=f'FWHM = {round(valores[0], 2)} eV',color='black')

# Exibir legenda
plt.legend()

# Salvar a figura com resolução alta (por exemplo, 300 DPI)
plt.savefig('Histogram.png',dpi=1000)   