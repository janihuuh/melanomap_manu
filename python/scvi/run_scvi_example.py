## Init
import os
import numpy as np
import numpy.random as random
import pandas as pd

from scvi.dataset.dataset import GeneExpressionDataset
from scvi.dataset.csv import CsvDataset
from scvi.inference import UnsupervisedTrainer
from scvi.models import SCANVI, VAE
from scvi.inference.autotune import auto_tune_scvi_model

from umap import UMAP

import torch
import scanpy as sc
import louvain

import logging
import pickle
from hyperopt import hp


# %matplotlib inline

use_cuda = True
n_epochs_all = None
save_path = ''
show_plot = True
os.chdir("/scratch/cs/csb/projects/lgll/")


## Download samples
huuhtanen_01_00 = CsvDataset(filename='results/scvi/input_files/huuhtanen_01_00.csv', save_path='', sep=',', new_n_genes=False)
huuhtanen_02_00 = CsvDataset(filename='results/scvi/input_files/huuhtanen_02_00.csv', save_path='', sep=',', new_n_genes=False)
huuhtanen_03_00 = CsvDataset(filename='results/scvi/input_files/huuhtanen_03_00.csv', save_path='', sep=',', new_n_genes=False)
huuhtanen_04_00 = CsvDataset(filename='results/scvi/input_files/huuhtanen_04_00.csv', save_path='', sep=',', new_n_genes=False)
huuhtanen_05_00 = CsvDataset(filename='results/scvi/input_files/huuhtanen_05_00.csv', save_path='', sep=',', new_n_genes=False)
huuhtanen_06_00 = CsvDataset(filename='results/scvi/input_files/huuhtanen_06_00.csv', save_path='', sep=',', new_n_genes=False)

huuhtanen_01_01 = CsvDataset(filename='results/scvi/input_files/huuhtanen_01_01.csv', save_path='', sep=',', new_n_genes=False)
huuhtanen_02_01 = CsvDataset(filename='results/scvi/input_files/huuhtanen_02_01.csv', save_path='', sep=',', new_n_genes=False)
huuhtanen_03_01 = CsvDataset(filename='results/scvi/input_files/huuhtanen_03_01.csv', save_path='', sep=',', new_n_genes=False)
huuhtanen_04_01 = CsvDataset(filename='results/scvi/input_files/huuhtanen_04_01.csv', save_path='', sep=',', new_n_genes=False)
huuhtanen_05_01 = CsvDataset(filename='results/scvi/input_files/huuhtanen_05_01.csv', save_path='', sep=',', new_n_genes=False)
huuhtanen_06_01 = CsvDataset(filename='results/scvi/input_files/huuhtanen_06_01.csv', save_path='', sep=',', new_n_genes=False)

huuhtanen_01_03 = CsvDataset(filename='results/scvi/input_files/huuhtanen_01_03.csv', save_path='', sep=',', new_n_genes=False)
huuhtanen_02_03 = CsvDataset(filename='results/scvi/input_files/huuhtanen_02_03.csv', save_path='', sep=',', new_n_genes=False)
huuhtanen_03_03 = CsvDataset(filename='results/scvi/input_files/huuhtanen_03_03.csv', save_path='', sep=',', new_n_genes=False)
huuhtanen_04_03 = CsvDataset(filename='results/scvi/input_files/huuhtanen_04_03.csv', save_path='', sep=',', new_n_genes=False)
huuhtanen_05_03 = CsvDataset(filename='results/scvi/input_files/huuhtanen_05_03.csv', save_path='', sep=',', new_n_genes=False)
huuhtanen_06_03 = CsvDataset(filename='results/scvi/input_files/huuhtanen_06_03.csv', save_path='', sep=',', new_n_genes=False)

all_dataset = GeneExpressionDataset()
all_dataset.populate_from_per_batch_list(Xs = [ huuhtanen_01_00.X,
                                                huuhtanen_02_00.X,
                                                huuhtanen_03_00.X,
                                                huuhtanen_04_00.X,
                                                huuhtanen_05_00.X,
                                                huuhtanen_06_00.X,

                                                huuhtanen_01_01.X,
                                                huuhtanen_02_01.X,
                                                huuhtanen_03_01.X,
                                                huuhtanen_04_01.X,
                                                huuhtanen_05_01.X,
                                                huuhtanen_06_01.X,

                                                huuhtanen_01_03.X,
                                                huuhtanen_02_03.X,
                                                huuhtanen_03_03.X,
                                                huuhtanen_04_03.X,
                                                huuhtanen_05_03.X,
                                                huuhtanen_06_03.X,
                                                ])


## Train, save and fin
vae      = VAE(all_dataset.nb_genes, n_batch=all_dataset.n_batches, n_labels=all_dataset.n_labels, n_hidden=128, n_latent=30, n_layers=2, dispersion='gene')
trainer  = UnsupervisedTrainer(vae, all_dataset, train_size=1.0)
trainer.train(n_epochs=100)
torch.save(trainer.model.state_dict(), 'results/scvi/results/lgll_oneshot.pkl')

## Sample posterior to get latent representation and save those embeddings
full = trainer.create_posterior(trainer.model, all_dataset, indices=np.arange(len(all_dataset)))

latent, batch_indices, labels = full.sequential().get_latent()
batch_indices = batch_indices.ravel()

np.savetxt("results/scvi/results/melanomap_oneshot_latent.csv", latent, delimiter=",")
np.savetxt("results/scvi/results/melanomap_oneshot_indices.csv", batch_indices, delimiter=",")
