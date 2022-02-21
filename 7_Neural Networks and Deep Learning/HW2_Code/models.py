
import numpy as np
import matplotlib.pyplot as plt

import torch
import torch.nn as nn
import torch.nn.functional as F
import torchmetrics

import pytorch_lightning as pl

#_______________________
#LIGHTNING AUTOENCODER
class LightningAutoEncoder(pl.LightningModule):
 
    def __init__(self, hyper):

        super().__init__()
                                            
        self.hyper = hyper
        
        self.encoder = nn.Sequential(
            nn.Conv2d(in_channels = 1,  out_channels = 8,  kernel_size = 3, stride = 2, padding = 1),   #first convolutional layer
            nn.ReLU(),                                                                              #activation
            nn.Conv2d(in_channels = 8,  out_channels = 16, kernel_size = 3, stride = 2, padding = 1),   #second convolutional layer
            nn.ReLU(),                                                                              #activation
            nn.Conv2d(in_channels=16, out_channels=32, kernel_size=3, stride=2, padding=0),             #third convolutional layer
            nn.ReLU(),                                                                              #activation
            nn.Flatten(start_dim=1),                                                                    #flatten layer
            nn.Linear(in_features=3*3*32, out_features=64),                                             #first linear layer
            nn.ReLU(),                                                                              #activation
            nn.Linear(in_features=64, out_features = self.hyper['encoded_space_dim'])                                   #second linear layer
        )

        self.decoder = nn.Sequential(
            nn.Linear(in_features = self.hyper['encoded_space_dim'], out_features = 64),
            nn.ReLU(),
            nn.Linear(in_features = 64, out_features = 288),
            nn.ReLU(),
            nn.Unflatten(dim=1, unflattened_size=(32, 3, 3)),
            nn.ConvTranspose2d(in_channels = 32, out_channels = 16, kernel_size = 3, stride = 2, output_padding = 0),
            nn.ReLU(),
            nn.ConvTranspose2d(in_channels = 16, out_channels = 8,  kernel_size = 3, stride = 2, padding = 1, output_padding = 1),
            nn.ReLU(),
            nn.ConvTranspose2d(in_channels = 8,  out_channels = 1,  kernel_size = 3, stride = 2, padding = 1, output_padding = 1),
            nn.Sigmoid()
        )

    def configure_loss(self, loss_func):
        self.loss_func = loss_func

    def configure_optimizers(self):

        opt      = self.hyper['opt']
        lr       = self.hyper['lr']
        reg      = self.hyper['reg']

        if opt == 'SGD':
            optimizer = torch.optim.SGD(self.parameters(),  lr=lr, momentum=0.9, weight_decay=reg)
        elif opt == 'Adam':
            optimizer = torch.optim.Adam(self.parameters(), lr=lr, weight_decay=reg)
        elif opt == 'Adadelta':
            optimizer = torch.optim.Adadelta(self.parameters(), lr=lr, weight_decay=reg)
        elif opt == 'Adagrad':
            optimizer = torch.optim.Adagrad(self.parameters(), lr=lr, weight_decay=reg)
        
        return optimizer
        
    def forward(self,  x):
        return self.encoder(x)

    def training_step(self, batch, batch_idx,  loss_name = 'train_loss'):
        x, _  = batch
        z     = self.encoder(x)
        x_hat = self.decoder(z)
        loss  = self.loss_func(x_hat, x) 
        self.log(loss_name, loss)
        return loss
    
    def validation_step(self, batch, batch_idx, loss_name = 'val_loss'):
        x, _  = batch
        z     = self.encoder(x)
        x_hat = self.decoder(z)
        loss  = self.loss_func(x_hat, x) 
        self.log(loss_name, loss, prog_bar = True) 
        return loss

    def test_step(self, batch, batch_idx):
        self.validation_step(batch, batch_idx, loss_name='test_loss')
#_______________________
#SUPERVISED MACHINE LEARNING CAE
class SupervisedCAE(pl.LightningModule):
 
    def __init__(self, hyper, pre_trained_CAE):

        super().__init__()
                                            
        self.hyper = hyper

        self.encoder = pre_trained_CAE.encoder

        self.fine_tuner =  nn.Sequential(nn.Linear(self.hyper['encoded_space_dim'], 64),
                                         nn.ReLU(True),
                                         nn.Linear(64, 10),
                                         nn.LogSoftmax()
                                         )
                                         
        self.accuracy = torchmetrics.Accuracy()

    def configure_loss(self, loss_func):
        self.loss_func = loss_func

    def configure_optimizers(self):

        opt      = self.hyper['opt']
        lr       = self.hyper['lr']
        reg      = self.hyper['reg']

        if opt == 'SGD':
            optimizer = torch.optim.SGD(self.parameters(),  lr=lr, momentum=0.9, weight_decay=reg)
        elif opt == 'Adam':
            optimizer = torch.optim.Adam(self.parameters(), lr=lr, weight_decay=reg)
        
        return optimizer
        
    def forward(self,  x):
        x = self.encoder(x)
        return self.fine_tuner(x)

    def training_step(self, batch, batch_idx,  loss_name = 'train_loss'):
        x, y = batch
        z = self.forward(x)
        loss = self.loss_func(z, y)
        self.log(loss_name, loss)
        return loss
    
    def validation_step(self, batch, batch_idx, loss_name = 'val_loss'):
        x, y = batch
        z = self.forward(x)
        loss = self.loss_func(z, y)
        self.log(loss_name, loss, prog_bar=True)
        return loss

    def test_step(self, batch, batch_idx):
        x, y = batch #now labels matter!
        z = self.forward(x)
        self.log('accuracy', self.accuracy(z, y), prog_bar=True)
        return self.accuracy(z, y)
# _______________________
# VARIATIONAL AUTOENCODER
class VariationalAutoEncoder(pl.LightningModule):
    
    def __init__(self, hyper):

        super().__init__()

        self.hyper = hyper

        self.encoder = nn.Sequential(
            nn.Conv2d(1, 8, kernel_size=3, padding=1, stride=2),
            nn.ReLU(),
            nn.Conv2d(8, 16, kernel_size=3, padding=1, stride=2),
            nn.ReLU(),
            nn.Conv2d(16, 32, kernel_size=3, padding=0, stride=2),
            nn.ReLU(),
            nn.Flatten(start_dim=1)
        )

        self.FCmu = nn.Sequential(
            nn.Linear(288, 64),
            nn.ReLU(),
            nn.Linear(64, hyper['encoded_space_dim'])
        )
        
        self.FCsi = nn.Sequential(
            nn.Linear(288, 64),
            nn.ReLU(),
            nn.Linear(64, hyper['encoded_space_dim'])
        )

        self.decoder =  nn.Sequential(
            nn.Linear(hyper['encoded_space_dim'], 64),
            nn.ReLU(True),
            nn.Linear(64, 288),
            nn.ReLU(True),
            nn.Unflatten(dim=1, unflattened_size=(32, 3, 3)),
            nn.ConvTranspose2d(32, 16, kernel_size=3, output_padding=0, stride=2),
            nn.ReLU(True),
            nn.ConvTranspose2d(16, 8, kernel_size=3, output_padding=1, padding=1, stride=2),
            nn.ReLU(True),
            nn.ConvTranspose2d(8, 1, kernel_size=3, output_padding=1, padding=1, stride=2),
            nn.Sigmoid()
        )

    def var_loss(self, x_hat, x_true, mu, log_var):
        reco_loss = F.mse_loss(x_hat, x_true, reduction='sum') 
        kl_div    = -0.5 * torch.sum(1. + log_var - mu**2 - torch.exp(log_var))
        return reco_loss + kl_div

    def sampling(self, mu, log_si):
        return mu + torch.randn_like(mu) * torch.exp(0.5 * log_si)

    def forward(self, x):
        
        internal_repr = self.encoder(x)

        pred_means  = self.FCmu(internal_repr)
        pred_log_si = self.FCsi(internal_repr) 
        sample = self.sampling(pred_means, pred_log_si)
        
        return sample, pred_means, pred_log_si

    def training_step(self, batch, batch_idx, log_name='train_loss', prog_bar=False):
        x, _ = batch
        sample, pred_means, pred_log_si = self.forward(x)
        reconstructed = self.decoder(sample)
        loss = self.var_loss(reconstructed, x, pred_means, pred_log_si)
        self.log(log_name, loss, prog_bar=prog_bar)
        
        return loss

    def configure_optimizers(self):
        opt      = self.hyper['opt']
        lr       = self.hyper['lr']
        reg      = self.hyper['reg']

        if opt == 'SGD':
            optimizer = torch.optim.SGD(self.parameters(),  lr=lr, momentum=0.9, weight_decay=reg)
        elif opt == 'Adam':
            optimizer = torch.optim.Adam(self.parameters(), lr=lr, weight_decay=reg)
        
        return optimizer

    def validation_step(self, batch, batch_idx, log_name='val_loss'):

        return self.training_step(batch, batch_idx, log_name=log_name, prog_bar=True)

    def test_step(self, batch, batch_idx, log_name='test_loss'):

        return self.training_step(batch, batch_idx, log_name=log_name, prog_bar=True)