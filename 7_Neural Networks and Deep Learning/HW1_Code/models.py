import torch.nn as nn
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import torch.nn.functional as F
import torch

from torch.autograd import Variable
from torch.nn import Linear, ReLU, CrossEntropyLoss, Sequential, Conv2d, MaxPool2d, Module, Softmax, BatchNorm2d, Dropout
from torch.optim import Adam, SGD

#_______________________
#FCN Simple
class FullyConnectedNeuralNetSimple(nn.Module):
    
    def __init__(self, net_par):
        '''

        -> This function initializes the network layers and activation functions

        Parameters:
        - Inside dictionary net_par:
            Ninp      - input size
            Nh1       - neurons in the hidden layer #1
            Nh2       - neurons in the hidden layer #2
            Nout      - output size
            drop_prob - probability of dropout

        '''
        super(FullyConnectedNeuralNetSimple, self).__init__()
        
        self.h1   = nn.Linear(net_par['Ninp'], net_par['Nh1'])  #first layer
        self.h2   = nn.Linear(net_par['Nh1'],  net_par['Nh2'])  #second layer
        self.out  = nn.Linear(net_par['Nh2'],  net_par['Nout']) #output layer
        self.act  = net_par['act_func']                         #activation function
    
    def forward(self, x):
        '''

        -> This function defines the forward pass of the network

        Parameters:
        x    - input of the network 
        
        '''
        x = self.h1(x)                                          #first layer
        x = self.act(x)                                         #first activation     
        x = self.h2(x)                                          #second layer
        x = self.act(x)                                         #second activation
        x = self.out(x)                                         #output layer

        return x

    def fit(self, train_dl, test_dl, hyp_par):
        '''

        -> This function fits the model to the data by doing
            backpropagation on the model weights with a custom loss
            function and a custom optimizer.

        Parameters:
        train_dl    - train DataLoader
        test_dl     - test DataLoader
        hyp_par     - other simulation hyperparameters:
                            num_epochs - number of epochs
                            device     - CPU or GPU (if possible)
                            opt        - optimizer
                            loss_fn    - loss function
        
        '''
        
        train_loss_log=[]
        val_loss_log=[]

        device=hyp_par['device']
        
        for epoch in tqdm(range(hyp_par['num_epochs'])):    #loop over the epochs
            
            self.train()                                    #TRAINING MODE
            epoch_tr_loss = []                              #training loss for each epoch
            for sample_batched in train_dl:                 #loop over the batches
                x_batch = sample_batched[0].to(device)      #batch separation: data
                y_batch = sample_batched[1].to(device)      #batch separation: labels
                out = self.forward(x_batch)                 #perform forward propagation
                tr_loss = hyp_par['loss_fn'](out, y_batch)  #compute loss function with prediction and label
                self.zero_grad()                            #reset accumulated gradients from previous batches 
                tr_loss.backward()                          #backpropagation (just compute derivatives given loss_fn)
                hyp_par['opt'].step()                       #update weights
                loss_batch = tr_loss.detach().cpu().numpy() #save train loss for this batch
                epoch_tr_loss.append(tr_loss.detach().cpu().numpy())
            train_loss_log.append(np.mean(epoch_tr_loss))   #append mean value of loss for the batches on each epoch

            self.eval()                                     #EVALUATION mode
            epoch_va_loss = []                              #validation loss for each epoch
            with torch.no_grad():                           #disable gradient tracking
                for sample_batched in test_dl:              #same as above
                    x_batch = sample_batched[0].to(device)
                    y_batch = sample_batched[1].to(device)
                    out = self.forward(x_batch)       
                    va_loss = hyp_par['loss_fn'](out, y_batch)     
                    loss_batch = va_loss.detach().cpu().numpy() 
                    epoch_va_loss.append(loss_batch)
            val_loss_log.append(np.mean(epoch_va_loss))    

        return train_loss_log, val_loss_log                 #return final losses for the simulation

    def reset(self):
        '''
        -> This function resets the parameters of each layer
           of self (by checking beforehand it is possible)
        
        '''
        for layer in self.children():
            if hasattr(layer, 'reset_parameters'):
                layer.reset_parameters()
#_______________________
#FCN Advanced
class FullyConnectedNeuralNet(nn.Module):
    
    def __init__(self, net_par):
        '''

        -> This function initializes the network layers and activation functions

        Parameters:
        - Inside dictionary net_par:
            Ninp      - input size
            Nh1       - neurons in the hidden layer #1
            Nh2       - neurons in the hidden layer #2
            Nh3       - neurons in the hidden layer #3
            Nout      - output size
            drop_prob - probability of dropout

        '''
        super(FullyConnectedNeuralNet, self).__init__()
        
        self.h1   = nn.Linear(net_par['Ninp'],  net_par['Nh1']) #1st hidden layer
        self.h2   = nn.Linear(net_par['Nh1'], net_par['Nh2'])   #2nd hidden layer
        self.h3   = nn.Linear(net_par['Nh2'], net_par['Nh3'])   #2nd hidden layer
        self.out  = nn.Linear(net_par['Nh3'], net_par['Nout'])  #output layer
        self.act  = net_par['act_func']                         #activation function
    
    def forward(self, x):
        '''

        -> This function defines the forward pass of the network

        Parameters:
        x    - input of the network 
        
        '''
        x = self.h1(x)                                          #first layer
        x = self.act(x)                                         #first activation     
        x = self.h2(x)                                          #secondl layer
        x = self.act(x)
        x = self.h3(x)
        x = self.act(x)                                         #second activation
        x = self.out(x)                                         #last layer

        return x

    def fit(self, train_dl, test_dl, hyp_par):
        '''

        -> This function fits the model to the data by doing
            backpropagation on the model weights with a custom loss
            function and a custom optimizer.

        Parameters:
        train_dl    - train DataLoader
        test_dl     - test DataLoader
        hyp_par     - other simulation hyperparameters:
                            num_epochs - number of epochs
                            device     - CPU or GPU (if possible)
                            opt        - optimizer
                            loss_fn    - loss function
        
        '''
        train_loss_log=[]
        val_loss_log=[]

        device=hyp_par['device']
        
        for epoch in tqdm(range(hyp_par['num_epochs'])):    #loop over the epochs
            
            self.train()                                    #TRAINING MODE
            epoch_tr_loss = []                              #training loss for each epoch
            for sample_batched in train_dl:                 #loop over the batches
                x_batch = sample_batched[0].to(device)      #batch separation: data
                y_batch = sample_batched[1].to(device)      #batch separation: labels
                out = self.forward(x_batch)                 #perform forward propagation
                tr_loss = hyp_par['loss_fn'](out, y_batch)  #compute loss function with prediction and label
                self.zero_grad()                            #reset accumulated gradients from previous batches 
                tr_loss.backward()                          #backpropagation (just compute derivatives given loss_fn)
                hyp_par['opt'].step()                       #update weights
                loss_batch = tr_loss.detach().cpu().numpy() #save train loss for this batch
                epoch_tr_loss.append(tr_loss.detach().cpu().numpy())
            train_loss_log.append(np.mean(epoch_tr_loss))   #append mean value of loss for the batches on each epoch

            self.eval()                                     #EVALUATION mode
            epoch_va_loss = []                              #validation loss for each epoch
            with torch.no_grad():                           #disable gradient tracking
                for sample_batched in test_dl:              #same as above
                    x_batch = sample_batched[0].to(device)
                    y_batch = sample_batched[1].to(device)
                    out = self.forward(x_batch)       
                    va_loss = hyp_par['loss_fn'](out, y_batch)     
                    loss_batch = va_loss.detach().cpu().numpy() 
                    epoch_va_loss.append(loss_batch)
            val_loss_log.append(np.mean(epoch_va_loss))    

        return train_loss_log, val_loss_log                 #return final losses for the simulation

    def reset(self):
        '''
        -> This function resets the parameters of each layer
           of self (by checking beforehand it is possible)
        
        '''
        for layer in self.children():
            if hasattr(layer, 'reset_parameters'):
                layer.reset_parameters()
#_______________________
#CCN     
class ConvolutionalNet(nn.Module):
    
    def __init__(self, prob_drop, prob_drop2, act_func, xavier):
        '''

            -> This function initializes the CNN architecture

            Parameters:
            - Inside dictionary net_par:
                channels
                prob_drop2D
                prob_drop1D

        '''
        verbose = False
        
        super().__init__()

        self.conv1 = nn.Conv2d(in_channels = 1,  out_channels=16, kernel_size=5, stride=1, padding=2)
        self.conv2 = nn.Conv2d(in_channels = 16, out_channels=8, kernel_size=5, stride=1, padding=2)
        self.mpool2d = nn.MaxPool2d(kernel_size = 2)
        self.actf = act_func
        self.dr  = nn.Dropout2d(p = prob_drop)
        self.dr2  = nn.Dropout2d(p = prob_drop2)
        self.fc1   = nn.Linear(in_features  = 8 * 7 * 7, out_features = 100)
        self.out   = nn.Linear(in_features  = 100, out_features = 10)
        
        if xavier:
            self.conv1.weight = nn.init.xavier_uniform_(self.conv1.weight, gain=nn.init.calculate_gain('relu'))
            self.conv2.weight = nn.init.xavier_uniform_(self.conv2.weight, gain=nn.init.calculate_gain('relu'))
            self.fc1.weight = nn.init.xavier_uniform_(self.fc1.weight, gain=nn.init.calculate_gain('relu'))

    def forward(self, x):

        x = self.conv1(x)
        x = self.actf(x)
        x = self.mpool2d(x)
        x = self.dr2(x)
        x = self.conv2(x)   
        x = self.actf(x)
        x = self.mpool2d(x)
        x = self.dr2(x)
        x = x.view(x.size(0),-1)
        x = self.fc1(x)
        x = self.actf(x)
        x = self.dr(x)
        x = self.out(x)
       
        return x
