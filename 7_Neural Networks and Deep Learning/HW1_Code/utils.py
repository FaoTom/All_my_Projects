import torch
from torch.utils.data import Dataset
import torch.optim as optim
import matplotlib.pyplot as plt
import numpy as np

def PlotLoss(loss, string_label):
    x = np.arange(1,len(loss)+1)
    if string_label == 'Train':
        plt.plot(x,loss, label = string_label, color = 'salmon')
    elif (string_label == 'Validation' or string_label == 'Test'):
        plt.semilogy(loss, label = string_label, color='firebrick')
    else: raise ValueError('[INVALID LABEL] Insert Train, Validation or Test as a label')

    title_string = string_label+' loss vs Epoch'

    plt.yscale('log')
    plt.xlabel('Epoch')
    plt.xlabel('Loss')
    plt.legend()
    plt.grid()
    plt.title(title_string)
    plt.tight_layout()
    plt.savefig(string_label+'_losses.pdf')

def PlotLosses(loss_vector, string_label):
    for i in range(len(loss_vector)):
        plt.plot(loss_vector[i])
    plt.xlabel('Epoch')
    plt.ylim(0,1)
    plt.savefig(string_label+'_losses.pdf')

class CsvDataset(Dataset):

    def __init__(self, df, transform=None):
        ''' (from the labs, slightly modified to have pd.DataFrame in input)

        -> This class represents a Python iterable over a dataset with
           support for custom operations on it.

        Parameters:
        df        - pd.DataFrame in input
        transform - optional transform to be applied on a sample

        '''
        self.transform = transform
        # Now self.data contains all our dataset.
        self.data = df
        self.names = np.array(df.columns)

    def __len__(self):
        # The length of the dataset is simply the length of the self.data list
        return len(self.data)

    def __getitem__(self, idx):
        # Our sample is the element idx of the list self.data
        sample = (self.data.iloc[idx][self.names[0]], self.data.iloc[idx][self.names[1]])
        if self.transform:
            sample = self.transform(sample)
        return sample

class ToTensor(object):

    def __call__(self, sample):
        ''' (from the labs)

        -> This class convert sample to torch.Tensor type

        Parameters:
        sample  - input sample to be converted

        '''
        x, y = sample
        return (torch.tensor([x]).float(),
                torch.tensor([y]).float())


