import numpy as np
import tensorflow as tf
import tensorboard
import torch
from torch.utils.tensorboard import SummaryWriter

class MLPNN(torch.nn.Module):

    def __init__(self, nFeatures):

        super().__init__()
        
        
        # Set MLP parameters
        self.nFeatures = nFeatures
        self.batchSize = 1
        self.learningRate = 0.001
        self.nEpochs = 500
        
        # MLP Architecture
        self.fc = torch.nn.Sequential(
            torch.nn.Linear(self.nFeatures, 100),
            torch.nn.ReLU(),
            torch.nn.Linear(100, 20),
            torch.nn.ReLU(),
            torch.nn.Linear(20, 1), # In the last layer we are estimating a probability, so need output of 1
            torch.nn.Sigmoid())

        
    def forward(self, features):
        '''
        Forward operation
        '''
        return self.fc(features)
    

    def trainModel(self, trainingData: list):
        '''
        Train model
        '''
        
        # Extract features and labels
        trainingFeatures, trainingLabels = trainingData
        
        # Training configuration
        nTrainingSamples = trainingFeatures.shape[0]
        nTrainingBatches = nTrainingSamples // self.batchSize
        if nTrainingBatches * self.batchSize < nTrainingSamples:
            nTrainingBatches += 1

        # Moving to CPU
        device = torch.device('cpu')
        model = self
        model.to(device=device)

        # Optimization
        loss = torch.nn.BCELoss() # Loss function
        optimizer = torch.optim.SGD(model.parameters(), lr=self.learningRate, momentum=0.5) # Optimization approach

        # Training
        writer = SummaryWriter()  # Initialize logger
        for epoch in range(self.nEpochs):
            
            # -> training mode
            model.train()

            epochLoss = 0
            epochAccuracy = 0
            for batch in range(nTrainingBatches):
                
                # Zeroing gradients
                #optimizer.zero_grad()

                # Numpy -> tensor data conversion
                x = torch.tensor(
                    trainingFeatures[batch*self.batchSize:(batch+1)*self.batchSize, :],
                    device=device,
                    dtype=torch.float32)

                y = torch.tensor(
                    trainingLabels[batch*self.batchSize:(batch+1)*self.batchSize].reshape((-1, 1)),
                    device=device,
                    dtype=torch.float32)

                # Forward
                y_pred = model(x)
                batchLoss = loss(y_pred, y)  # Loss computation
                batchLoss.backward()         # Backpropagate
                optimizer.step()             # Update parameters
                optimizer.zero_grad()
                epochLoss += batchLoss.data * x.shape[0] / nTrainingSamples

                labels_pred = torch.round(y_pred)
                correct = (y == labels_pred).float()
                batchAccuracy = correct.sum() / correct.numel()
                epochAccuracy += batchAccuracy * 100 * x.shape[0] / nTrainingSamples
                
                # Write to tensorboard
                #writer.add_scalar('Batch Loss', batchLoss, batch)
                writer.add_scalar('Batch Accuracy', batchAccuracy, batch)

            model.train(False)
            
            # Validation after each epoch
            # Numpy -> tensor data conversion for validation set
            #x = torch.tensor(
            #    validationFeatures,
            #    device=device,
            #    dtype=torch.float32)

            #y = torch.tensor(
            #    validationLabels.reshape((-1, 1)),
            #    device=device,
            #    dtype=torch.float32)

            # Forward
            #y_pred = model(x)
            #labels_pred = torch.round(y_pred)
            #correct = (y == labels_pred).float()
            #validationAccuracy = correct.sum() / correct.numel()

            # Write to tensorboard
            writer.add_scalar('Epoch Loss', epochLoss, epoch)
            writer.add_scalar('Epoch Accuracy', epochAccuracy, epoch)
            #writer.add_scalar('Validation Accuracy', validationAccuracy, epoch)

        # Stop logger
        writer.close()

    def save(self, path):
        torch.save(self, path)

    def predict(self, x, proba=False):
        device = torch.device('cpu')
        model = self
        model.to(device=device)

        # Numpy -> tensor data conversion for validation set
        x = torch.tensor(x,
            device=device,
            dtype=torch.float32)

        # Get predictions
        y_predProb = model(x)
        y_predClass = torch.round(y_predProb)
        
        # Return predictions
        return y_predProb if proba else y_predClass

