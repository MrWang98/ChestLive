import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pylab as pl
import scipy.io as scio

ChestTrainResultPath = 'ChestLiveTrainingSet.mat'
ChestTrainResult = scio.loadmat(ChestTrainResultPath)
print(ChestTrainResult.keys())
User1 = ChestTrainResult['U1']
User2 = ChestTrainResult['U2']
User3 = ChestTrainResult['U3']
User4 = ChestTrainResult['U4']
User5 = ChestTrainResult['U5']
Times = ChestTrainResult['times']
plt.plot(Times, User1) #, color="r", linewidth=3.0
plt.plot(Times, User2) #, color="m", linewidth=3.0
plt.plot(Times, User3) #, color="y", linewidth=3.0
plt.plot(Times, User4) #, color="g", linewidth=3.0
plt.plot(Times, User5) #, color="c", linewidth=3.0
plt.show()

