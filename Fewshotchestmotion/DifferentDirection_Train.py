import numpy as np
import scipy.io as scio

def ReadData1(data_All1):
    data_cup1_tmp = []
    data_cup1_tmp0 = np.concatenate((data_All1[0],
                                     np.zeros((2000 - data_All1[0].shape[0], 39))))
    Range1 = range(1, 70)
    for i in Range1:
        data_cup1_tmp = np.concatenate((data_All1[i],
                                        np.zeros((2000 - data_All1[i].shape[0], 39))))
        if i == 1:
            data_cup1_tmp_a = np.array(data_cup1_tmp0)
            data_cup1_tmp_b = np.array(data_cup1_tmp)
            data_cup11 = np.array([data_cup1_tmp_a, data_cup1_tmp_b])
        else:
            data_cup1_com_1 = np.append(data_cup11, data_cup1_tmp)
            data_cup1_dim = data_cup11.shape
            data_cup11 = data_cup1_com_1.reshape(data_cup1_dim[0] + 1,
                                                 data_cup1_dim[1],
                                                 data_cup1_dim[2])

    return data_cup11

def LoadData():

    path_base_cat1 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/XueMeng/Matlab/HeySiri'
    # KuangRuiLin XueMeng
    path_base_cat2 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/QinDang/Matlab/HeySiri'

    cupLabels = []

    Range1 = range(1,71)

    # Get the data from person one.
    for i in Range1:
        path_tmp = path_base_cat1 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)
        if i == 1:
            data_All1 = [together_data]
            cupLabels.append('A')  # label person 1
        else:
            data_All1.append(together_data)
            cupLabels.append('A')  # label person 1
    data_cup1 = ReadData1(data_All1)

    # Get the data from person two.
    for i in Range1:
        path_tmp = path_base_cat2 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)
        if i == 1:
            data_All2 = [together_data]
            cupLabels.append('B')  # label person 2
        else:
            data_All2.append(together_data)
            cupLabels.append('B')  # label person 2
    data_cup2 = ReadData1(data_All2)

    # Concatenate the data off three people.
    data_cup = np.concatenate((data_cup1, data_cup2), axis=0)

    return data_cup, cupLabels
