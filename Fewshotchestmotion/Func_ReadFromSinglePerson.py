import numpy as np
import scipy.io as scio

def ReadData1(data_All1):
    data_cup1_tmp = []
    data_cup1_tmp0 = np.concatenate((data_All1[0]['cor1_ccc2'],
                                     np.zeros((2000 - data_All1[0]['cor1_ccc2'].shape[0], 39))))
    for i in range(1, 100):
        data_cup1_tmp = np.concatenate((data_All1[i]['cor1_ccc2'],
                                        np.zeros((2000 - data_All1[i]['cor1_ccc2'].shape[0], 39))))
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

def ReadData2(data_All1):
    data_cup1_tmp = []
    data_cup1_tmp0 = np.concatenate((data_All1[0]['ccc2xue'],
                                     np.zeros((2000 - data_All1[0]['ccc2xue'].shape[0], 39))))
    for i in range(1, 100):
        data_cup1_tmp = np.concatenate((data_All1[i]['ccc2xue'],
                                        np.zeros((2000 - data_All1[i]['ccc2xue'].shape[0], 39))))
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
    path_base_cat1 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/TheFirstData/ZZA/CompareTohimself4161/CUP'
    path_base_cat2 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/TheFirstData/ZZA/CompareTohimself4162/CUP'

    cupLabels = []

    # Get the data from person one.
    for i in range(1, 101):
        path_tmp = path_base_cat1 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        if i == 1:
            data_All1 = [data_tmp]
            cupLabels.append('A')  # label person 1
        else:
            data_All1.append(data_tmp)
            cupLabels.append('A')  # label person 1
    data_cup1 = ReadData1(data_All1)
    # print("Person1 shape", data_cup1.shape, '\t', "Person1 data num", len(data_All1))

    # Get the data from person two.
    DataStr = 'ccc2xue'
    for i in range(1, 101):
        path_tmp = path_base_cat2 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        if i == 1:
            data_All2 = [data_tmp]
            cupLabels.append('B')  # label person 2
        else:
            data_All2.append(data_tmp)
            cupLabels.append('B')  # label person 2
    data_cup2 = ReadData2(data_All2)
    # print("Person2 shape", data_cup2.shape, '\t', "Person2 data num", len(data_All2))

    # Concatenate the data off three people.
    data_cup = np.concatenate((data_cup1, data_cup2), axis=0)

    return data_cup, cupLabels
