import numpy as np
import scipy.io as scio

def ReadData1(data_All1):
    data_cup1_tmp = []
    data_cup1_tmp0 = np.concatenate((data_All1[0],
                                     np.zeros((2000 - data_All1[0].shape[0], 39))))
    for i in range(1, 100):
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
    path_base_cat1 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/GQY/Matlab/HeySiri'
    path_base_cat2 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/DengYangTao/Matlab/HeySiri'
    path_base_cat3 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/WuYuan/Matlab/HeySiri'
    path_base_cat4 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/HuHaiYan/Matlab/HeySiri'
    path_base_cat5 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/OuRunMin/Matlab/HeySiri'


    cupLabels = []

    # Get the data from person one.
    for i in range(1, 101):
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
    # print("Person1 shape", data_cup1.shape, '\t', "Person1 data num", len(data_All1))

    # Get the data from person two.
    for i in range(1, 101):
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
    # print("Person2 shape", data_cup2.shape, '\t', "Person2 data num", len(data_All2))

    # # Get the data from person Three.
    for i in range(1, 101):
        path_tmp = path_base_cat3 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)

        if i == 1:
            data_All3 = [together_data]
            cupLabels.append('C')  # label person 3
        else:
            data_All3.append(together_data)
            cupLabels.append('C')  # label person 3
    data_cup3 = ReadData1(data_All3)
    # print("Person3 shape", data_cup3.shape, '\t', "Person3 data num", len(data_All3))

    # Get the data from person Four.
    for i in range(1, 101):
        path_tmp = path_base_cat4 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)

        if i == 1:
            data_All4 = [together_data]
            cupLabels.append('D')  # label person 3
        else:
            data_All4.append(together_data)
            cupLabels.append('D')  # label person 3
    data_cup4 = ReadData1(data_All4)
    # print("Person4 shape", data_cup4.shape, '\t', "Person4 data num", len(data_All4))

    # Get the data from person Five.
    for i in range(1, 101):
        path_tmp = path_base_cat5 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)

        if i == 1:
            data_All5 = [together_data]
            cupLabels.append('E')  # label person 3
        else:
            data_All5.append(together_data)
            cupLabels.append('E')  # label person 3
    data_cup5 = ReadData1(data_All5)
    # print("Person5 shape", data_cup5.shape, '\t', "Person5 data num", len(data_All5))

    # Concatenate the data off three people.
    data_cup = np.concatenate((data_cup1, data_cup2, data_cup3, data_cup4, data_cup5), axis=0)  #

    return data_cup, cupLabels
