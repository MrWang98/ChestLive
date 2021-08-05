import numpy as np
import scipy.io as scio
# skiplist = [19, 20, 22, 23, 24, 26, 45, 59, 69]

def ReadData(data_All1):
    data_cup1_tmp = [] #HiSiri
    data_cup1_tmp0 = np.concatenate((data_All1[0]['cor1_ccc2'],
                                     np.zeros((2000 - data_All1[0]['cor1_ccc2'].shape[0], 39))))
    for i in range(1, 100):
        # if i in skiplist:
        #     continue
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

def LoadData():
    # path_base_cat1 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/TheFirstData/Guanqianyun/Matlab20200407/CUP'
    # path_base_cat2 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/TheFirstData/ZZA/MimicGQY/CUP'
    path_base_cat3 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/TheFirstData/XZQ/MimicGQY/CUP'
    path_base_cat4 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/TheFirstData/Lirouting/MimicGQY/CUP'
    path_base_cat5 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/TheFirstData/Xusixiao/MimicGQY/CUP'


    path_base_cat1 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/Attack/Matlab/OuRunMin/HeySiri'
    path_base_cat2 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/OuRunMin/Matlab/HeySiri'
    # path_base_cat2 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/Attackmimic/WithBreath/XueMeng/Matlab20200511/HiSiri'

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
    data_cup1 = ReadData(data_All1)
    # print("Person1 shape", data_cup1.shape, '\t', "Person1 data num", len(data_All1))

    # # Get the data from person two.
    for i in range(1, 101):
        path_tmp = path_base_cat2 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        if i == 1:
            data_All2 = [data_tmp]
            cupLabels.append('B')  # label person 2
        else:
            data_All2.append(data_tmp)
            cupLabels.append('B')  # label person 2
    data_cup2 = ReadData(data_All2)
    # # print("Person2 shape", data_cup2.shape, '\t', "Person2 data num", len(data_All2))

    # # Get the data from person third.
    # for i in range(1, 101):
    #     path_tmp = path_base_cat3 + str(i) + '.mat'
    #     data_tmp = scio.loadmat(path_tmp)
    #     if i == 1:
    #         data_All3 = [data_tmp]
    #         cupLabels.append('B')  # label person 3
    #     else:
    #         data_All3.append(data_tmp)
    #         cupLabels.append('B')  # label person 3
    # data_cup3 = ReadData(data_All3)
    # # print("Person3 shape", data_cup3.shape, '\t', "Person3 data num", len(data_All3))
    #
    # # Get the data from person four.
    # for i in range(1, 101):
    #     path_tmp = path_base_cat4 + str(i) + '.mat'
    #     data_tmp = scio.loadmat(path_tmp)
    #     if i == 1:
    #         data_All4 = [data_tmp]
    #         cupLabels.append('B')  # label person 3
    #     else:
    #         data_All4.append(data_tmp)
    #         cupLabels.append('B')  # label person 3
    # data_cup4 = ReadData(data_All4)
    # # print("Person4 shape", data_cup4.shape, '\t', "Person4 data num", len(data_All4))
    #
    # # Get the data from person five.
    # for i in range(1, 101):
    #     path_tmp = path_base_cat5 + str(i) + '.mat'
    #     data_tmp = scio.loadmat(path_tmp)
    #     if i == 1:
    #         data_All5 = [data_tmp]
    #         cupLabels.append('B')  # label person 3
    #     else:
    #         data_All5.append(data_tmp)
    #         cupLabels.append('B')  # label person 3
    # data_cup5 = ReadData(data_All5)
    # # print("Person5 shape", data_cup5.shape, '\t', "Person5 data num", len(data_All5))

    # Concatenate the data off three people.
    # cupLabels = np.delete(cupLabels, skiplist)
    # cupLabels = np.delete(cupLabels, [i + 91 for i in skiplist]) %, data_cup3, data_cup4, data_cup5
    data_cup = np.concatenate((data_cup1, data_cup2), axis=0)

    return data_cup, cupLabels
