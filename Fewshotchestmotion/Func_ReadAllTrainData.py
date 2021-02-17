import numpy as np
import scipy.io as scio

def ReadData1(data_All1):
    data_cup1_tmp = []
    data_cup1_tmp0 = np.concatenate((data_All1[0],
                                     np.zeros((2000 - data_All1[0].shape[0], 39))))
    Range1 = range(1, 100)
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
    path_base_cat2 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/LPY/Matlab/HeySiri'
    path_base_cat3 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/ZhuXiaoTian/Matlab/HeySiri'
    path_base_cat4 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/HuangJiaQian/Matlab/HeySiri'
    path_base_cat5 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/LiLinWei/Matlab/HeySiri'

    path_base_cat6 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/GQY/Matlab/HeySiri'
    path_base_cat7 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/ZhuTianLin/Matlab/HeySiri'
    path_base_cat8 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/XZQ/Matlab/HeySiri'
    path_base_cat9 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/QinDang/Matlab/HeySiri'
    path_base_cat10 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/OuRunMin/Matlab/HeySiri'

    path_base_cat11 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/WuYuan/Matlab/HeySiri'
    path_base_cat12 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/BiHongLiang/Matlab/HeySiri'
    path_base_cat13 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/HuHaiYan/Matlab/HeySiri'
    path_base_cat14 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/DengYangTao/Matlab/HeySiri'
    path_base_cat15 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/KuangRuiLin/Matlab/HeySiri'

    cupLabels = []

    Range1 = range(1,101)

    # Get the data from person one.
    # 4 6 13 14 15. multi-verify.
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

    # # Get the data from person Three.
    for i in Range1:
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

    # Get the data from person Four.
    for i in Range1:
        path_tmp = path_base_cat4 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)

        if i == 1:
            data_All4 = [together_data]
            cupLabels.append('D')  # label person 4
        else:
            data_All4.append(together_data)
            cupLabels.append('D')  # label person 4
    data_cup4 = ReadData1(data_All4)

    # Get the data from person Five.
    for i in Range1:
        path_tmp = path_base_cat5 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)

        if i == 1:
            data_All5 = [together_data]
            cupLabels.append('E')  # label person 5
        else:
            data_All5.append(together_data)
            cupLabels.append('E')  # label person 5
    data_cup5 = ReadData1(data_All5)

    # Get the data from person one.
    for i in Range1:
        path_tmp = path_base_cat6 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)
        if i == 1:
            data_All6 = [together_data]
            cupLabels.append('F')  # label person 6
        else:
            data_All6.append(together_data)
            cupLabels.append('F')  # label person 6
    data_cup6 = ReadData1(data_All6)

    # Get the data from person two.
    for i in Range1:
        path_tmp = path_base_cat7 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)
        if i == 1:
            data_All7 = [together_data]
            cupLabels.append('G')  # label person 7
        else:
            data_All7.append(together_data)
            cupLabels.append('G')  # label person 7
    data_cup7 = ReadData1(data_All7)

    # # Get the data from person Three.
    for i in Range1:
        path_tmp = path_base_cat8 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)

        if i == 1:
            data_All8 = [together_data]
            cupLabels.append('H')  # label person 8
        else:
            data_All8.append(together_data)
            cupLabels.append('H')  # label person 8
    data_cup8 = ReadData1(data_All8)

    # Get the data from person Four.
    for i in Range1:
        path_tmp = path_base_cat9 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)

        if i == 1:
            data_All9 = [together_data]
            cupLabels.append('I')  # label person 9
        else:
            data_All9.append(together_data)
            cupLabels.append('I')  # label person 9
    data_cup9 = ReadData1(data_All9)

    # Get the data from person Five.
    for i in Range1:
        path_tmp = path_base_cat10 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)

        if i == 1:
            data_All10 = [together_data]
            cupLabels.append('J')  # label person 10
        else:
            data_All10.append(together_data)
            cupLabels.append('J')  # label person 10
    data_cup10 = ReadData1(data_All10)


    # Get the data from person one.
    for i in Range1:
        path_tmp = path_base_cat11 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)
        if i == 1:
            data_All11 = [together_data]
            cupLabels.append('K')  # label person 11
        else:
            data_All11.append(together_data)
            cupLabels.append('K')  # label person 11
    data_cup11 = ReadData1(data_All11)

    # Get the data from person two.
    for i in Range1:
        path_tmp = path_base_cat12 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)
        if i == 1:
            data_All12 = [together_data]
            cupLabels.append('L')  # label person 12
        else:
            data_All12.append(together_data)
            cupLabels.append('L')  # label person 12
    data_cup12 = ReadData1(data_All12)

    # # Get the data from person Three.
    for i in Range1:
        path_tmp = path_base_cat13 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)

        if i == 1:
            data_All13 = [together_data]
            cupLabels.append('M')  # label person 13
        else:
            data_All13.append(together_data)
            cupLabels.append('M')  # label person 13
    data_cup13 = ReadData1(data_All13)

    # Get the data from person Four.
    for i in Range1:
        path_tmp = path_base_cat14 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)

        if i == 1:
            data_All14 = [together_data]
            cupLabels.append('N')  # label person 14
        else:
            data_All14.append(together_data)
            cupLabels.append('N')  # label person 14
    data_cup14 = ReadData1(data_All14)

    # Get the data from person Five.
    for i in Range1:
        path_tmp = path_base_cat15 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)

        if i == 1:
            data_All15 = [together_data]
            cupLabels.append('O')  # label person 15
        else:
            data_All15.append(together_data)
            cupLabels.append('O')  # label person 15
    data_cup15 = ReadData1(data_All15)

    # Concatenate the data off three people.
    data_cup = np.concatenate((data_cup1, data_cup2, data_cup3, data_cup5, data_cup7,
                                 data_cup8, data_cup9, data_cup10, data_cup11,
                                data_cup6, data_cup4, data_cup13, data_cup14, data_cup15,
                                data_cup12,), axis=0)  #

    return data_cup, cupLabels
