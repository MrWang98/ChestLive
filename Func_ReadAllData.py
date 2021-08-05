import numpy as np
import scipy.io as scio

def ReadData1(data_All1):
    data_cup1_tmp = []
    data_cup1_tmp0 = np.concatenate((data_All1[0],
                                     np.zeros((2000 - data_All1[0].shape[0], 39))))
    nda = np.array(data_All1).shape[0]
    Range1 = range(1, nda)
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

def LoadData(training=False):
    path_base_cat1 = 'Final/XueMeng/Matlab/HeySiri'
    path_base_cat2 = 'Final/LPY/Matlab/HeySiri'
    path_base_cat3 = 'Final/ZhuXiaoTian/Matlab/HeySiri'
    path_base_cat4 = 'Final/HuangJiaQian/Matlab/HeySiri'
    path_base_cat5 = 'Final/LiLinWei/Matlab/HeySiri'

    path_base_cat6 = 'Final/GQY/Matlab/HeySiri'
    path_base_cat7 = 'Final/ZhuTianLin/Matlab/HeySiri'
    path_base_cat8 = 'Final/XZQ/Matlab/HeySiri'
    path_base_cat9 = 'Final/QinDang/Matlab/HeySiri'
    path_base_cat10 = 'Final/OuRunMin/Matlab/HeySiri'

    path_base_cat11 = 'Final/WuYuan/Matlab/HeySiri'
    path_base_cat12 = 'Final/BiHongLiang/Matlab/HeySiri'
    path_base_cat13 = 'Final/HuHaiYan/Matlab/HeySiri'
    path_base_cat14 = 'Final/DengYangTao/Matlab/HeySiri'
    path_base_cat15 = 'Final/KuangRuiLin/Matlab/HeySiri'
    cupLabels = []

    Range2 = range(1, 101) #71
    # Used to divide data in a 7:3 ratio
    # DataSet_Split = "train" if training else "test"
    # if DataSet_Split == "train":
    #     Range2 = range(1, 71)
    # else:
    #     Range2 = range(71, 101)

    data_All1 = []
    data_All2 = []
    data_All3 = []
    data_All4 = []
    data_All5 = []

    data_All6 = []
    data_All7 = []
    data_All8 = []
    data_All9 = []
    data_All10 = []

    data_All11 = []
    data_All12 = []
    data_All13 = []
    data_All14 = []
    data_All15 = []

    # Used for labels
    label_All1 = []
    label_All2 = []
    label_All3 = []
    label_All4 = []
    label_All5 = []

    label_All6 = []
    label_All7 = []
    label_All8 = []
    label_All9 = []
    label_All10 = []

    label_All11 = []
    label_All12 = []
    label_All13 = []
    label_All14 = []
    label_All15 = []



    # Get the data from person one.
    for i in Range2:
        path_tmp = path_base_cat1 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)
        if i == 1:
            data_All1 = [together_data]
            label_All1.append('A')  # label person 1
        else:
            data_All1.append(together_data)
            label_All1.append('A')  # label person 1
    data_cup1 = ReadData1(data_All1)

    # Get the data from person two.
    # for i in Range2:
    #     path_tmp = path_base_cat2 + str(i) + '.mat'
    #     data_tmp = scio.loadmat(path_tmp)
    #     together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)
    #     if i == 1:
    #         data_All2 = [together_data]
    #         label_All2.append('B')  # label person 2
    #     else:
    #         data_All2.append(together_data)
    #         label_All2.append('B')  # label person 2
    # data_cup2 = ReadData1(data_All2)

    # # Get the data from person Three.
    for i in Range2:
        path_tmp = path_base_cat3 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)

        if i == 1:
            data_All3 = [together_data]
            label_All3.append('C')  # label person 3
        else:
            data_All3.append(together_data)
            label_All3.append('C')  # label person 3
    data_cup3 = ReadData1(data_All3)

    # Get the data from person Four.
    for i in Range2:
        path_tmp = path_base_cat4 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)

        if i == 1:
            data_All4 = [together_data]
            label_All4.append('D')  # label person 4
        else:
            data_All4.append(together_data)
            label_All4.append('D')  # label person 4
    data_cup4 = ReadData1(data_All4)

    # Get the data from person Five.
    for i in Range2:
        path_tmp = path_base_cat5 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)

        if i == 1:
            data_All5 = [together_data]
            label_All5.append('E')  # label person 5
        else:
            data_All5.append(together_data)
            label_All5.append('E')  # label person 5
    data_cup5 = ReadData1(data_All5)

    # Get the data from person one.
    for i in Range2:
        path_tmp = path_base_cat6 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)
        if i == 1:
            data_All6 = [together_data]
            label_All6.append('F')  # label person 6
        else:
            data_All6.append(together_data)
            label_All6.append('F')  # label person 6
    data_cup6 = ReadData1(data_All6)

    # Get the data from person two.
    for i in Range2:
        path_tmp = path_base_cat7 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)
        if i == 1:
            data_All7 = [together_data]
            label_All7.append('G')  # label person 7
        else:
            data_All7.append(together_data)
            label_All7.append('G')  # label person 7
    data_cup7 = ReadData1(data_All7)

    # # Get the data from person Three.
    for i in Range2:
        path_tmp = path_base_cat8 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)

        if i == 1:
            data_All8 = [together_data]
            label_All8.append('H')  # label person 8
        else:
            data_All8.append(together_data)
            label_All8.append('H')  # label person 8
    data_cup8 = ReadData1(data_All8)

    # Get the data from person Four.
    for i in Range2:
        path_tmp = path_base_cat9 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)

        if i == 1:
            data_All9 = [together_data]
            label_All9.append('I')  # label person 9
        else:
            data_All9.append(together_data)
            label_All9.append('I')  # label person 9
    data_cup9 = ReadData1(data_All9)

    # Get the data from person Five.
    for i in Range2:
        path_tmp = path_base_cat10 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)

        if i == 1:
            data_All10 = [together_data]
            label_All10.append('J')  # label person 10
        else:
            data_All10.append(together_data)
            label_All10.append('J')  # label person 10
    data_cup10 = ReadData1(data_All10)


    # Get the data from person one.
    # for i in Range2:
    #     path_tmp = path_base_cat11 + str(i) + '.mat'
    #     data_tmp = scio.loadmat(path_tmp)
    #     together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)
    #     if i == 1:
    #         data_All11 = [together_data]
    #         label_All11.append('K')  # label person 11
    #     else:
    #         data_All11.append(together_data)
    #         label_All11.append('K')  # label person 11
    # data_cup11 = ReadData1(data_All11)

    # # Get the data from person two.
    # for i in Range2:
    #     path_tmp = path_base_cat12 + str(i) + '.mat'
    #     data_tmp = scio.loadmat(path_tmp)
    #     together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)
    #     if i == 1:
    #         data_All12 = [together_data]
    #         label_All12.append('L')  # label person 12
    #     else:
    #         data_All12.append(together_data)
    #         label_All12.append('L')  # label person 12
    # data_cup12 = ReadData1(data_All12)

    # # Get the data from person Three.
    # for i in Range2:
    #     path_tmp = path_base_cat13 + str(i) + '.mat'
    #     data_tmp = scio.loadmat(path_tmp)
    #     together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)
    #
    #     if i == 1:
    #         data_All13 = [together_data]
    #         label_All13.append('M')  # label person 13
    #     else:
    #         data_All13.append(together_data)
    #         label_All13.append('M')  # label person 13
    # data_cup13 = ReadData1(data_All13)

    # Get the data from person Four.
    for i in Range2:
        path_tmp = path_base_cat14 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)

        if i == 1:
            data_All14 = [together_data]
            label_All14.append('N')  # label person 14
        else:
            data_All14.append(together_data)
            label_All14.append('N')  # label person 14
    data_cup14 = ReadData1(data_All14)

    # Get the data from person Five.
    for i in Range2:
        path_tmp = path_base_cat15 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)

        if i == 1:
            data_All15 = [together_data]
            label_All15.append('O')  # label person 15
        else:
            data_All15.append(together_data)
            label_All14.append('O')  # label person 15
    data_cup15 = ReadData1(data_All15)

    # Concatenate the data off three people.
    # data_cup = np.concatenate((data_cup1, data_cup3, data_cup5,
    #                            data_cup6, data_cup9, data_cup10,
    #                            data_cup7, data_cup8, data_cup14,
    #                            data_cup12, data_cup15), axis=0)
    DataSet_Split = "train" if training else "test"
    if DataSet_Split == "train":
        data_cup = np.concatenate((data_cup1, data_cup3, data_cup5,
                                   data_cup6, data_cup9, data_cup10), axis=0)
        cupLabels = np.concatenate((label_All1, label_All3, label_All5,
                                   label_All6, label_All9, label_All10), axis=0)
    elif DataSet_Split == "test":
        data_cup = np.concatenate((data_cup7, data_cup8, data_cup14,
                                data_cup4,  data_cup15), axis=0)  #data_cup4, data_cup13,data_cup11,data_cup2,
        cupLabels = np.concatenate((label_All7, label_All8, label_All14,
                                   label_All4, label_All15), axis=0)
    else:
        print("Dataload False")

    return data_cup, cupLabels



