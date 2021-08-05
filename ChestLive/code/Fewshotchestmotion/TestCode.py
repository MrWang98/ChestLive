import numpy as np
import scipy.io as scio

def ReadData1(data_All1):
    data_cup1_tmp = []
    data_cup1_tmp0 = np.concatenate((data_All1[0],
                                     np.zeros((2000 - data_All1[0].shape[0], 39))))
    # Used to get the dimension of data.
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

def LoadData(training):

    path_base_cat1 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/HuangJiaQian/Matlab/HeySiri'
    # path_base_cat1 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/GQY/Matlab/HeySiri'
    # path_base_cat1 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/DengYangTao/Matlab/HeySiri'
    path_base_cat2 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/BiHongLiang/Matlab/HeySiri'
    path_base_cat6 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/GQY/Matlab/HeySiri'
    path_base_cat7 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/ZhuTianLin/Matlab/HeySiri'
    path_base_cat8 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/XZQ/Matlab/HeySiri'

    cupLabels = []

    DataSet_Split = "train" if training else "test"
    if DataSet_Split == "train":
        Range1 = range(1, 71)
    else:
        Range1 = range(71, 101)
    print(Range1, Range1[0])
    # Get the data from person one.
    for i in Range1:
        path_tmp = path_base_cat1 + str(i) + '.mat'
        data_tmp = scio.loadmat(path_tmp)
        together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)
        if i == Range1[0]:
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
        if i == Range1[0]:
            data_All2 = [together_data]
            cupLabels.append('B')  # label person 2
        else:
            data_All2.append(together_data)
            cupLabels.append('B')  # label person 2
    data_cup2 = ReadData1(data_All2)

    # # Get the data from person one.
    # for i in Range1:
    #     path_tmp = path_base_cat6 + str(i) + '.mat'
    #     data_tmp = scio.loadmat(path_tmp)
    #     together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)
    #     if i == 1:
    #         data_All6 = [together_data]
    #         cupLabels.append('B')  # label person 6
    #     else:
    #         data_All6.append(together_data)
    #         cupLabels.append('B')  # label person 6
    # data_cup6 = ReadData1(data_All6)
    #
    # # Get the data from person two.
    # for i in Range1:
    #     path_tmp = path_base_cat7 + str(i) + '.mat'
    #     data_tmp = scio.loadmat(path_tmp)
    #     together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)
    #     if i == 1:
    #         data_All7 = [together_data]
    #         cupLabels.append('B')  # label person 7
    #     else:
    #         data_All7.append(together_data)
    #         cupLabels.append('B')  # label person 7
    # data_cup7 = ReadData1(data_All7)
    #
    # # Get the data from person Three.
    # for i in Range1:
    #     path_tmp = path_base_cat8 + str(i) + '.mat'
    #     data_tmp = scio.loadmat(path_tmp)
    #     together_data = np.concatenate((data_tmp['cor1_ccc2'], data_tmp['cor2_ccc2']), axis=0)
    #
    #     if i == 1:
    #         data_All8 = [together_data]
    #         cupLabels.append('B')  # label person 8
    #     else:
    #         data_All8.append(together_data)
    #         cupLabels.append('B')  # label person 8
    # data_cup8 = ReadData1(data_All8)

    # Concatenate the data off three people.
    data_cup = np.concatenate((data_cup1, data_cup2), axis=0)  # data_cup7,, data_cup6, data_cup8

    return data_cup, cupLabels

# data, labels = LoadData(training=True)

class Dataset:
    def __init__(self, training):
        self.data = {}
        data_train, labels_train = LoadData(training=training)
        for data_tmp, label_tmp in zip(data_train, labels_train):
            # data_tmp1 = np.array(data_tmp)
            data_tmp = np.expand_dims(data_tmp, axis=2)
            print(data_tmp.shape)
            if label_tmp not in self.data:
                self.data[label_tmp] = []
            self.data[label_tmp].append(data_tmp)
            labels = list(self.data.keys())

train_dataset = Dataset(training=True)
test_dataset = Dataset(training=False)




