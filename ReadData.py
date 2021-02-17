import os
import scipy.io as scio
import numpy as np

def readMatData(data):
    data_cup1_tmp = []
    a = data[0]
    b = np.zeros((2000 - data[0].shape[0], 39))
    data_cup1_tmp0 = np.concatenate((data[0],
                                     np.zeros((2000 - data[0].shape[0], 39))))
    nda = np.array(data).shape[0]
    Range1 = range(1, nda)
    for i in Range1:
        data_cup1_tmp = np.concatenate((data[i],
                                        np.zeros((2000 - data[i].shape[0], 39))))
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

def readData(path,startIndex):
    dataRange = range(1, 101)
    allData = {}  # 所有处理后的数据，key是人名，value是对应的所有的数据
    index = startIndex  # person编码，第一个人是0，第二个人是1，……
    allLabels = {}  # 所有标签，key是人名，value是对应的所有标签

    list = os.listdir(path)  # 获取所有person的name
    print("personlist:{}".format(list))
    for dir in list:
        print(dir)
        personPath = os.path.join(path, dir)
        if os.path.exists(personPath):
            fileList = os.listdir(personPath)
            # if len(dir.split("_"))==2 and dir.split("_")[1]=="attacker":
            #     fileList.sort(key=lambda x:int(x[7:-4]))
            data = []
            labels = []
            for file in fileList:
                if(file.split('.')[-1]=="mat"):
                    tmpPath = os.path.join(personPath,file)
                    # print(tmpPath)
                    tmpData = scio.loadmat(tmpPath)
                    together_data = np.concatenate((tmpData['cor1_ccc2'], tmpData['cor2_ccc2']), axis=0)
                    if(together_data.shape[0]>2000):
                        os.remove(tmpPath)
                    if(together_data.shape[0]<2000):
                        data.append(together_data)
                        labels.append(dir)  # 将person编码使用ascii转成字母，0->A,1->B,……
            tmpDataCup = readMatData(data)  # 处理数据
            allData[dir] = tmpDataCup
            allLabels[dir] = labels
            index += 1

    dataValue = []
    labelsValue = []
    for key in allData.keys():
        dataValue.append(allData[key])
        labelsValue.append(allLabels[key])
    dataValue = np.array(dataValue)
    labelsValue = np.array(labelsValue)
    return dataValue,labelsValue

def loadData(startPath, training = True):
    trainSetPath = os.path.join(startPath,"train")
    testSetPath = os.path.join(startPath, "test")
    length = len(os.listdir(trainSetPath))

    dataSetSplit = "train" if training else "test"          #根据training选择训练集或测试集
    if dataSetSplit == "train":
        dataValue,labelsValue = readData(trainSetPath,0)
        finalData = dataValue[0]
        cupLabels = labelsValue[0]
        for i in range(1,len(dataValue)):
            finalData = np.concatenate((finalData, dataValue[i]), axis=0)
            cupLabels = np.concatenate((cupLabels,labelsValue[i]), axis=0)
    elif dataSetSplit=="test":
        dataValue, labelsValue = readData(testSetPath,length)
        finalData = dataValue[0]
        cupLabels = labelsValue[0]
        for i in range(1,len(dataValue)):
            finalData = np.concatenate((finalData, dataValue[i]), axis=0)
            cupLabels = np.concatenate((cupLabels,labelsValue[i]), axis=0)

    return finalData,cupLabels

if __name__ == '__main__':
    train_data, train_labes = loadData(startPath="split2",training=True)
    test_data, test_labels = loadData(startPath="split2",training=False)
    print()


