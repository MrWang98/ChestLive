import os
import scipy.io as scio
import numpy as np

def readMatData(data):
    data_cup1_tmp = []
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

def loadAllData(finalDir,scaleOfTrainSet,training = True):
    dataRange = range(1,71)

    allData = {}    #所有处理后的数据，key是人名，value是对应的所有的数据
    index = 0       #person编码，第一个人是0，第二个人是1，……
    allLabels = {}  #所有标签，key是人名，value是对应的所有标签

    list = os.listdir(finalDir) #获取所有person的name
    for dir in list:
        matlabPath = os.path.join(finalDir, dir, "Matlab")  # 检查dir对应person文件夹下Matlab路径是否存在
        if os.path.exists(matlabPath):
            data = []
            labels = []
            pathBaseCat = os.path.join(matlabPath,"HeySiri") #获取dir对应person文件夹下HeySiri的路径
            for i in dataRange:
                tmpPath = pathBaseCat + str(i) + '.mat'
                tmpData = scio.loadmat(tmpPath)
                together_data = np.concatenate((tmpData['cor1_ccc2'], tmpData['cor2_ccc2']), axis=0)
                data.append(together_data)
                labels.append(chr(index+65))            # 将person编码使用ascii转成字母，0->A,1->B,……
            tmpDataCup = readMatData(data)              # 处理数据
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

    middleNum = round(dataValue.shape[0]*scaleOfTrainSet)   #数据分离的中间线，前scaleOfTrainSet比例的数据做训练集，剩下的做测试集

    dataSetSplit = "train" if training else "test"          #根据training选择训练集或测试集
    if dataSetSplit == "train":
        finalData = dataValue[0]
        cupLabels = labelsValue[0]
        for i in range(1,middleNum):
            finalData = np.concatenate((finalData, dataValue[i]), axis=0)
            cupLabels = np.concatenate((cupLabels,labelsValue[i]), axis=0)
    else:
        finalData = dataValue[middleNum]
        cupLabels = labelsValue[middleNum]
        for i in range(middleNum+1,len(dataValue)):
            finalData = np.concatenate((finalData, dataValue[i]), axis=0)
            cupLabels = np.concatenate((cupLabels,labelsValue[i]), axis=0)

    return finalData,cupLabels

def LoadData(training):
    finalDir = "../Train" #Final文件夹的路径，即所有收集以人名命名的文件夹所在的目录的路径
    scaleOfTrainSet = 2 / 3   #可在这里调整训练集占总数据的比例
    return loadAllData(finalDir,scaleOfTrainSet,training)  #第二个参数默认值是True，意思是取训练集


def loadTestData(testPath,scale):
    trainNum = round(70 * scale)
    trainDataRange = range(1, trainNum)
    testDataRange = range(trainNum,71)
    list = os.listdir(testPath)  # 获取所有person的name
    finalDir = testPath
    data = []
    labels = []
    trainData = []
    trainLabels = []
    testData = []
    testLabels = []
    for dir in list:
        matlabPath = os.path.join(finalDir, dir, "Matlab")  # 检查dir对应person文件夹下Matlab路径是否存在
        if os.path.exists(matlabPath):
            pathBaseCat = os.path.join(matlabPath, "HeySiri")  # 获取dir对应person文件夹下HeySiri的路径
            for i in trainDataRange:
                tmpPath = pathBaseCat + str(i) + '.mat'
                tmpData = scio.loadmat(tmpPath)
                together_data = np.concatenate((tmpData['cor1_ccc2'], tmpData['cor2_ccc2']), axis=0)
                trainData.append(together_data)
                trainLabels.append(dir)
    for dir in list:
        matlabPath = os.path.join(finalDir, dir, "Matlab")  # 检查dir对应person文件夹下Matlab路径是否存在
        if os.path.exists(matlabPath):
            pathBaseCat = os.path.join(matlabPath, "HeySiri")  # 获取dir对应person文件夹下HeySiri的路径
            for i in testDataRange:
                tmpPath = pathBaseCat + str(i) + '.mat'
                tmpData = scio.loadmat(tmpPath)
                together_data = np.concatenate((tmpData['cor1_ccc2'], tmpData['cor2_ccc2']), axis=0)
                testData.append(together_data)
                testLabels.append(dir)
    trainData = readMatData(trainData)
    testData = readMatData(testData)

    trainData = np.array(trainData)
    trainLabels = np.array(trainLabels)
    testData = np.array(testData)
    testLabels = np.array(testLabels)

    return trainData,trainLabels,testData,testLabels


if __name__ == '__main__':
    finalDir = "Test"
    scaleOfTrainSet = 2/3
    trainData, trainLabes, testData, testLabels = loadTestData(finalDir, scaleOfTrainSet)
    print("imageData:\t{},{}".format(trainData.shape,testData.shape))
    print("labels:\t\t{},{}".format(trainLabes.shape,testLabels.shape))

