# *ChestLive*: Fortifying Voice-based Authentication with Chest MotionBiometric on Smart Devices

## 环境

1、python 3.8.8

2、tensorflow-gpu 2.4.1

3、cuda 11.2

4、cudatoolkit 10.1.243

5、cudnn 7.6.5

6、sklearn 0.0

7、numpy  1.19.2

## 使用方法

### 训练

前往 [Fewshotchestmotion](https://github.com/MrWang98/ChestLive/tree/main/Fewshotchestmotion)，在 [Train.ipynb](https://github.com/MrWang98/ChestLive/blob/main/Fewshotchestmotion/Train.ipynb) 中设置相应参数，然后运行 [Train.ipynb](https://github.com/MrWang98/ChestLive/blob/main/Fewshotchestmotion/Train.ipynb).。一个h5格式的文件将在 [checkpoint](https://github.com/MrWang98/ChestLive/tree/main/checkpoint) 中创建，该文件是运行后保存的模型。

Dataset是训练模型的类。在初始化函数中ReadData.loadData()指定存放数据的文件，选择读取训练或测试数据。通过调用get_mini_dataset()函数获得。

在训练代码的最后面可以更改保存训练模型的文件夹及文件名。

### 测试

前往 [Fewshotchestmotion](https://github.com/MrWang98/ChestLive/tree/main/Fewshotchestmotion)，在 [test.py](https://github.com/MrWang98/ChestLive/blob/main/Fewshotchestmotion/test.py) 中设置与训练代码相同的参数。最后会将测试结果保存在指定文件里，可更改。

testDataset是测试的类。在初始化函数中会读取测试数据。

每次对比不同情况的数据时需要更改testDataset中init()第四行读取不同场景数据的代码，ReadData.loadData的第一个参数是读取数据所在的文件夹，需要和原始数据不同文件但是同样的训练集或者测试集。







