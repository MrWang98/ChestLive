# ChestLive



运行[Fewshotchestmotion](https://github.com/MrWang98/ChestLive/tree/main/Fewshotchestmotion)文件夹中的[Train.ipynb](https://github.com/MrWang98/ChestLive/blob/main/Fewshotchestmotion/Train.ipynb)文件

note：Dataset是训练模型用的类，后面的testDataset是测试用的。每次对比不同情况的数据时需要更改testDataset中init()第四行读取不同场景数据的代码，ReadData.loadData的第一个参数是读取数据所在的文件夹，需要和原始数据不同文件但是同样的训练集或者测试集.

预训练好的模型在[checkpoint](https://github.com/MrWang98/ChestLive/tree/main/checkpoint)文件夹中：

1、其中model1是用split/train数据进行训练，test进行测试

2、model2与1的训练和测试相反

3、model_2way是数据训练的2分类模型
