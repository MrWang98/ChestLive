# ChestLive



运行[Fewshotchestmotion](https://github.com/MrWang98/ChestLive/tree/main/Fewshotchestmotion)文件夹中的[Train.ipynb](https://github.com/MrWang98/ChestLive/blob/main/Fewshotchestmotion/Train.ipynb)文件

note：Dataset是训练模型用的类，后面的testDataset是测试用的。每次对比不同情况的数据时需要更改testDataset中init()第四行读取不同场景数据的代码，ReadData.loadData的第一个参数是读取数据所在的文件夹，需要和原始数据不同文件但是同样的训练集或者测试集，例如某人的原始数据在set1/train，该原始数据对应的不同场景的数据在set2/train。init第三行第一个参数是set1，第四行第一个参数是set2，第二个参数都相同。还需更改get_mini_dataset()的if split段代码，该段代码的意思是，若抽取到某个人的数据用于训练，则将测试的数据更换为该用户在不同场景下的数据。

预训练好的模型在[checkpoint](https://github.com/MrWang98/ChestLive/tree/main/checkpoint)文件夹中：

1、其中model1是用split/train数据进行训练，test进行测试

2、model2与1的训练和测试相反

3、model_2way是用split/train中数据训练的2分类模
