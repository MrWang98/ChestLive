import torch
import torch.nn as nn	# 各种层类型的实现
import torch.nn.functional as F	# 各中层函数的实现，与层类型对应，如：卷积函数、池化函数、归一化函数等等
import torch.optim as optim	# 实现各种优化算法的包
from torchvision import datasets, transforms

class ChestNet(nn.Module):
    def __init__(self):
        super(ChestNet, self).__init__()
        self.conv = torch.nn.Sequential()
        self.conv.add_module("conv1", torch.nn.Conv2d(3, 32, 3, 1, 1))
        self.conv.add_module("relu1", torch.nn.ReLU())
        self.conv.add_module("pool1", torch.nn.MaxPool2d(2))

        self.dense = torch.nn.Sequential()
        self.dense.add_module("dense1", torch.nn.Linear(32 * 3 * 3, 128))
        self.dense.add_module("relu2", torch.nn.ReLU())
        self.dense.add_module("dense2", torch.nn.Linear(128, 10))

    def forward(self, x):
        conv_out = self.conv1(x)
        res = conv_out.view(conv_out.size(0), -1)
        out = self.dense(res)
        return out

print("Method 3:")
model3 = ChestNet()
print(model3)
