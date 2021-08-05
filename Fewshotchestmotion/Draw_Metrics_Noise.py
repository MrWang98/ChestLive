import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pylab as pl


Noisedata = {'User1': {'50db': 100, '55db': 97.5, '65db': 97.5},
             'User2': {'50db': 100, '55db': 97.5, '65db': 95}}

Devicedata = {'Registered': {'4E:P30': 98.57, 'Mi10:S30': 100, 'Mi6:4E': 100},
              'Signin': {'4E:P30': 97.85, 'Mi10:S30': 97.85, 'Mi6:4E': 98.57}}

Directiondata = {'User1': {'Center': 96.66, 'Left': 96.66, 'Right': 94.28},
                 'User2': {'Center': 98.33, 'Left': 98.33, 'Right': 96.66}}

Arousedata = {'Aourse Words': {'Alexa': 97.5, 'HeySiri': 98.21, 'Guess': 97.49,  'OkGoogle': 99.28, 'Choice': 97.85, 'Describe': 95}}

DistanceData = pd.DataFrame([92.50, 97.50, 97.50, 95.00, 92.50, 89.99],
                            columns=['Distance'],
                            index=np.arange(5, 35, 5))
MimicData = {'User1': {'False accept rate': 0.83, 'False reject rate': 0.83},
             'User2': {'False accept rate': 2.73, 'False reject rate': 2.73},
             'User3': {'False accept rate': 1.36, 'False reject rate': 1.36}}
# 1 different direction.
# 2 different noise.
# 3 different device.
# 4 different arousing words.
# 5 different arousing distance.
# 6 MimicData

storePic = 1
flatui = ["#386AB1", "#7293CA", "#CBE0F7"]
palette1 = sns.color_palette(flatui)

different_figure = 0

if different_figure == 1 or different_figure == 0:
    figure1_palette = sns.color_palette(palette=sns.color_palette("Blues_r", 3))
    figure1_data = pd.DataFrame(Directiondata)
    with sns.color_palette(palette1):  # figure2_palette
        fig1 = plt.figure(1)
        ax1 = figure1_data.plot.bar()
    plt.xlabel("Different direction", fontproperties='Times New Roman', size=18)
    pl.xticks(rotation=360)
    plt.ylabel("Accuracy (%)", fontproperties='Times New Roman', size=18)
    plt.yticks(fontproperties='Times New Roman', size=18)
    plt.xticks(fontproperties='Times New Roman', size=18)
    font = {'family': 'Times New Roman', 'weight': 'normal', 'size': 18}
    plt.legend(['User1', 'User2'], prop=font, loc='lower right', handletextpad=0.05)
    plt.grid(axis='y', ls='--')
    ax1.set_axisbelow(True)
    if (storePic == 1):
        plt.savefig('./Figures/ImpactofDirection.png', dpi=400, bbox_inches='tight')
    plt.show()
#
if different_figure == 2 or different_figure == 0:
    figure2_palette = sns.color_palette(palette=sns.color_palette("Blues_r", 3))
    figure2_data = pd.DataFrame(Noisedata)
    # figure2_data.index = ['50db', '55db', '65db']
    with sns.color_palette(palette1):  # figure2_palette
        fig2 = plt.figure(2)
        ax2 = figure2_data.plot.bar()
    plt.xlabel("Ambient noise level", fontproperties='Times New Roman', size=18)
    pl.xticks(rotation=360)
    plt.ylabel("Accuracy (%)", fontproperties='Times New Roman', size=18)
    plt.yticks(fontproperties='Times New Roman', size=18)
    plt.xticks(fontproperties='Times New Roman', size=18)
    font = {'family': 'Times New Roman', 'weight': 'normal', 'size': 18}
    plt.legend(['User1', 'User2'], prop=font, loc='lower right', handletextpad=0.05)
    plt.grid(axis='y', ls='--')
    ax2.set_axisbelow(True)
    if (storePic == 1):
        plt.savefig('./Figures/ImpactofNoise.png', dpi=400, bbox_inches='tight')
    plt.show()
#
if different_figure == 3 or different_figure == 0:
    figure3_palette = sns.color_palette(palette=sns.color_palette("Blues_r", 3))
    figure3_data = pd.DataFrame(Devicedata)
    with sns.color_palette(palette1):  # figure2_palette
        fig3 = plt.figure(3)
        ax3 = figure3_data.plot.bar()

    plt.xlabel("Different device", fontproperties='Times New Roman', size=18)
    plt.ylabel("Accuracy (%)", fontproperties='Times New Roman', size=18)
    pl.xticks(rotation=360)
    plt.xticks(fontproperties='Times New Roman', size=18)
    plt.yticks(fontproperties='Times New Roman', size=18)
    font = {'family': 'Times New Roman', 'weight': 'normal', 'size': 18}
    plt.legend(['Registered device', 'Sign in device'], prop=font, loc='lower right', handletextpad=0.05)
    plt.grid(axis='y', ls='--')
    ax3.set_axisbelow(True)
    if(storePic == 1):
        plt.savefig('./Figures/ImpactofDevice.png', dpi=400, bbox_inches='tight')
    plt.show()

if different_figure == 4 or different_figure == 0:
    figure4_palette = sns.color_palette(palette=sns.color_palette("Blues_r", 3))
    figure4_data = pd.DataFrame(Arousedata)
    with sns.color_palette(palette1):  # figure2_palette
        fig4 = plt.figure(4)
        ax4 = figure4_data.plot.bar()
    ax4.legend_.remove()
    plt.xlabel("Different arouse words", fontproperties='Times New Roman', size=18)
    pl.xticks(rotation=360)
    plt.ylabel("Accuracy (%)", fontproperties='Times New Roman', size=18)
    plt.yticks(fontproperties='Times New Roman', size=18)
    plt.xticks(fontproperties='Times New Roman', size=16)
    font = {'family': 'Times New Roman', 'weight': 'normal', 'size': 18}
    plt.grid(axis='y', ls='--')
    ax4.set_axisbelow(True)
    if (storePic == 1):
        plt.savefig('./Figures/ImpactofArouseword.png', dpi=400, bbox_inches='tight')
    plt.show()

if different_figure == 5 or different_figure == 0:
    figure5_palette = sns.color_palette(palette=sns.color_palette("Blues_r", 3))
    figure5_data = DistanceData
    with sns.color_palette(palette1):  # figure2_palette
        fig5 = plt.figure(5)
        ax5 = figure5_data.plot.bar()
    pl.xticks(rotation=360)
    ax5.legend_.remove()
    plt.xlabel("Distance (cm)", fontproperties='Times New Roman', size=18)
    plt.ylabel("Accuracy (%)", fontproperties='Times New Roman', size=18)
    plt.yticks(fontproperties='Times New Roman', size=18)
    plt.xticks(fontproperties='Times New Roman', size=18)
    font = {'family': 'Times New Roman', 'weight': 'normal', 'size': 18}
    plt.grid(axis='y', ls='--')
    ax5.set_axisbelow(True)
    if (storePic == 1):
        plt.savefig('./Figures/ImpactofDistance.png', dpi=400, bbox_inches='tight')
    plt.show()

if different_figure == 6 or different_figure == 0:
    figure6_palette = sns.color_palette(palette=sns.color_palette("Blues_r", 3))
    figure6_data = pd.DataFrame(MimicData)
    with sns.color_palette(palette1):  # figure2_palette
        fig6 = plt.figure(6)
        ax6 = figure6_data.plot.bar()
    pl.xticks(rotation=360)
    plt.ylabel("Percentage (%)", fontproperties='Times New Roman', size=18)
    plt.xticks(fontproperties='Times New Roman', size=18)
    plt.yticks(fontproperties='Times New Roman', size=18)
    font = {'family': 'Times New Roman', 'weight': 'normal', 'size': 18}
    plt.legend(['Authenticated u1', 'Authenticated u2', 'Authenticated u3'], prop=font, loc='upper center', handletextpad=0.05)
    plt.grid(axis='y', ls='--')
    ax6.set_axisbelow(True)
    if (storePic == 1):
        plt.savefig('./Figures/MimicAttack.png', dpi=300, bbox_inches='tight')
    plt.show()