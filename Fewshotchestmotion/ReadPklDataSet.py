"""
    Load the dataset with formats of pkl.
    type(),
"""

import pickle
# The path of DataSet.
path_train = '/Users/mengxue/PycharmProjects/AIModels/DataSet/omniglot/train_vinyals_aug90.pkl'
path_test = '/Users/mengxue/PycharmProjects/AIModels/DataSet/omniglot/test_vinyals_aug90.pkl'
path_val = '/Users/mengxue/PycharmProjects/AIModels/DataSet/omniglot/val_vinyals_aug90.pkl'
# # Load train data.
# f_train = open(path_train, 'rb')
# train_data = pickle.load(f_train, encoding='iso-8859-1')
# f_train.close()
# # Load test data.
# f_test = open(path_test, 'rb')
# test_data = pickle.load(f_test, encoding='iso-8859-1')
# f_test.close()
# Load val data.

# f_val = open(path_val, 'rb')
# val_data = pickle.load(f_val, encoding='iso-8859-1')
# f_val.close()

# Store data.
# path_val_trash = '/Users/mengxue/PycharmProjects/AIModels/DataSet/omniglot/val_Test_Trash.pkl'
# df_trash = open(path_val_trash,'wb')
# pickle.dump(val_data, df_trash)
# df_trash.close()
# Read data
path_val_trash = '/Users/mengxue/PycharmProjects/AIModels/DataSet/omniglot/val_Test_Trash.pkl'
df_trash = open(path_val_trash,'rb')
val_data = pickle.load(df_trash, encoding='iso-8859-1')
df_trash.close()



