# from keras import backend as K
import tensorflow as tf

def cal_base(y_true, y_pred):
    y_pred_positive = tf.round(y_pred)
    y_pred_negative = 1 - y_pred_positive

    y_positive = tf.round(y_true)
    y_negative = 1 - y_positive

    TP = tf.reduce_sum(y_positive * y_pred_positive)
    TN = tf.reduce_sum(y_negative * y_pred_negative)

    FP = tf.reduce_sum(y_negative * y_pred_positive)
    FN = tf.reduce_sum(y_positive * y_pred_negative)

    return TP, TN, FP, FN

def Modle_TP(y_true, y_pred):
    TP, TN, FP, FN = cal_base(y_true, y_pred)
    return TP
def Modle_TN(y_true, y_pred):
    TP, TN, FP, FN = cal_base(y_true, y_pred)
    return TN
def Modle_FP(y_true, y_pred):
    TP, TN, FP, FN = cal_base(y_true, y_pred)
    return FP
def Modle_FN(y_true, y_pred):
    TP, TN, FP, FN = cal_base(y_true, y_pred)
    return FN

def acc(y_true, y_pred):
    TP, TN, FP, FN = cal_base(y_true, y_pred)
    ACC = (TP + TN) / (TP + FP + FN + TN ) #+ K.epsilon()
    return ACC


def sensitivity(y_true, y_pred):
    """ recall """
    TP, TN, FP, FN = cal_base(y_true, y_pred)
    SE = TP/(TP + FN ) #+ K.epsilon()
    return SE


def precision(y_true, y_pred):
    TP, TN, FP, FN = cal_base(y_true, y_pred)
    PC = TP/(TP + FP ) #+ K.epsilon()
    return PC


def specificity(y_true, y_pred):
    TP, TN, FP, FN = cal_base(y_true, y_pred)
    SP = TN / (TN + FP ) #+ K.epsilon()
    return SP


def f1_socre(y_true, y_pred):
    SE = sensitivity(y_true, y_pred)
    PC = precision(y_true, y_pred)
    F1 = 2 * SE * PC / (SE + PC ) # + K.epsilon()
    return F1
