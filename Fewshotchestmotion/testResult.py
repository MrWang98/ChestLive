2021-01-13 21:17:05.446035: I tensorflow/core/platform/cpu_feature_guard.cc:143] Your CPU supports instructions that this TensorFlow binary was not compiled to use: AVX2 AVX512F FMA
2021-01-13 21:17:05.473564: I tensorflow/core/platform/profile_utils/cpu_utils.cc:102] CPU Frequency: 2400000000 Hz
2021-01-13 21:17:05.476689: I tensorflow/compiler/xla/service/service.cc:168] XLA service 0x7f3c64000b60 initialized for platform Host (this does not guarantee that XLA will be used). Devices:
2021-01-13 21:17:05.476741: I tensorflow/compiler/xla/service/service.cc:176]   StreamExecutor device (0): Host, Default Version
2021-01-13 21:17:05.481764: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcuda.so.1
2021-01-13 21:17:06.623682: I tensorflow/compiler/xla/service/service.cc:168] XLA service 0x55621b539120 initialized for platform CUDA (this does not guarantee that XLA will be used). Devices:
2021-01-13 21:17:06.623764: I tensorflow/compiler/xla/service/service.cc:176]   StreamExecutor device (0): GeForce RTX 3080, Compute Capability 8.6
2021-01-13 21:17:06.625739: I tensorflow/core/common_runtime/gpu/gpu_device.cc:1561] Found device 0 with properties: 
pciBusID: 0000:d5:00.0 name: GeForce RTX 3080 computeCapability: 8.6
coreClock: 1.71GHz coreCount: 68 deviceMemorySize: 9.77GiB deviceMemoryBandwidth: 707.88GiB/s
2021-01-13 21:17:06.626173: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcudart.so.10.1
2021-01-13 21:17:06.629810: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcublas.so.10
2021-01-13 21:17:06.633153: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcufft.so.10
2021-01-13 21:17:06.633658: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcurand.so.10
2021-01-13 21:17:06.637224: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcusolver.so.10
2021-01-13 21:17:06.639328: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcusparse.so.10
2021-01-13 21:17:06.645702: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcudnn.so.7
2021-01-13 21:17:06.647718: I tensorflow/core/common_runtime/gpu/gpu_device.cc:1703] Adding visible gpu devices: 0
2021-01-13 21:17:06.647772: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcudart.so.10.1
2021-01-13 21:17:06.648912: I tensorflow/core/common_runtime/gpu/gpu_device.cc:1102] Device interconnect StreamExecutor with strength 1 edge matrix:
2021-01-13 21:17:06.648932: I tensorflow/core/common_runtime/gpu/gpu_device.cc:1108]      0 
2021-01-13 21:17:06.648943: I tensorflow/core/common_runtime/gpu/gpu_device.cc:1121] 0:   N 
2021-01-13 21:17:06.650920: I tensorflow/core/common_runtime/gpu/gpu_device.cc:1247] Created TensorFlow device (/job:localhost/replica:0/task:0/device:GPU:0 with 4003 MB memory) -> physical GPU (device: 0, name: GeForce RTX 3080, pci bus id: 0000:d5:00.0, compute capability: 8.6)
2021-01-13 21:17:15.875723: I tensorflow/core/common_runtime/gpu/gpu_device.cc:1561] Found device 0 with properties: 
pciBusID: 0000:d5:00.0 name: GeForce RTX 3080 computeCapability: 8.6
coreClock: 1.71GHz coreCount: 68 deviceMemorySize: 9.77GiB deviceMemoryBandwidth: 707.88GiB/s
2021-01-13 21:17:15.875823: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcudart.so.10.1
2021-01-13 21:17:15.875839: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcublas.so.10
2021-01-13 21:17:15.875854: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcufft.so.10
2021-01-13 21:17:15.875866: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcurand.so.10
2021-01-13 21:17:15.875893: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcusolver.so.10
2021-01-13 21:17:15.875905: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcusparse.so.10
2021-01-13 21:17:15.875917: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcudnn.so.7
2021-01-13 21:17:15.876960: I tensorflow/core/common_runtime/gpu/gpu_device.cc:1703] Adding visible gpu devices: 0
2021-01-13 21:17:15.878105: I tensorflow/core/common_runtime/gpu/gpu_device.cc:1561] Found device 0 with properties: 
pciBusID: 0000:d5:00.0 name: GeForce RTX 3080 computeCapability: 8.6
coreClock: 1.71GHz coreCount: 68 deviceMemorySize: 9.77GiB deviceMemoryBandwidth: 707.88GiB/s
2021-01-13 21:17:15.878143: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcudart.so.10.1
2021-01-13 21:17:15.878155: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcublas.so.10
2021-01-13 21:17:15.878166: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcufft.so.10
2021-01-13 21:17:15.878175: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcurand.so.10
2021-01-13 21:17:15.878185: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcusolver.so.10
2021-01-13 21:17:15.878195: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcusparse.so.10
2021-01-13 21:17:15.878205: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcudnn.so.7
2021-01-13 21:17:15.879288: I tensorflow/core/common_runtime/gpu/gpu_device.cc:1703] Adding visible gpu devices: 0
2021-01-13 21:17:15.879329: I tensorflow/core/common_runtime/gpu/gpu_device.cc:1102] Device interconnect StreamExecutor with strength 1 edge matrix:
2021-01-13 21:17:15.879335: I tensorflow/core/common_runtime/gpu/gpu_device.cc:1108]      0 
2021-01-13 21:17:15.879341: I tensorflow/core/common_runtime/gpu/gpu_device.cc:1121] 0:   N 
2021-01-13 21:17:15.880393: I tensorflow/core/common_runtime/gpu/gpu_device.cc:1247] Created TensorFlow device (/job:localhost/replica:0/task:0/device:GPU:0 with 4003 MB memory) -> physical GPU (device: 0, name: GeForce RTX 3080, pci bus id: 0000:d5:00.0, compute capability: 8.6)
2021-01-13 21:23:17.339809: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcublas.so.10
2021-01-13 21:25:06.041939: I tensorflow/stream_executor/platform/default/dso_loader.cc:44] Successfully opened dynamic library libcudnn.so.7
2021-01-13 21:43:30.245495: W tensorflow/stream_executor/gpu/asm_compiler.cc:81] Running ptxas --version returned 256
2021-01-13 21:43:30.322815: W tensorflow/stream_executor/gpu/redzone_allocator.cc:314] Internal: ptxas exited with non-zero error code 256, output: 
Relying on driver to perform ptx compilation. 
Modify $PATH to customize ptxas location.
This message will be only logged once.
Using TensorFlow backend.
/home/xuemeng/Proj/ChestMotion/Fewshotchestmotion/Func_ReadAllData.py:9: VisibleDeprecationWarning: Creating an ndarray from ragged nested sequences (which is a list-or-tuple of lists-or-tuples-or ndarrays with different lengths or shapes) is deprecated. If you meant to do this, you must specify 'dtype=object' when creating the ndarray
  nda = np.array(data).shape[0]
test.py:70: DeprecationWarning: elementwise comparison failed; this will raise an error in the future.
  accAll+=(preds==labels).sum()
开始测试：
K…………
Traceback (most recent call last):
  File "test.py", line 70, in <module>
    accAll+=(preds==labels).sum()
AttributeError: 'bool' object has no attribute 'sum'
