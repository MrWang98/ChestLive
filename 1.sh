`rm -rf Set/test/*`
$(cp -r process_data/angle/$1/* Set/test/)
$(cp -r process_data/distance/$1/* Set/test/)
$(cp -r process_data/noise/$1/* Set/test/)
$(cp -r process_data/pose/$1/* Set/test/)
$(cp -r process_data/devices/$1* Set/test/)
$(cp -r process_data/wakeword/$1* Set/test/)
