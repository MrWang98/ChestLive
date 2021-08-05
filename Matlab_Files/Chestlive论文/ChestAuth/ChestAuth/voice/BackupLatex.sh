#!/bin/bash

# `这个是用来执行终端命令的`  tar是打包命令 -zcvf是参数选项 打包后的文件名 后面可以以空格接一串文件或者文件夹。
tar -zcvf Match_VoiceandUndulate.`date +%Y%m%d`.back.tar.gz ./speech/*