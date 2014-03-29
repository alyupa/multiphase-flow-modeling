#!/usr/bin/pvpython

import paraview.simple
import os
import tempfile
import shutil
import sys

SRC_DIR = sys.argv[1]
DST_DIR = sys.argv[2]

if not os.path.isdir(DST_DIR):
    os.makedirs(DST_DIR)

for filename in os.listdir(SRC_DIR):
	print SRC_DIR + "/" + filename
	fileName, fileExtension = os.path.splitext(filename)

	temp_path = tempfile.mkdtemp()

	reader = paraview.simple.OpenDataFile(SRC_DIR + "/" + filename)
	writer = paraview.simple.CreateWriter(temp_path + "/" + fileName + ".vtm", reader)
	writer.UpdatePipeline()

	shutil.move(temp_path + "/" + fileName + "/" + fileName + "_0_0.vts", DST_DIR + "/" + fileName + ".vts")
	shutil.rmtree(temp_path)
