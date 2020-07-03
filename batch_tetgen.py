import subprocess
import os

from os import listdir
from os.path import isfile, join

inputFolderPath = os.path.realpath('.') + '/input/triMeshes/verschoor/'
progPath = os.path.realpath('.') + '/src/Projects/MeshProcessing/meshprocessing'

onlyfiles = [
    f for f in listdir(inputFolderPath) if isfile(join(inputFolderPath, f))
]
for inputModelNameI in onlyfiles:
    runCommand = progPath + ' 0 ' + inputFolderPath + inputModelNameI + ' 3 1e10 0'
    if subprocess.call([runCommand], shell=True):
        continue
