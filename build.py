import subprocess

# build ADD
runCommand = (
    'mkdir build\ncd build\nrm CMakeCache.txt\ncmake -DCMAKE_BUILD_TYPE=Release ..\nmake -j 12')
subprocess.call([runCommand], shell=True)
