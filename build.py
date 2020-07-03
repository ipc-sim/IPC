import subprocess

# build ADD
runCommand = (
    'mkdir build\ncd build\ncmake -DCMAKE_BUILD_TYPE=Release ..\nmake -j 12')
subprocess.call([runCommand], shell=True)
