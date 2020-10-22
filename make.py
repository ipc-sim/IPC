import subprocess

# build ADD
runCommand = (
    'cd build\nmake -j 12')
subprocess.call([runCommand], shell=True)
