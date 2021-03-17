'''
This build helper is modified to build IPC in c3d docker environment

To run this, you firstly have to log into c3d docker
'''

import subprocess

# build ADD
runCommand = (
    'mkdir build\ncd build\ncmake -DCMAKE_BUILD_TYPE=Release ..\nmake -j 12')
c3d_path = '$HOME/c3d/devtools/c3d'
store_base = '$HOME/c3d/bazel-out/k8-opt/bin'
cmds = ['cd $HOME/c3d\n',
        '{} build third_party/gmp_test\n'.format(c3d_path),
        'cd -\n',
        'rm -rf build\n',
        'mkdir build\n',
        'cd build\n',
        'cmake',
        '-DCMAKE_CXX_FLAGS=-std=c++17',
        '-DCMAKE_BUILD_TYPE=Release',
        '-DGMP_LIBRARIES={}/third_party/gmp/lib/libgmp.a'.format(store_base),
        '-DGMP_INCLUDE_DIRS={}/third_party/gmp/include'.format(store_base),
        '..\n',
        'make -j16']
#subprocess.call([runCommand], shell=True)
subprocess.call(' '.join(cmds), shell=True)
