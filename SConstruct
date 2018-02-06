# Copyright (c) 2016, Kazuki Yoshizoe
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# AREDISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

from __future__ import print_function

import os
import commands
import ConfigParser

# todo: set compiler name (CXX) by command line options
# default values for CXX are g++ and mpicxx

config = ConfigParser.SafeConfigParser()
config.read('local.cfg')

orig_env = Environment(ENV=os.environ)

if config.has_option('compilers', 'option'):
    print( config.get('compilers', 'option') )
    orig_env.Append(CXXFLAGS=config.get('compilers', 'option'))
if config.has_option('compilers', 'libs'):
    print( config.get('compilers', 'libs') )
    libs_str_list = config.get('compilers', 'libs').split()
    for lib in libs_str_list:
        orig_env.Append(LIBS=[lib])
if config.has_option('paths', 'include'):
    print( config.get('paths', 'include') )
    inc_str_list = config.get('paths', 'include').split()
    for inc in inc_str_list:
        orig_env.Append(CPPPATH=[Dir([inc])])
if config.has_option('paths', 'library'):
    print( config.get('paths', 'library') )
    libs_str_list = config.get('paths', 'library').split()
    for lib in libs_str_list:
        orig_env.Append(LIBPATH=[Dir([lib])])

debug = ARGUMENTS.get('debug', 0)
log = ARGUMENTS.get('log', 0)

orig_env['ARFLAGS'] = 'crsv'
build_dir_base = 'build'


def AddDefEnv(e):
    e.Append(CXXFLAGS='-Wall -Wextra')
    e.Append(LIBS=['gflags'])
    e.Append(LIBS=['pthread'])


def Build(mode, compile_flags):
    build_dir = os.path.join(build_dir_base, mode)
    single_env_base = orig_env.Clone()
    if config.has_option('compilers', 'single'):
        single_env_base['CXX'] = config.get('compilers', 'single')
    else:
        single_env_base['CXX'] = 'g++'

    # src
    env = single_env_base.Clone()
    env.Append(CXXFLAGS=compile_flags)
    AddDefEnv(env)
    Export('env')
    VariantDir(os.path.join(build_dir, 'src'), 'src', duplicate=0)
    SConscript(os.path.join(build_dir, 'src', 'SConscript'))

    # tests
    test_dir = os.path.join(build_dir, 'tests')
    env = single_env_base.Clone()
    env.Append(CXXFLAGS=compile_flags)
    env.Append(LIBPATH=['#' + os.path.join(build_dir, 'src'),
                        '#' + test_dir])
    env.Append(CPPPATH=['#src', '#tests'])
    env.Append(LIBS=['lampsearch'])

    # make gtest.a and gtest_main.a
    env.StaticLibrary('#' + os.path.join(test_dir, 'gtest'),
                      ['#' + os.path.join(test_dir, 'gtest', 'gtest-all.cc')])
    env.StaticLibrary('#' + os.path.join(test_dir, 'gtest_main'),
                      ['#' + os.path.join(test_dir, 'gtest', 'gtest_main.cc')])
    env.Append(LIBS=['gtest', 'gtest_main'])

    AddDefEnv(env)
    Export('env')
    VariantDir(os.path.join(build_dir, 'tests'),
               'tests', duplicate=0)
    SConscript(os.path.join(build_dir, 'tests', 'SConscript'))

    # main
    env = single_env_base.Clone()
    env.Append(CXXFLAGS=compile_flags)
    env.Append(LIBPATH=['#' + os.path.join(build_dir, 'src')])
    env.Append(CPPPATH=['#src'])
    env.Append(LIBS=['lampsearch'])
    AddDefEnv(env)
    Export('env')
    VariantDir(os.path.join(build_dir, 'main'), 'main', duplicate=0)
    SConscript(os.path.join(build_dir, 'main', 'SConscript'))


def MPIBuild(mode, compile_flags):
    mpi_build_dir = os.path.join(mpi_build_dir_base, mode)
    mpi_env_base = orig_env.Clone()
    if config.has_option('compilers', 'parallel'):
        mpi_env_base['CXX'] = config.get('compilers', 'parallel')
    else:
        mpi_env_base['CXX'] = 'mpicxx'

    # src
    env = mpi_env_base.Clone()
    env.Append(CXXFLAGS=compile_flags)
    AddDefEnv(env)
    Export('env')
    VariantDir(os.path.join(mpi_build_dir, 'src'), 'src', duplicate=0)
    SConscript(os.path.join(mpi_build_dir, 'src', 'SConscript'))

    # mp-src
    env = mpi_env_base.Clone()
    env.Append(CXXFLAGS=compile_flags)
    env.Append(LIBPATH=['#' + os.path.join(mpi_build_dir, 'src')])
    env.Append(CPPPATH=['#mp-src', '#src'])
    env.Append(LIBS=['lampsearch'])

    AddDefEnv(env)
    Export('env')
    VariantDir(os.path.join(mpi_build_dir, 'mp-src'), 'mp-src', duplicate=0)
    SConscript(os.path.join(mpi_build_dir, 'mp-src', 'SConscript'))

    # tests
    test_dir = os.path.join(mpi_build_dir, 'tests')
    env = mpi_env_base.Clone()
    env.Append(CXXFLAGS=compile_flags)
    env.Append(LIBPATH=['#' + os.path.join(mpi_build_dir, 'mp-src'),
                        '#' + os.path.join(mpi_build_dir, 'src'),
                        '#' + os.path.join(mpi_build_dir, 'tests')])
    env.Append(CPPPATH=['#mp-src', '#src', '#tests'])
    env.Append(LIBPATH=['#' + os.path.join(mpi_build_dir, 'src')])
    env.Append(LIBS=['mplampsearch', 'lampsearch'])

    # make gtest.a and gtest_main.a
    env.StaticLibrary('#' + os.path.join(test_dir, 'gtest'),
                      ['#' + os.path.join(test_dir, 'gtest', 'gtest-all.cc')])
    env.StaticLibrary('#' + os.path.join(test_dir, 'gtest_main'),
                      ['#' + os.path.join(test_dir, 'gtest', 'gtest_main.cc')])
    env.Append(LIBS=['gtest', 'gtest_main'])

    AddDefEnv(env)
    Export('env')
    VariantDir(os.path.join(mpi_build_dir, 'tests'), 'tests', duplicate=0)
    SConscript(os.path.join(mpi_build_dir, 'tests', 'SConscript'))

    # mp-tests
    test_dir = os.path.join(mpi_build_dir, 'mp-tests')
    env = mpi_env_base.Clone()
    env.Append(CXXFLAGS=compile_flags)
    env.Append(LIBPATH=['#' + os.path.join(mpi_build_dir, 'mp-src'),
                        '#' + os.path.join(mpi_build_dir, 'src'),
                        '#' + os.path.join(mpi_build_dir, 'mp-tests'),
                        '#' + os.path.join(mpi_build_dir, 'tests')])
    env.Append(CPPPATH=['#mp-src', '#src', '#tests'])
    env.Append(LIBS=['mplampsearch', 'lampsearch'])
    env.Append(LIBS=['gtest', 'gtest_main'])

    AddDefEnv(env)
    Export('env')
    VariantDir(os.path.join(mpi_build_dir, 'mp-tests'),
               'mp-tests', duplicate=0)
    SConscript(os.path.join(mpi_build_dir, 'mp-tests', 'SConscript'))

    # mp-main
    env = mpi_env_base.Clone()
    env.Append(CXXFLAGS=compile_flags)
    env.Append(LIBPATH=['#' + os.path.join(mpi_build_dir, 'mp-src'),
                        '#' + os.path.join(mpi_build_dir, 'src')])
    env.Append(CPPPATH=['#mp-src', '#src'])
    env.Append(LIBS=['mplampsearch', 'lampsearch'])
    AddDefEnv(env)
    Export('env')
    VariantDir(os.path.join(mpi_build_dir, 'mp-main'), 'mp-main', duplicate=0)
    SConscript(os.path.join(mpi_build_dir, 'mp-main', 'SConscript'))

    if mode == 'opt':
        Command("mp-lamp", os.path.join(mpi_build_dir, "mp-main",
                                        "mp-lamp"), Copy("$TARGET", "$SOURCE"))


# opt build
Build(mode='opt', compile_flags='-O3 -DNDEBUG')
# dbg build
if int(debug):
    Build(mode='dbg', compile_flags='-g')


# profile build for gprof? google profiler?

# --------
# mpi build settings

mpi_build_dir_base = 'mp_build'

# mpi opt build
MPIBuild(mode='opt', compile_flags='-O3 -DNDEBUG -DOPTLOG')

# mpi dbg build
if int(debug):
    MPIBuild(mode='dbg', compile_flags='-g')

# mpi log build
if int(log):
    MPIBuild(mode='log', compile_flags='-O3 -DNDEBUG')
