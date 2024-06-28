from setuptools.command.build_ext import build_ext
from setuptools.command.build_py import build_py
from setuptools import setup
from setuptools import find_packages
from setuptools import Extension

try:
    import configparser
except ImportError:
    import ConfigParser as configparser

import os
import glob
import sys
import platform
import re

from tempfile import mkstemp
from shutil import move, copymode
from os import fdopen, remove


vrna_flags = [
    'with-svm',
    'with-naview',
    'with-gsl',
    'with-mpfr',
    'with-openmp',
    'with-simd',
]

vrna_options = [
    ('svm_macro', None,
     "Preprocessor macro activating SVM support."),
    ('svm_libs', None,
     "Additional libraries required when linking with SVM support."),
    ('gsl_macro', None,
     "Preprocessor macro activating GSL support."),
    ('gsl_libs', None,
     "Additional libraries required when linking with GSL support."),
    ('mpfr_macro', None,
     "Preprocessor macro activating MPFR support."),
    ('mpfr_libs', None,
     "Additional libraries required when linking with MPFR support."),
    ('naview_macro', None,
     "Preprocessor macro activating NAVIEW layout support."),
    ('openmp_macro', None,
     "Preprocessor macro activating OpenMP support."),
    ('openmp_cflags', None,
     "OpenMP compiler/linker flags."),
]

vrna_simd_cflags = {
    'SSE41': '-msse4.1',
    'AVX512': '-mavx512f',
}

vrna_simd_files = {
    'SSE41' : [
        'src/ViennaRNA/utils/higher_order_functions_sse41.c',
    ],
    'AVX512' : [
        'src/ViennaRNA/utils/higher_order_functions_avx512.c',
    ],
}

# readme
THIS_DIRECTORY = os.path.abspath(os.path.dirname(__file__))
README_FILE = os.path.join(THIS_DIRECTORY, "README.md")
LONG_DESCRIPTION = ""

remove_sections = [
  "Table of Contents",
  "Executable Programs",
  "Configuration",
  "Energy Parameters"
]

header_pat = re.compile(r'^(\#+)\s*(.*)$')

with open(README_FILE, encoding="utf-8") as readme:
    within_section_remove = False
    section_level = 0
    for l in readme:
        m = header_pat.match(l)
        if m:
            if within_section_remove:
                if len(m.group(1)) <= section_level:
                    if m.group(2) in remove_sections:
                        section_level = len(m.group(1))
                        continue
                    else:
                        within_section_remove = False
                else:
                    continue
            elif m.group(2) in remove_sections:
                section_level = len(m.group(1))
                within_section_remove = True
                continue
        elif within_section_remove:
            continue

        LONG_DESCRIPTION += l


def comment_lines(file_path, patterns, comment_char = "//"):
    #Create temp file
    fh, abs_path = mkstemp()
    with fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                replaced = False
                for p in patterns:
                    if line.startswith(p):
                        new_file.write(comment_char + line)
                        replaced = True
                        break
                if not replaced:
                    new_file.write(line)
    #Copy the file permissions from the old file to the new file
    copymode(file_path, abs_path)
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)


class vrna_build_ext(build_ext):
    user_options = build_ext.user_options[:]
    user_options.extend(vrna_options)

    boolean_options = build_ext.boolean_options[:]
    boolean_options.extend(vrna_flags)

    def __init__(self, *args, **kwargs):
        build_ext.__init__(self, *args, **kwargs)

    def run(self):
        build_ext.run(self)

    def initialize_options(self):
        self.parse_options()
        build_ext.initialize_options(self)

    def compiler_is_msvc(self):
        return self.get_compiler_name().lower().startswith('msvc')

    def compiler_is_mingw(self):
        return self.get_compiler_name().lower().startswith('mingw')

    def get_compiler_name(self):
        """Return the name of the C compiler used to compile extensions.

        If a compiler was not explicitly set (on the command line, for
        example), fall back on the default compiler.
        """
        if self.compiler:
            # distutils doesn't keep the type of self.compiler uniform; we
            # compensate:
            if isinstance(self.compiler, str):
                name = self.compiler
            else:
                name = self.compiler.compiler_type
        else:
            name = get_default_compiler()
        return name

    def build_extensions(self):
        for ext in self.extensions:
            ext.define_macros += [('VRNA_VERSION', "\"2.6.4\"")]
            # we include config.h, so let's tell the source files about it
            if not self.compiler_is_msvc():
                ext.define_macros += [('HAVE_CONFIG_H', None)]

            ext.libraries += ['m']

            if sys.platform == 'darwin':
                # always disable OpenMP on Mac OS X (until we resolve a good way to activate it)
                self.with_openmp = False
                ext.extra_compile_args += ["-stdlib=libc++",
                  "-mmacosx-version-min=10.14",
                  "-std=c++1z"]  # c++17 for dlib 19.24

                # deactivate SIMD implementations for ARM builds
                if platform.machine() == "arm64" or "arm64" in os.getenv("ARCHFLAGS", ""):
                    self.with_simd = False
                    # collect all SIMD specific files
                    simd_files = [f for tmp in vrna_simd_files.values() for f in tmp]
                    ext.sources = [s for s in ext.sources if s not in simd_files]
                else:
                    ext.define_macros += [('VRNA_WITH_SIMD_AVX512', None)]
                    ext.define_macros += [('VRNA_WITH_SIMD_SSE41', None)]
            elif self.compiler_is_msvc():
                self.with_openmp = False
                ext.define_macros += [('VRNA_WITH_SIMD_AVX512', None)]
                ext.define_macros += [('VRNA_WITH_SIMD_SSE41', None)]
                ext.sources.append("src/dlib-19.24/dlib/all/source.cpp")
            else:
                if (platform.machine() == "arm64" or "arm64" in os.getenv("ARCHFLAGS", "")) or \
                   (platform.machine() == "aarch64" or "aarch64" in os.getenv("ARCHFLAGS", "")):
                    self.with_simd = False
                    # collect all SIMD specific files
                    simd_files = [f for tmp in vrna_simd_files.values() for f in tmp]
                    ext.sources = [s for s in ext.sources if s not in simd_files]
                else:
                    ext.define_macros += [('VRNA_WITH_SIMD_AVX512', None)]
                    ext.define_macros += [('VRNA_WITH_SIMD_SSE41', None)]


            for sw in [f.replace('with-', '') for f in vrna_flags]:
                if getattr(self, 'with_' + sw):
                    if getattr(self, sw + '_macro', None):
                        ext.define_macros += [(getattr(self, sw + '_macro'), None)]
                    if getattr(self, sw + '_libs', None):
                        ext.libraries += [l for l in getattr(self, sw + '_libs').split()]
                    if getattr(self, sw + '_cflags', None):
                        ext.extra_compile_args += [f for f in getattr(self, sw + '_cflags').split()]
                        ext.extra_link_args += [f for f in getattr(self, sw + '_cflags').split()]

            if self.compiler_is_msvc():
                ext.libraries.remove("m")
                try:
                    ext.libraries.remove("stdc++")
                except:
                    pass


        build_ext.build_extensions(self)

    def build_extension(self, ext):
        # overwrite build_ext._compile() method to change compilee flags
        # at a per-file base. This supposedly doesn't work for MSVC, so
        # we only do it if we are sure that a ._compile() method exists
        #
        # Idea taken from https://stackoverflow.com/a/49772928/18609162
        if not self.compiler_is_msvc():
            original__compile = self.compiler._compile
            if self.with_simd:
                ext.extra_compile_args.extend([f for f in vrna_simd_cflags.values()])

                def new__compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
                    for simd_opt, simd_flag in vrna_simd_cflags.items():
                        if src not in vrna_simd_files[simd_opt]:
                            extra_postargs = [arg for arg in extra_postargs if arg != simd_flag]

                    # remove c++ flags for .c files
                    if src.lower().endswith('.c'):
                      extra_postargs = [arg for arg in extra_postargs if not arg.lower().startswith('-std')]

                    return original__compile(obj, src, ext, cc_args, extra_postargs, pp_opts)
            else:
                def new__compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
                    # remove c++ flags for .c files
                    if src.lower().endswith('.c'):
                      extra_postargs = [arg for arg in extra_postargs if not arg.lower().startswith('-std')]

                    return original__compile(obj, src, ext, cc_args, extra_postargs, pp_opts)

            self.compiler._compile = new__compile

            try:
                build_ext.build_extension(self, ext)
            finally:
                del self.compiler._compile
            return

        return build_ext.build_extension(self, ext)

    def parse_options(self, configfile = "setup.cfg"):
        pyproject_entry = 'build_ext'
        parser = configparser.ConfigParser()
        parser.read(configfile)

        for f in [f.replace('with-','') for f in vrna_flags]:
            switch = 'with_' + f

            # all options default to False
            setattr(self, switch, False)

            if parser.has_option(pyproject_entry, switch):
                setattr(self, switch, parser.getboolean(pyproject_entry, switch))

            if getattr(self, switch):
                # Check if macro name is given
                macro = f + '_macro'
                if parser.has_option(pyproject_entry, macro):
                    m = parser.get(pyproject_entry, macro)
                    if m:
                        setattr(self, macro, m)

                # Check if libraries are given
                libs = f + '_libs'
                if parser.has_option(pyproject_entry, libs):
                    l = [ll for ll in parser.get(pyproject_entry, libs).split()]
                    if l:
                        setattr(self, libs, l)

                # Check compile flags
                cflags = f + '_cflags'
                if parser.has_option(pyproject_entry, cflags):
                    fl = parser.get(pyproject_entry, cflags)
                    if fl:
                        setattr(self, cflags, fl)


def vrna_ext(name = "RNA._RNA"):
    define_macros = []
    libraries = []
    extra_compile_args = []
    extra_link_args = []
    includeDirs = []
    srcFiles = []

    # collect all RNAlib source files
    srcDir = 'src/ViennaRNA'
    globStr = "%s/*.c*" % srcDir
    files = glob.glob(globStr)
    srcFiles += files

    srcFiles += ['src/json/json.c']
    srcFiles += ['src/cephes/kn.c',
                 'src/cephes/const.c',
                 'src/cephes/mtherr.c',
                 'src/cephes/expn.c']

    for root, dirnames, filenames in os.walk(srcDir):
      for dirname in dirnames:
        absPath = os.path.join(root, dirname)
        #print('adding dir to path: %s' % absPath)
        globStr = "%s/*.c*" % absPath
        files = glob.glob(globStr)
        #print(files)
        srcFiles += files

    srcFiles += ['interfaces/Python/RNA_wrap.cpp']

    # append include paths required to build the Python interface
    includeDirs.append('.')
    includeDirs.append('src')
    includeDirs.append('src/ViennaRNA')
    includeDirs.append('src/dlib-19.24')
    includeDirs.append('src/cephes')

    if "src/libsvm-3.31" != "src/":
        srcFiles += ['src/libsvm-3.31/svm.cpp']
        includeDirs.append('src/libsvm-3.31')

    #print("includeDirs:")
    #print(includeDirs)
    #print("srcFiles:")
    #print(srcFiles)

    # comment out every cpp define in config.h that we want to handle
    # from within this setup.py script
    defines = ['#define VRNA_WITH_GSL',
               '#define VRNA_NR_SAMPLING_MPFR',
               '#define VRNA_WITH_SVM',
               '#define VRNA_WITH_NAVIEW_LAYOUT',
               '#define VRNA_WITH_OPENMP',
               '#define VRNA_WITH_SIMD_AVX512',
               '#define VRNA_WITH_SIMD_SSE41'
              ]

    comment_lines('config.h', defines)

    return Extension(
        name = name,
        sources = srcFiles,
        include_dirs = includeDirs,
        define_macros = define_macros,
        extra_compile_args = extra_compile_args,
        libraries = libraries,
        language="c++",
        extra_link_args = extra_link_args,
    )


# Build extensions before python modules,
# or the generated SWIG python files will be missing.
#class BuildPy(build_py):
#    def run(self):
#        self.run_command('build_ext')
#        super(build_py, self).run()



setup(
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    packages=['ViennaRNA', 'RNA'],
    ext_modules=[vrna_ext()],
    cmdclass={'build_ext': vrna_build_ext},
    py_modules = ['RNA.RNA'],
)
