#!/usr/bin/env python


from setuptools import Extension, setup
from Cython.Build import cythonize


extensions = [
    Extension("klib.kseq",
              ["src/kseq.pyx"],
              include_dirs=['src'],
              libraries=["z"]),
]

extensions = cythonize(extensions, include_path=['src'])

version = "0.6.6"
metadata = dict(
    name='klib',
    version=version,
    description="Python wrapping of Heng Li's klib",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Cython",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    url="https://github.com/tanghaibao/bio-pipeline/tree/master/klib",
    author=("Haibao Tang"),
    author_email="tanghaibao@gmail.com",
    ext_modules=extensions,
    license="BSD",
    packages=["klib"],
    package_dir={"klib": "src"},
    requires=['cython (>=0.21)'],
)

if __name__ == '__main__':
    setup(**metadata)
