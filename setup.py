import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name = 'Lmixing',
    version = '1.0',
    author = 'Daniel Vrinceanu',
    author_email = 'daniel.vrinceanu@tsu.edu',
    url = 'https://github.com/Lmixing',
    license = 'Creative Commons Attribution-Noncommercial-Share Alike license',
    description = "A package to calculate L-mixing rate coefficients",
    long_description = open('README.md').read(),
    long_description_content_type = "text/markdown",
    packages = setuptools.find_packages(),
    install_requires = ['scipy', 'mpmath'],
    python_requires='>=3',
    classifiers = [
    'Topic :: Scientific/Engineering :: Astronomy',
    'Programming Language :: Phyton :: 3',
    'Topic :: Scientific/Engineering :: Physics',
    'Operating System :: OS Independent',
    ]
)
