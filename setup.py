from setuptools import setup, find_packages


setup(
    name='pypolycontain',
    author='Sadra Sadraddini',
    description='A python package for polytopic objects, operations, and containment encodings',
    author_email='sadra@mit.edu',
    version='1.3',
    packages=['pypolycontain'], #fix
    long_description='A python package for polytopic objects, operations, and containment encodings' 
    +'\nPlease refer to https://pypolycontain.readthedocs.io/en/latest/ for documentation' ,
#    package_dir={'pypolycontain': 'lib'},
    license='MIT',
    url="https://github.com/sadraddini/pypolycontain",
    zip_safe=False,
)
