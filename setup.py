from setuptools import setup

setup(
    name='g3readers',
    version='1.0',
    description='Read Gadget2/3 snapshots, large simulations and catalogues',
    url='https://aragagnin.github.io',
    author='Antonio Ragagnin',
    author_email='antonio.ragagnin@inaf.it',
    license='GPLv3+',
    packages=['g3readers'],
    package_dir={'g3readers': '.'},
    zip_safe=False
)



