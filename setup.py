import fnmatch
from setuptools import setup
from setuptools.command.build_py import build_py as build_py_orig

excluded =  ['g3pp.py','g3maps.py','test_g3read_ids.py','test_g3read.py', 'setup.py','c2pap_batch.py','g3read_units']

class build_py(build_py_orig):
    def find_package_modules(self, package, package_dir):
        modules = super().find_package_modules(package, package_dir)
        return [
            (pkg, mod, file)
            for (pkg, mod, file) in modules
            if not any(fnmatch.fnmatchcase(file, pat=pattern) for pattern in excluded)
        ]


setup(
    cmdclass={'build_py': build_py},
    name='g3read',
    version='1.0',
    description='Read Gadget2/3 snapshots, large simulations and catalogues',
    url='https://aragagnin.github.io',
    author='Antonio Ragagnin',
    author_email='antonio.ragagnin@inaf.it',
    license='GPLv3+',
    packages=['.'],
    zip_safe=False
)



