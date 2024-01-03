from skbuild import setup
from setuptools import find_packages
import versioneer

setup(
    name="pcsr",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    zip_safe=False,
    packages=find_packages(),
    package_data={'pcsr': ["models/data/*.obj"], 'pcsr.models.data': ["*.obj"]},
    cmake_with_sdist=True
)
