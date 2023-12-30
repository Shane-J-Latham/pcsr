from skbuild import setup
import versioneer

setup(
    name="pcsr",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    zip_safe=False,
    packages=["pcsr",],
    package_data={'pcsr.models': ["data/*.obj"]},
    cmake_with_sdist=True
)
