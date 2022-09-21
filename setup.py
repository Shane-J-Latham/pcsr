import os
from skbuild import setup
import versioneer

setup(
    name="pcsr",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    zip_safe=False,
    packages=["pcsr",],
    package_data={'pcsr': [os.path.join("models", "*.obj")]},
    install_requires=[
        "numpy>=1.7",
    ],
    setup_requires=[
        "numpy>=1.7"
    ]
)
