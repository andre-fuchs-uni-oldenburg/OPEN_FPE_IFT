import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="turbtools",
    version="22.6a3",
    description="Python Package for Statistical Analysis of Turbulence Data",
    author="Lille Turbulence Program 2022",
    maintainer = "Aakash Patil",
    license="GPL",
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords="turbulence fluids analysis stats",
    url="https://github.com/aakash30jan/turbtools",
    packages=setuptools.find_packages(),
    install_requires=['requests','scipy','numpy','h5py','matplotlib'],
    platforms=['any'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    entry_points={"console_scripts": ["turbtools = turbtools.turbtools:main"]},
)
 
