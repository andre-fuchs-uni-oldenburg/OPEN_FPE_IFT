import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name="turbtools",
    version="22.8a2",
    description="An Open Source Package for Basic and Advanced Statistical Analysis of Turbulence Data and Other Complex Systems",
    author="Turbulence, Wind energy and Stochastics Group",
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
 
