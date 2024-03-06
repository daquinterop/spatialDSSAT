import setuptools
import os 
import pdb

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
if os.path.exists('requirements.txt'):
    with open('requirements.txt', 'r') as f:
        requirements = [i.split('=')[0] for i in f.readlines()]
else:
    with open('/home/diego/Py_DSSATTools/requirements.txt', 'r') as f:
        requirements = [i.split('=')[0] for i in f.readlines()]

setuptools.setup(
    name="spatialDSSAT",
    version="0.0b",
    author="Diego Quintero",
    author_email="daquinterop@gmail.com",
    description="A library to run DSSAT spatial mode using python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    # url="https://github.com/daquinterop/Py_DSSATTools",
    # project_urls={
    #     "Bug Tracker": "https://github.com/daquinterop/Py_DSSATTools/issues",
    #     "Documentation": "https://py-dssattools.readthedocs.io/en/latest/"
    #       },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "."},
    python_requires=">=3.6",
    license='MIT',
    install_requires=requirements,
)