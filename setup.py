from setuptools import setup, find_packages

setup(
    name="geobeam",
    version="3.0.0",
    author='Alexandre JANIN',
    author_email='alexandre.janin@protonmail.com',
    url='https://github.com/AlexandrePFJanin/geoBEAM',
    description='A 3D Boundary-Element Algorithm for Modelling in Geodynamics',
    long_description=open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    license='Apache License 2.0',
    packages=find_packages(),
    package_data={
        "geobeam": ["_all3Ddef*.so"],
    },
    zip_safe=False,
    install_requires=[
        'ipython>=8.15.0',
        'numpy>=1.12',
        'matplotlib>=3.0',
        'termcolor',
        'h5py'
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
)
