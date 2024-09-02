from setuptools import setup, find_packages
import os

def readme():
    with open('README.md') as f:
        return f.read()

setup(
    name='chimera',
    version='1.0.0',
    description='A versatile tool for metagenomic classification',
    long_description=readme(),
    long_description_content_type='text/markdown',
    author='Qinzhong Tian',
    author_email='tianqinzhong@qq.com',
    url='https://loadstar822.github.io/',
    packages=find_packages(),
    include_package_data=True,
    scripts=['chimera.py'],
    install_requires=[
        'multitax',
        'pandas'
    ],
    entry_points={
        'console_scripts': [
            'chimera=chimera:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',  # 最低Python版本要求
)
