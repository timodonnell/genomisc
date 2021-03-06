# Copyright (c) 2014. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function
import os

from setuptools import setup

current_directory = os.path.dirname(__file__)
readme_filename = 'README.md'
readme_path = os.path.join(current_directory, readme_filename)

readme = ""
try:
    with open(readme_path, 'r') as f:
        readme = f.read()
except IOError as e:
    print(e)
    print("Failed to open %s" % readme_path)

try:
    import pypandoc
    readme = pypandoc.convert(readme, to='rst', format='md')
except ImportError as e:
    print(e)
    print("Failed to convert %s to reStructuredText", readme_filename)
    pass

if __name__ == '__main__':
    setup(
        name='genomisc',
        version="0.0.1",
        description="Collection of scripts for DNA/RNA-seq analysis",
        author="Tim O'Donnell",
        author_email="tim {dot} odonnell {at} mssm {dot} edu",
        url="https://github.com/timodonnell/genomisc",
        license="http://www.apache.org/licenses/LICENSE-2.0.html",
        entry_points={
            'console_scripts': [
                'genomisc-tcga-make-sefara-collection = '
                    'genomisc.tcga.make_sefara_collection:run',
            ]
        },
        package_data={'genomisc': ['data/tcga-code-tables/*.csv']},
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Operating System :: OS Independent',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: Apache Software License',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
        install_requires=[
            "typechecks>=0.0.2",
            "matplotlib>=1.4.3",
            "scipy>=0.15.1",
            "pandas>=0.16.1",
            "lxml>=3.4.4",
        ],
        long_description=readme,
        packages=['genomisc'],            
    )
