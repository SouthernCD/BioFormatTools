# coding utf8
import setuptools
from bioformattools.versions import get_versions

with open('README.md') as f:
    LONG_DESCRIPTION = f.read()

#
setuptools.setup(
    name="BioFormatTools",
    version=get_versions(),
    author="Yuxing Xu",
    author_email="xuyuxing@mail.kib.ac.cn",
    description="A small set of tools for working with biological file formats.",
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    url="https://github.com/SouthernCD/BioFormatTools",

    entry_points={
        "console_scripts": ["BioFormatTools = bioformattools.cli:main"]
    },  

    packages=setuptools.find_packages(),

    install_requires=[
        "toolbiox>=0.0.38",
    ],

    include_package_data=True,

    python_requires='>=3.5',

)