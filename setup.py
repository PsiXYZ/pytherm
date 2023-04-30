from setuptools import setup, find_packages
import pytherm

with open("README.rst", "r") as readme_file:
    readme = readme_file.read()

requirements = ['numpy', ]

setup(
    name='pytherm',
    # packages = ['pytherm'],
    packages=find_packages(),
    license='MIT',
    version=pytherm.__version__,
    description='Pytherm is an open-source scientific Python package for thermodynamic modeling.',
    author='Ignaty Efimov',
    install_requires=requirements,
    long_description=readme,
    author_email='efimov.ignaty@gmail.com',
    url='https://github.com/PsiXYZ/pytherm',
    python_requires='>=3.9',
    # package_data={'pytherm': ['activity/db/*',]},
    include_package_data=True,

)
