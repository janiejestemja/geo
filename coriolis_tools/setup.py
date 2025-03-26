from setuptools import setup, Extension

# Define your extension module
coriolis_module = Extension(
    'coriolis_module',  # The name of the module
    sources=['coriolis_module.cpp'],  # The C++ source file(s)
)

# Call the setup function
setup(
    name='coriolis_module',
    version='1.0',
    description='A simple C++ extension module',
    ext_modules=[coriolis_module],  # List of extension modules
)