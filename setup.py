from setuptools import setup, find_packages, Extension

astro_core = Extension(
    name="islamic_times.astro_core",
    sources=[
        "src/astro_core.c",
        "src/_sun_equations.c",
        "src/_time_equations.c",
        "src/_calculation_equations.c"
    ],
    include_dirs=[
        "include",
        "C:/Program Files/Python313/include"
    ],
    library_dirs=["C:/Program Files/Python313/libs"],
    extra_compile_args=["/0x"],
)

setup(
    name='islamic_times',
    version='1.8.0',
    description='Various calculations for Islamic purposes',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Hassan Tahan',
    author_email='contact@hassantahan.com',
    url='https://github.com/hassantahan/islamic_times',
    packages=find_packages(),
    install_requires=[
        'numpy', 'timezonefinder', 'pytz', 'matplotlib', 'shapely', 'geopandas'
    ],
    ext_modules=[astro_core],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    python_requires='>=3.7',
)
