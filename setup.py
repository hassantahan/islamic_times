from setuptools import setup, find_packages, Extension
import numpy
import sys

compile_args = []
link_args = []

if sys.platform == "win32":
    compile_args = ["/Od", "/Zi"]
    link_args = ["/DEBUG"]
else:
    compile_args = ["-O3"]
    link_args = []

astro_core = Extension(
    name="islamic_times.astro_core",
    sources=[
        "src/astro_core.c",
        "src/c_visibilities.c",
        "src/c_moon_equations.c",
        "src/c_sun_equations.c",
        "src/c_time_equations.c",
        "src/c_calculation_equations.c",
        "src/c_datetime.c"
    ],
    include_dirs=[
        "include",
        numpy.get_include()
    ],
    extra_compile_args=compile_args,
    extra_link_args=link_args
)

setup(
    name='islamic_times',
    description='Various calculations for Islamic purposes',
    use_scm_version=True,
    setup_requires=["setuptools_scm"],
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
    license="MIT",
)
