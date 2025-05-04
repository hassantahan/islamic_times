from setuptools import setup, find_packages, Extension
import numpy
import sys

compile_args = []
link_args = []

if sys.platform == "win32":
    compile_args = ["/Od", "/Zi", "/fsanitize=address,undefined"]
    link_args = ["/DEBUG", "/fsanitize=address,undefined"]
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
    define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
    extra_compile_args=compile_args,
    extra_link_args=link_args
)

setup(
    name='islamic_times',
    description='Various calculations for Islamic purposes',
    use_scm_version=True,
    setup_requires=["setuptools_scm"],
    long_description=open('README.md', encoding='utf-8').read(),
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
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
    ],
    python_requires='>=3.10',
    license="MIT",
)
