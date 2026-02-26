from setuptools import Extension, setup
import numpy
import sys

compile_args = []
link_args = []

if sys.platform == "win32":
    compile_args = [
        "/O2",
        "/GL",
        "/GS",
        "/guard:cf",
        "/Zc:inline",
        "/Oi",
        "/Gy",
        "/DNDEBUG",
        "/MD"
    ]
    link_args = [
        "/LTCG",
        "/OPT:REF",
        "/OPT:ICF",
        "/DYNAMICBASE",
        "/NXCOMPAT",
        "/DEBUG:NONE"
    ]
else:
    compile_args = [
        "-O3",
        "-flto",
        "-fPIC",
        "-fvisibility=hidden",
        "-fno-strict-aliasing",
        "-fwrapv",
        "-fstack-protector-strong",
        "-D_FORTIFY_SOURCE=3",
        "-DNDEBUG"
    ]
    link_args = [
        "-flto",
        "-Wl,-O3",
    ]

astro_core = Extension(
    name="islamic_times.astro_core",
    sources=[
        "src/native/astro_core.c",
        "src/native/c_event_solver.c",
        "src/native/c_visibilities.c",
        "src/native/c_moon_equations.c",
        "src/native/c_sun_equations.c",
        "src/native/c_time_equations.c",
        "src/native/c_calculation_equations.c",
        "src/native/c_datetime.c"
    ],
    include_dirs=[
        "src/native/include",
        numpy.get_include()
    ],
    define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
    extra_compile_args=compile_args,
    extra_link_args=link_args
)

setup(
    ext_modules=[astro_core],
)
