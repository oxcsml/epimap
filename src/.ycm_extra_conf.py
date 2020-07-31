import os


SOURCE_EXTENSIONS = [ '.cpp', '.cxx', '.cc', '.h', '.hpp']


def Settings(**kwargs):
    environ = os.environ["CONDA_PREFIX"]
    return {
        "flags": [
            "-x",
            "c++",
            "-Wall",
            "-Wextra",
            "-Wpedantic",
            "-std=c++17",
            "-isystem", f"{environ}",
            "-I",
            f"{environ}/include",
        ]
    }
