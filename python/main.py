"""
Main entry point of the fluorescence analysis application.

Usage from the command line:
    cd "d:\\Lyon thèse\\soft\\Fluorescence_treatment"
    python python/main.py

Usage from another Python script:
    import sys
    sys.path.insert(0, r"d:\\Lyon thèse\\soft\\Fluorescence_treatment\\python")
    from gui.main_window import run_fluorescence_app
    run_fluorescence_app()

Or using the individual functions:
    from functions.load_fluo import load_fluo
    from functions.corrected_fluo import corrected_fluo_ls_wl_1

Required dependencies:
    pip install numpy scipy matplotlib PyQt5
"""

import sys
import os

# Add the python/ folder to the Python search path
# so that "from functions.xxx import xxx" imports work
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)


def main():
    from gui.main_window import run_fluorescence_app
    return run_fluorescence_app(sys.argv)


if __name__ == '__main__':
    sys.exit(main())
