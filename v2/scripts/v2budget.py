"""ch4_budget_core for the V2 scripts, made safe to use while a global run is
still writing: HDF5 file locking is off in this process, and the grid areas
come from a year the caller names instead of the newest history year on disk
(the newest one may be the file the model has open; reading it at 01:20 on
2026-09-25 preceded the V-3 abort at 01:24). Callers must only pass years
whose history files are complete."""
import os

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

import sys                                                          # noqa: E402

sys.path.insert(0, '/share/home/dq076/mode/Methane/scripts/figures/budget')
import ch4_budget_core as B                                          # noqa: E402,F401


def load_areas(case, year):
    """(A_land, A_active, A_lake) read from the history file of `year`."""
    B._area_source_year = lambda c: year
    return B.load_areas(case)
