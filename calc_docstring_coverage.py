import glob
import os
import sys
import json
from docstr_coverage import get_docstring_coverage
from pathlib import Path


def make_file_list():

    p = os.path.join(d, "matchms")

    return [f for f in glob.glob(p + "**/**/*.py", recursive=True)]


def make_ignore_names():
    ignore_names = ()
    ignore_names_file = Path(d, ".docstr_coverage")

    if os.path.isfile(ignore_names_file):
        ignore_names = tuple([line.split() for line in open(ignore_names_file).readlines() if ' ' in line])

    return ignore_names


def calc_coverage():
    files = make_file_list()
    ignore_names = make_ignore_names()
    result = get_docstring_coverage(files, ignore_names=ignore_names)
    print(json.dumps(result[0], indent=4))
    return result[1]["coverage"]


if __name__ == "__main__":
    d = os.path.join(os.path.dirname(__file__))
    threshold = 39
    coverage = calc_coverage()
    if coverage > threshold:
        print("Docstring coverage ({0:.0f}%) satisfies the threshold ({1:.0f}%).".format(coverage, threshold))
        sys.exit(0)
    else:
        print("Docstring coverage ({0:.0f}%) smaller than threshold ({1:.0f}%).".format(coverage, threshold))
        sys.exit(1)
