import os
import re

import numpy as np
import pandas as pd


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def check_header(path, n=None, debug=False, more_than_one_row=False):
    """
    This monstrosity returns what line a header is at
    """
    try:
        # attempt importing the datafile with the header "n"
        df = pd.read_csv(path, delimiter=r"\t+|\s+", header=n, engine="python")
    except pd.errors.ParserError as pe:
        if debug:
            print("ParserError:", pe)

        # if fail, recursively try again with the next row as the header
        n = -1 if n is None else n
        return check_header(path, n=n + 1, more_than_one_row=True)
    except pd.errors.EmptyDataError:

        if more_than_one_row:
            return None
        else:
            print(os.path.split(path)[-1], "is empty?")
            return -1

    h = df.columns

    # todo: if header is Int64Index, check the 2nd row (i.e. first row of data for the not isfloat)
    # ... so maybe change the h in [not isfloat(x) for x in h] to the second row???
    if isinstance(h, pd.Index) or sum(isfloat(x) for x in h) >= 0.3 * len(h) // 1:
        if debug:
            print("Some are floats...")

        # recursively try again with the next row as the header
        n = -1 if n is None else n
        return check_header(path, n=n + 1, more_than_one_row=True)
    else:
        return n  # <-- the final stop in our recursion journey


def check_datatype(filename):
    """
    Given a filename, try and guess what dataset the data comes from
    (e.g., Si, Kann, Oates, etc.)

    For example, a file named '*_Oates.txt' will be interpreted as data in the same format
    as Sam Oates' data.

    """
    check = lambda x, *args: any([f in filename.lower() for f in [x, *args]])

    if check("zaninoni"):
        datatype = "zaninoni"

    elif check("si", "gendre", "tarot"):
        datatype = "si"

    elif check("liang"):
        datatype = "liang"

    # regex here is for 'GRBid.txt' files
    elif check("kann") or re.search(r"(?<!.)\d+[A-Z]?\.txt", filename):
        datatype = "kann"

    elif re.search(r"combined(?!rest)", filename) or re.search(
        r"comb(?!ined)", filename
    ):
        datatype = "combined"

    elif check("combinedrest"):
        datatype = "combinedrest"

    elif check("oates"):
        datatype = "oates"

    elif check("wczytywanie") or check("block"):
        datatype = "wczytywanie"
    
    elif check("converted"):
        datatype = "GCN"
    
    else:
        datatype = "si"

    return datatype


def read_data(path, datatype="", debug=False):
    data = {}

    if debug:
        print("First 10 Lines:\n", "".join(open(path).readlines()[:10]))

    header = check_header(path)

    if header == -1:
        return

    df = pd.read_csv(path, delimiter=r"\t+|\s+", header=header, engine="python")
    header = h = df.columns

    filename = os.path.split(path)[-1].lower()
    datatype = datatype.lower() if datatype else check_datatype(filename)
    if datatype in ["si", "liang", "combinedrest"]:

        time = df[h[0]]
        flux = df[h[1]]
        fluxerr = df[h[2]]

    elif datatype == "zaninoni":

        time = df[h[0]]
        flux = df[h[2]]
        fluxerr = df[h[3]]

    elif datatype == "kann":

        time = df[h[0]]
        flux = df[h[1]]
        posfluxerr, negfluxerr = df[h[2]], df[h[3]]
        fluxerr = (posfluxerr + negfluxerr) / 2

    elif datatype == "oates":

        time = df[h[0]]
        flux = df[h[2]]
        maxflux = df[h[3]]
        minflux = df[h[4]]
        fluxerr = (maxflux - minflux) / (2 * 1.65)

    elif datatype in ["combined", "comb"]:
        z = df[h[3]]
        beta = df[h[4]]
        time = df[h[0]] * (1 + z)
        flux = df[h[1]] * (1 + z) ** (1 - beta)
        fluxerr = df[h[2]] * (1 + z) ** (1 - beta)

    elif datatype == "wczytywanie":
        time = df[h[0]]
        flux = df[h[3]]
        maxflux = df[h[4]]
        minflux = df[h[5]]
        fluxerr = (maxflux - minflux) / (2 * 1.65)

    elif datatype == 'GCN':
        time = df[h[0]]
        flux = df[h[1]]
        fluxerr = df[h[2]]
        band = df[h[3]]
        source = df[h[4]]
        
    else:
        # if debug:
        # print('No datatype found. Assuming format:\n| time | flux | fluxerr |')
        # read_data(path, datatype='si', debug=debug)
        time = np.array([1])
        flux = np.array([1])
        fluxerr = np.array([0])

    try:
        logtime = np.log10(time)
    except Exception as e:
        print("Issue with logT calculation:", time, e)

    logflux = np.log10(flux)
    logfluxerr = fluxerr / (flux * np.log(10))

    if all(logtime > 0):
        data["time_sec"] = time
        data["flux"] = flux
        data["flux_err"] = fluxerr
        try :
            data["band"] = band
        except NameError:
            data["band"] = ["r" for _ in logtime]
    else:
        raise ImportError("Some logT's are < 0... Ahh!")

    return pd.DataFrame(data)


def readin(directory="."):
    import glob2

    return np.asarray(glob2.glob(directory + "/*.txt"))
