
import os
import os.path as op
import logging
logging.basicConfig(level=logging.DEBUG)

from subprocess import call
from optparse import OptionParser
from multiprocessing import Process, Lock


class Jobs (list):
    """
    Runs multiple funcion calls on the SAME computer, using multiprocessing.
    """
    def __init__(self, target, args):

        for x in args:
            self.append(Process(target=target, args=x))

    def run(self):
        for pi in self:
            pi.start()

        for pi in self:
            pi.join()


def sh(cmd, grid=False, infile=None, outfile=None, errfile=None,
        background=False):
    """
    simple wrapper for system calls
    """
    if grid:
        return 0  # A fake retcode
    else:
        if infile:
            cmd += " < {0} ".format(infile)
        if outfile and outfile != "stdout":
            cmd += " > {0} ".format(outfile)
        if errfile:
            cmd += " 2> {0} ".format(errfile)
        if background:
            cmd += " & "

        logging.debug(cmd)
        return call(cmd, shell=True)


def set_outfile(instance, outfile="stdout"):
    """
    Add --outfile options to print out to filename.
    """
    assert isinstance(instance, OptionParser)

    instance.add_option("-o", "--outfile", default=outfile,
            help="Outfile name [default: %default]")


def is_newer_file(a, b):
    """
    Check if the file a is newer than file b
    """
    if not (op.exists(a) and op.exists(b)):
        return False
    am = os.stat(a).st_mtime
    bm = os.stat(b).st_mtime
    return am > bm


def need_update(a, b):
    """
    Check if file a is newer than file b and decide whether or not to update
    file b. Can generalize to two lists.
    """
    if isinstance(a, basestring):
        a = [a]
    if isinstance(b, basestring):
        b = [b]

    return any((not op.exists(x)) for x in b) or \
           any(is_newer_file(x, y) for x in a for y in b)


def depends(func):
    """
    Decorator to perform check on infile and outfile. When infile is not present, issue
    warning, and when outfile is present, skip function calls.
    """
    infile = "infile"
    outfile = "outfile"
    def wrapper(*args, **kwargs):
        assert outfile in kwargs, \
            "You need to specify `outfile=` on function call"
        if infile in kwargs:
            infilename = kwargs[infile]
            if isinstance(infilename, basestring):
                infilename = [infilename]
            for x in infilename:
                assert op.exists(x), \
                    "The specified infile `{0}` does not exist".format(x)

        outfilename = kwargs[outfile]
        if need_update(infilename, outfilename):
            return func(*args, **kwargs)
        else:
            msg = "File `{0}` exists. Computation skipped." \
                .format(outfilename)
            logging.debug(msg)

        if isinstance(outfilename, basestring):
            outfilename = [outfilename]

        for x in outfilename:
            assert op.exists(x), \
                    "Something went wrong, `{0}` not found".format(x)

        return outfilename

    return wrapper


def must_open(filename, mode="r", checkexists=False, skipcheck=False):
    """
    Accepts filename and returns filehandle.

    Checks on multiple files, stdin/stdout/stderr, .gz or .bz2 file.
    """
    if isinstance(filename, list):
        assert "r" in mode

        import fileinput
        return fileinput.input(filename)

    if filename in ("-", "stdin"):
        assert "r" in mode
        fp = sys.stdin

    elif filename == "stdout":
        assert "w" in mode
        fp = sys.stdout

    elif filename == "stderr":
        assert "w" in mode
        fp = sys.stderr

    elif filename == "tmp" and mode == "w":
        from tempfile import NamedTemporaryFile
        fp = NamedTemporaryFile(delete=False)

    elif filename.endswith(".gz"):
        import gzip
        fp = gzip.open(filename, mode)

    elif filename.endswith(".bz2"):
        import bz2
        fp = bz2.BZ2File(filename, mode)

    else:
        if checkexists:
            assert mode == "w"
            overwrite = (not op.exists(filename)) if skipcheck \
                        else check_exists(filename)
            if overwrite:
                fp = open(filename, "w")
            else:
                logging.debug("File `{0}` already exists. Skipped."\
                        .format(filename))
                return None
        else:
            fp = open(filename, mode)

    return fp


