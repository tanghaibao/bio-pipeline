import sys


cdef int BUF_SIZE = 0x10000


cdef class Kopen:
    cdef FILE* fp
    cdef void* x

    def __cinit__(self, char *filename):
        cdef int fd
        self.x = kopen(filename, &fd)
        fp = fdopen(fd, "r")
        if <long>fp == 0:
            print >> sys.stderr, "ERROR: fail to open the input"
            sys.exit(1)
        self.fp = fp

    def read(self):
        cdef int l
        cdef unsigned char buf[0x10000]
        while True:
            l = fread(buf, 1, BUF_SIZE, self.fp)
            if l != 0:
                fwrite(buf, 1, l, stdout)
            else:
                break
            if l == BUF_SIZE:
                break
        fclose(self.fp)
        kclose(self.x)
