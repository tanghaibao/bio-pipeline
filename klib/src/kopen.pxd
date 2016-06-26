from libc.stdio cimport FILE, fclose, fdopen, fread, fwrite, stdout


cdef extern from "klib/kopen.c":
    void *kopen(const char *FILENAME, int *fd)
    int kclose(void *a)
