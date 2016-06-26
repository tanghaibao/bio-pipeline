cdef extern from *:
    ctypedef char const_char "const char"
    ctypedef void const_void "const void"

cdef extern from "zlib.h":
     ctypedef void *gzFile

     gzFile  *gzopen   (const_char *FILENAME, const_char  *OPENTYPE)
     int     gzclose   (gzFile *STREAM)

     int     gzread    (gzFile gzfile, void *buf, unsigned long len)
     int     gzrewind  (gzFile gzfile)
