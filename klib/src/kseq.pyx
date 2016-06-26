"""
For wrapping Heng Li's kseq.h FASTA/FASTQ parser.

see http://lh3lh3.users.sourceforge.net/parsefastq.shtml
"""

import sys
cimport kseq
from libc.stdio cimport FILE, fopen, fclose, fread, const_char
from zlib cimport gzFile, gzopen, gzclose, gzrewind


cdef class Kseq:

    def __cinit__(self, char *filename):
        self._c_file = gzopen(filename, 'r')
        self.seq = kseq.kseq_init(self._c_file); # initialize seq

    def __dealloc__(self):
        if self._c_file is not NULL:
            gzclose(self._c_file)

    def read_sequence(self):
        """
        read a single sequence and return a tuple:
        (name, sequence)
        """
        cdef int l
        l = kseq.kseq_read(self.seq)
        if l == -1:
            return None  # EOF
        elif l == -2:
            raise ValueError, "Kseq raised exception for read  %s"%(self.seq.name.s)
        else:
            return (self.seq.name.s, self.seq.seq.s)

    cdef kseq.kseq_t* c_read_sequence(self):
        """
        read a single sequence, and return pointer to active kseq or null if no more seqs to read.
        ctypedef struct  kseq_t:
            kstring_t name
            kstring_t comment
            kstring_t seq
            kstring_t qual
            int last_char

        cdef struct __kstring_t:
            int l               # string length
            int m               # buffer size ?
            char *s
        """
        cdef int l
        l = kseq.kseq_read(self.seq)
        if l < 0:  # NO ERROR CHECKING HERE YET.
            return NULL
        else:
            return self.seq

    def read_sequence_and_quals(self):
        """
        read a single sequence and return a tuple:
        (name, sequence, qualstring)
        """
        cdef int l
        l = kseq.kseq_read(self.seq)
        if l == -1:
            return None  # EOF
        elif l == -2:
            raise ValueError, "Kseq raised exception for read  %s"%(self.seq.name.s)
        else:
            return (self.seq.name.s, self.seq.seq.s, self.seq.qual.s)

    def get_seq_tuples(self):
        """
        returns all sequences as a list of tuples:
        [(name, sequence), (name, sequence), ...]
        """
        self.rewind()
        l = []
        while (kseq.kseq_read(self.seq) >= 0):  # returns length, so as long as >= 0, we're still getting seqs
            l.append((self.seq.name.s, self.seq.seq.s))
        return l

    def get_name_length_tuples(self):
        """
        returns all sequences as a list of tuples:
        [(name, len(sequence)), (name, len(sequence)), ...]
        """
        self.rewind()
        l = []
        while (kseq.kseq_read(self.seq) >= 0):  # returns length, so as long as >= 0, we're still getting seqs
            l.append((self.seq.name.s, self.seq.seq.l))
        return l

    cpdef rewind(self):
        gzrewind(self._c_file)
        kseq.kseq_rewind(self.seq)
