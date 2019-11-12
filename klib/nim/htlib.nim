import zlib

type
  kstring_t = object
    ll: int
    m: int
    s: ptr char
  kstream_t = object
    begin: int
    endd: int
    is_eof: int
  kseq_t = object
    name: kstring_t
    comment: kstring_t
    seq: kstring_t
    qual: kstring_t
    last_char: int
    f: ptr kstream_t
  gzFile = pointer

proc kseq_init(fp: gzFile): ptr kseq_t {.header: "kseq.h",
    importc: "kseq_init".}
proc kseq_rewind(seq: ptr kseq_t) {.header: "kseq.h", importc: "kseq_rewind".}
proc kseq_read(seq: ptr kseq_t): int {.header: "kseq.h", importc: "kseq_read".}


proc main() =
  let fp = gzopen("./test.fasta.gz", "r")
  let seq = kseq_init(fp)
  var ll: int;
  while true:
    if kseq_read(seq) < 0:
      break
    echo seq.name.s
    echo seq.seq.s
  echo "Hello world!"

when isMainModule:
  main()
