# Package

version       = "0.1.0"
author        = "Haibao Tang"
description   = "Heng Li klib in nim"
license       = "MIT"

# Dependencies

requires "nim >= 0.17.2"

task test, "run the tests":
  exec "nim c -r -d:release -l:-lz htlib"
