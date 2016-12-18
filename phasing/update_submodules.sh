#!/bin/bash

# Source: https://github.com/tj/git-extras/pull/80
git pull --recurse-submodules
git submodule update --remote --recursive
