#!/bin/bash

# Downloads and unpacks libfmt (https://github.com/fmtlib/fmt),
# a modern c++ formatting library.
if [ ! -f fmt-6.2.1.zip ]; then
  wget https://github.com/fmtlib/fmt/releases/download/6.2.1/fmt-6.2.1.zip
fi
unzip -oq fmt-*.*.*.zip && rm fmt-*.*.*.zip
