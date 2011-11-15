#! /bin/sh

# Run this script to generate the configure script and other files that will
# be included in the distribution.  These files are not checked in because they
# are automatically generated.

set -e

# Check that gtest is present.
if test ! -e gtest; then
  echo "Google Test not present.  Fetching gtest-1.5.0 from the web..."
  curl http://googletest.googlecode.com/files/gtest-1.5.0.tar.bz2 | tar jx
  mv gtest-1.5.0 gtest
fi

set -ex

autoreconf -f -i

rm -rf autom4te.cache config.h.in~
exit 0
