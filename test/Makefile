GTEST_DIR = ../gtest

CPPFLAGS += -I$(GTEST_DIR)/include -I../include -L../lib

CXXFLAGS += -g -Wall -Wextra

TESTS = geometry_test

GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

all : $(TESTS)

clean :
	rm -f $(TESTS)

geometry_test : geometry_test.cc lib/gtest_main.a ../lib/libgmml.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -lpthread $^ -o $@ -lgmml
