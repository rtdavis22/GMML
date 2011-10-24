../lib/libgmml.a:
	g++ -g -Wall -Wextra -c -Wno-write-strings -I../include *.cc
	ar rs libgmml.a *.o
	mkdir -p ../lib
	mv libgmml.a ../lib
	rm *.o
