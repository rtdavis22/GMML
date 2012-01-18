TARGETS = \
  amber_topology_file_ex \
  parameter_set_ex

all: $(TARGETS)

./%: %.cc
	g++ -I../include -L.. -o $@ $< -lgmml

clean:
	rm -rf $(TARGETS)
