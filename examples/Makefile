TARGETS = \
  amber_topology_file_ex \
  carbohydrates_ex \
  minimization_ex \
  parameter_set_ex \
  pdb_ex \
  solvation_ex


all: $(TARGETS)

./%: %.cc
	g++ -I../include -L.. -o $@ $< -lgmml

clean:
	rm -f $(TARGETS) *\.top *\.rst
