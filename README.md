##### The GLYCAM Molecular Modelling Library (GMML) assists in a variety of molecular modelling tasks, with a particular focus on carbohydrates. It works particularly well with the AMBER suite of molecular simulation programs.

Users of GMML can:

- Read AMBER prep files
- Read AMBER parameter sets
- Read AMBER library (OFF) files
- Read and write PDB files
- Read and write AMBER topology files
- Read and write AMBER coordinate files
- Build residues and molecules using data from input files
- Build residues and molecules manually
- Build oligosaccharides easily with GLYCAM's simple carbohydrate nomenclature
- Generate all possible conformations of a structure that adhere to a user-defined list of possible glycosidic torsions
- Geometrically modify structures in a variety of ways
- Solvate structures
- Minimize structures through an interface with SANDER, AMBER's primary simulation program
- GMML is the engine behind the new and improved GLYCAM web interface, currently being developed.

### Getting Started
To build:
```
    ./configure
    make
```

To run tests:
```
    make check
```

To install:
```
    make install
```

Building requires:
```
    boost headers (for shared_ptr)
```

If you don't have a configure file, run 'sh autogen.sh' first. This requires:
```
    autoconf
    automake
    libtool
```
