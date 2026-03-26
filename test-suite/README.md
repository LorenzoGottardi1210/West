# WEST test-suite

## List of tests

1. SiH4 molecule, GW, Gamma only
2. CH4 molecule, GW, Gamma only, JSON input
3. CH4 molecule, spin-polarized GW, Gamma only
4. Si bulk, GW, k-mesh 1x1x2
5. SiH4 molecule, westpp output PDEP, Gamma only
6. SiH4 molecule, GW hybrid ACE, Gamma only
7. Si bulk, GW hybrid ACE, k-mesh 1x1x2
8. SiC, westpp localization factor, Wannier localization, dipole moment, Gamma only
9. O2 molecule, GW fractional occupation, Gamma only
10. SiH4 molecule, GW `l_off_diagonal`, Gamma only
11. SiH4 molecule, GW `qp_bands`, Gamma only
12. SiH4 molecule, QDET, Gamma only
13. SiH4 molecule, QDET verbosity, Gamma only
14. NV- diamond, spin-polarized QDET, Gamma only
15. SiH4 molecule, TDDFT (PBE0) Lanczos, Gamma only
16. Formaldehyde molecule, TDDFT (PBE) forces, geNACs and eeNACs, Gamma only
17. NV- diamond spin-polarized TDDFT (PBE) forces, geNACs and eeNACs, Gamma only
18. Formaldehyde molecule, TDDFT (PBE0) forces and eeNACs, Gamma only
19. NV- diamond spin-polarized TDDFT (DDH) forces and eeNACs, Gamma only
20. O2 molecule, spin-flip TDDFT (LDA) forces, Gamma only
21. O2 molecule, spin-flip TDDFT (PBE) forces, Gamma only
22. O2 molecule, spin-flip TDDFT (PBE0) forces, Gamma only
23. O2 molecule, spin-flip TDDFT (HSE) forces, Gamma only
24. Formaldehyde molecule, westpp excited state decomposition, Gamma only
25. SiH4 molecule, BSE (PBE) Lanczos, Gamma only
26. Formaldehyde molecule, BSE (PBE) forces, Gamma only
27. NV- diamond spin-polarized BSE (PBE) forces, Gamma only
28. Formaldehyde molecule, BSE (PBE0) forces, Gamma only
29. NV- diamond spin-polarized BSE (DDH) forces, Gamma only
30. O2 molecule, spin-flip BSE (PBE) forces, Gamma only
31. O2 molecule, spin-flip BSE (HSE) forces, Gamma only
32. Pb atom, GW, spin orbit
33. NV- diamond spin-flip TDDFT (PBE) forces and eeNACs, Gamma only
34. NV- diamond spin-flip TDDFT (PBE0) forces and eeNACs, Gamma only

## Executing tests

In order to manually run the test, simply run `make` or `make all` in the test-suite.

## Adding new tests

Additional tests can be added with the following procedure:

1) Create a new test directory with:
  - `prepare_inputs.sh`, to generate all inputs (including links to download pseudopotential files)
  - `ref` subdirectory, that contains reference data
  - `Makefile`, that contains instructions about the execution of the codes

We recommend that you copy and adjust an existing test script.

2) Add the new subdirectory name to the list in the `test-suite/Makefile`.

3) Customize the file `main_test.py`.

4) Add a short description of the new test to the list above.
