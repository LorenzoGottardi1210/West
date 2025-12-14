#!/bin/bash

${WGET} http://www.quantum-simulation.org/potentials/sg15_oncv/upf/Pb_ONCV_PBE_FR-1.0.upf

cat > pw.in << EOF
&control
calculation  = 'scf'
restart_mode = 'from_scratch'
pseudo_dir   = './'
outdir       = './'
prefix       = 'test'
/
&system
ibrav     = 1
celldm(1) = 20
nat       = 1
ntyp      = 1
ecutwfc   = 25
nbnd      = 20
noncolin  = .true.
lspinorb  = .true.
/
&electrons
diago_full_acc = .true.
/
ATOMIC_SPECIES
Pb  207.2  Pb_ONCV_PBE_FR-1.0.upf
ATOMIC_POSITIONS bohr
Pb  0.0000  0.0000  0.0000
K_POINTS automatic
1 1 1 0 0 0
EOF


cat > wstat.in << EOF
input_west:
  qe_prefix: test
  west_prefix: test
  outdir: ./

wstat_control:
  wstat_calculation: S
  n_pdep_eigen: 20
EOF


cat > wfreq.in << EOF
input_west:
  qe_prefix: test
  west_prefix: test
  outdir: ./

wstat_control:
  wstat_calculation: S
  n_pdep_eigen: 20

wfreq_control:
  wfreq_calculation: XWGQ
  macropol_calculation: N
  n_pdep_eigen_to_use: 20
  qp_bandrange: [13,18]
EOF
