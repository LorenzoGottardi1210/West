#!/bin/bash

${WGET} http://www.quantum-simulation.org/potentials/sg15_oncv/upf/Zn_ONCV_PBE_FR-1.0.upf

cat > pw.in << EOF
&control
calculation  = 'scf'
restart_mode = 'from_scratch'
pseudo_dir   = './'
outdir       = './'
prefix       = 'test'
/
&system
ibrav                     = 1
celldm(1)                 = 20
nat                       = 1
ntyp                      = 1
ecutwfc                   = 25
nbnd                      = 30
assume_isolated           = 'mp'
input_dft                 = 'lda'
nspin                     = 4
noncolin                  = .true
lspinorb                  = .true
starting_magnetization(1) = 1
/
&electrons
diago_full_acc = .true.
/
ATOMIC_SPECIES
Zn 65.38  Zn_ONCV_PBE_FR-1.0.upf
ATOMIC_POSITIONS crystal
Zn        0.500000000   0.500000000   0.500000000
K_POINTS automatic
1 1 1 0 0 0
EOF


cat > wbse.in << EOF
input_west:
  qe_prefix: test
  west_prefix: test
  outdir: ./

wbse_init_control:
  wbse_init_calculation: S
  solver: TDDFT

wbse_control:
  wbse_calculation: D
  n_liouville_eigen: 10
  n_liouville_times: 4
  trev_liouville: 0.00000001
  trev_liouville_rel: 0.000001
  l_pre_shift: True
EOF
