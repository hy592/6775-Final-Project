### You can insert your own directives here ###

set_directive_top -name ahc_top ahc_top

# Set array partition directives for AHC
set_directive_array_partition -type complete -dim 0 AHC x
set_directive_array_partition -type complete -dim 0 AHC J
set_directive_array_partition -type complete -dim 0 AHC MVM_out
set_directive_array_partition -type complete -dim 0 AHC e 

# Set unroll directives for AHC initialization functions
set_directive_unroll AHC::AHC/initialize_vectors
set_directive_unroll AHC::AHC/initialize_matrix

# square_single
set_directive_bind_op square_single tmp_xx -op mul -impl fabric

# AHC::square()
set_directive_inline -off AHC::square
set_directive_pipeline AHC::square
set_directive_unroll -factor 16 AHC::square/square_loop

# AHC::setSpins
set_directive_unroll AHC::setSpins/setSpins_loop

# AHC::matmul
set_directive_inline -off AHC::matmul
set_directive_pipeline AHC::matmul/MVM_outer
set_directive_latency -min=4 -max=5 AHC::matmul/MVM_outer
set_directive_unroll -factor 8 AHC::matmul/MVM_inner

# AHC::IsingEnergy
set_directive_pipeline AHC::IsingEnergy/MVM_inner 
set_directive_unroll -factor 8 AHC::IsingEnergy/MVM_inner 
set_directive_pipeline AHC::IsingEnergy/IsingEnergy_loop

# AHC::update
set_directive_inline -off AHC::update
set_directive_pipeline AHC::update/update_spin_and_error
# set_directive_unroll -factor 8 AHC::update/update_spin_and_error 
set_directive_latency -min=4 -max=10 AHC::update

# AHC::reset
set_directive_inline -off AHC::reset
set_directive_unroll AHC::reset/reset_MVM

# AHC::updateSpins

# AHC::ahc_solver
set_directive_array_partition -type complete -dim 0 AHC::ahc_solver xx
set_directive_array_partition -type complete -dim 0 AHC::ahc_solver dx
set_directive_array_partition -type complete -dim 0 AHC::ahc_solver dx2
set_directive_array_partition -type complete -dim 0 AHC::ahc_solver de

# ahc_top
# set_directive_interface -mode bram ahc_top bestSpinsOut
