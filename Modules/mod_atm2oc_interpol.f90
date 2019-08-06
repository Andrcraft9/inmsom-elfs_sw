module atm2oc_interpol
    use mpi_parallel_tools
    implicit none

    !parameters for matrix interpolation from atm to ocean grid
    type(block3D_real8), dimension(:), pointer :: wght_mtrx_a2o !nonzero matrix elements for a2o interpolation

    type(block3D_real4), dimension(:), pointer :: i_input_a2o,   &    !i-numbers of matrix elements
                                                  j_input_a2o         !j-numbers of matrix elements

endmodule atm2oc_interpol
