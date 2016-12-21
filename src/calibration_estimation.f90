module  calibration_estimation_mod

    use global_variables_mod
    use functions_mod
    use grids_mod
    use LinInterpModule
    use omp_lib
    use solve_values_ss_mod
    use policy_to_mm    

    implicit none 

contains      

    subroutine  Main_ss_cali
        integer :: i


        CALL RANDOM_SEED ()   

        sig_e=0.15+0.08*(4-1);
        call make_grids 
        call solve_v2 
        call moment_iterate_ss
        call panel_mm_calculate
        call mm_to_print


    end subroutine  Main_ss_cali

end module  calibration_estimation_mod
