AddTest(
    NAME 2D_ComponentTransport_ConcentrationDiffusionOnly
    PATH Parabolic/ComponentTransport/
    EXECUTABLE ogs
    EXECUTABLE_ARGS ConcentrationDiffusionOnly.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-7 RELTOL 1e-10
    DIFF_DATA
    DiffusionOnly_pcs_0_ts_1_t_1.000000_expected.vtu DiffusionOnly_pcs_0_ts_1_t_1.000000.vtu linear_1_to_0 concentration
    DiffusionOnly_pcs_0_ts_1_t_1.000000_expected.vtu DiffusionOnly_pcs_0_ts_1_t_1.000000.vtu zero pressure
    DiffusionOnly_pcs_0_ts_1_t_1.000000_expected.vtu DiffusionOnly_pcs_0_ts_1_t_1.000000.vtu zero darcy_velocity_x
    DiffusionOnly_pcs_0_ts_1_t_1.000000_expected.vtu DiffusionOnly_pcs_0_ts_1_t_1.000000.vtu zero darcy_velocity_y
    VIS DiffusionOnly_pcs_0_ts_1_t_1.000000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_ConcentrationDiffusionAndStorage
    PATH Parabolic/ComponentTransport/
    EXECUTABLE ogs
    EXECUTABLE_ARGS ConcentrationDiffusionAndStorage.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-7 RELTOL 1e-10
    DIFF_DATA
    DiffusionAndStorage_pcs_0_ts_100_t_0.150000_expected.vtu DiffusionAndStorage_pcs_0_ts_100_t_0.150000.vtu concentration concentration
    DiffusionAndStorage_pcs_0_ts_100_t_0.150000_expected.vtu DiffusionAndStorage_pcs_0_ts_100_t_0.150000.vtu pressure pressure
    DiffusionAndStorage_pcs_0_ts_100_t_0.150000_expected.vtu DiffusionAndStorage_pcs_0_ts_100_t_0.150000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorage_pcs_0_ts_100_t_0.150000_expected.vtu DiffusionAndStorage_pcs_0_ts_100_t_0.150000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorage_pcs_0_ts_134_t_1.500000_expected.vtu DiffusionAndStorage_pcs_0_ts_134_t_1.500000.vtu concentration concentration
    DiffusionAndStorage_pcs_0_ts_134_t_1.500000_expected.vtu DiffusionAndStorage_pcs_0_ts_134_t_1.500000.vtu zero pressure
    DiffusionAndStorage_pcs_0_ts_134_t_1.500000_expected.vtu DiffusionAndStorage_pcs_0_ts_134_t_1.500000.vtu zero darcy_velocity_x
    DiffusionAndStorage_pcs_0_ts_134_t_1.500000_expected.vtu DiffusionAndStorage_pcs_0_ts_134_t_1.500000.vtu zero darcy_velocity_y
    VIS DiffusionAndStorage_pcs_0_ts_134_t_1.500000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_DiffusionAndStorageAndAdvection
    PATH Parabolic/ComponentTransport/
    EXECUTABLE ogs
    EXECUTABLE_ARGS DiffusionAndStorageAndAdvection.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-7 RELTOL 1e-10
    DIFF_DATA
    DiffusionAndStorageAndAdvection_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_100_t_5.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvection_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_200_t_35.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvection_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_300_t_155.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvection_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_400_t_315.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvection_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_500_t_495.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvection_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_600_t_720.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvection_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_672_t_900.000000.vtu concentration concentration
    DiffusionAndStorageAndAdvection_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_100_t_5.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvection_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_200_t_35.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvection_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_300_t_155.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvection_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_400_t_315.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvection_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_500_t_495.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvection_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_600_t_720.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvection_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_672_t_900.000000.vtu pressure pressure
    DiffusionAndStorageAndAdvection_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_100_t_5.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvection_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_200_t_35.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvection_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_300_t_155.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvection_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_400_t_315.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvection_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_500_t_495.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvection_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_600_t_720.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvection_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_672_t_900.000000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvection_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_100_t_5.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvection_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_200_t_35.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvection_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_300_t_155.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvection_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_400_t_315.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvection_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_500_t_495.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvection_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_600_t_720.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvection_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_672_t_900.000000.vtu darcy_velocity_y darcy_velocity_y
    VIS DiffusionAndStorageAndAdvection_pcs_0_ts_672_t_900.000000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_DiffusionAndStorageAndGravityAndDispersionHalf
    PATH Parabolic/ComponentTransport/
    EXECUTABLE ogs
    EXECUTABLE_ARGS DiffusionAndStorageAndGravityAndDispersionHalf.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-5 RELTOL 1e-10
    DIFF_DATA
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_100_t_5.700000.vtu concentration concentration
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_200_t_35.700000.vtu concentration concentration
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_300_t_155.700000.vtu concentration concentration
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_400_t_315.700000.vtu concentration concentration
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_500_t_495.700000.vtu concentration concentration
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_600_t_720.700000.vtu concentration concentration
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_672_t_900.000000.vtu concentration concentration
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_100_t_5.700000.vtu pressure pressure
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_200_t_35.700000.vtu pressure pressure
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_300_t_155.700000.vtu pressure pressure
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_400_t_315.700000.vtu pressure pressure
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_500_t_495.700000.vtu pressure pressure
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_600_t_720.700000.vtu pressure pressure
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_672_t_900.000000.vtu pressure pressure
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_100_t_5.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_200_t_35.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_300_t_155.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_400_t_315.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_500_t_495.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_600_t_720.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_672_t_900.000000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_100_t_5.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_200_t_35.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_300_t_155.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_400_t_315.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_500_t_495.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_600_t_720.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_672_t_900.000000.vtu darcy_velocity_y darcy_velocity_y
    VIS DiffusionAndStorageAndAdvection_pcs_0_ts_672_t_900.000000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_DiffusionAndStorageAndAdvectionAndDispersion
    PATH Parabolic/ComponentTransport/
    EXECUTABLE ogs
    EXECUTABLE_ARGS DiffusionAndStorageAndAdvectionAndDispersion.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-7 RELTOL 1e-10
    DIFF_DATA
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_100_t_5.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_200_t_35.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_300_t_155.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_400_t_315.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_500_t_495.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_600_t_720.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_672_t_900.000000.vtu concentration concentration
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_100_t_5.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_200_t_35.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_300_t_155.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_400_t_315.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_500_t_495.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_600_t_720.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_672_t_900.000000.vtu pressure pressure
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_100_t_5.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_200_t_35.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_300_t_155.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_400_t_315.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_500_t_495.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_600_t_720.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_672_t_900.000000.vtu darcy_velocity_x darcy_velocity_x
DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_100_t_5.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_200_t_35.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_300_t_155.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_400_t_315.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_500_t_495.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_600_t_720.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_672_t_900.000000.vtu darcy_velocity_y darcy_velocity_y
    VIS DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_672_t_900.000000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_DiffusionAndStorageAndAdvectionAndDecay
    PATH Parabolic/ComponentTransport/
    EXECUTABLE ogs
    EXECUTABLE_ARGS DiffusionAndStorageAndAdvectionAndDecay.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-7 RELTOL 1e-10
    DIFF_DATA
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_100_t_5.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_200_t_35.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_300_t_155.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_400_t_315.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_500_t_495.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_600_t_720.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_672_t_900.000000.vtu concentration concentration
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_100_t_5.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_200_t_35.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_300_t_155.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_400_t_315.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_500_t_495.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_600_t_720.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_672_t_900.000000.vtu pressure pressure
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_100_t_5.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_200_t_35.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_300_t_155.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_400_t_315.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_500_t_495.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_600_t_720.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_672_t_900.000000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_100_t_5.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_200_t_35.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_300_t_155.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_400_t_315.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_500_t_495.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_600_t_720.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_672_t_900.000000.vtu darcy_velocity_y darcy_velocity_y
    VIS DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_672_t_900.000000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_DiffusionAndStorageAndAdvectionAndDispersionHalf
    PATH Parabolic/ComponentTransport/
    EXECUTABLE ogs
    EXECUTABLE_ARGS DiffusionAndStorageAndAdvectionAndDispersionHalf.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-7 RELTOL 1e-10
    DIFF_DATA
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_100_t_5.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_200_t_35.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_300_t_155.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_400_t_315.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_500_t_495.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_600_t_720.700000.vtu concentration concentration
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_672_t_900.000000.vtu concentration concentration
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_100_t_5.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_200_t_35.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_300_t_155.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_400_t_315.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_500_t_495.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_600_t_720.700000.vtu pressure pressure
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_672_t_900.000000.vtu pressure pressure
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_100_t_5.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_200_t_35.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_300_t_155.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_400_t_315.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_500_t_495.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_600_t_720.700000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_672_t_900.000000.vtu darcy_velocity_x darcy_velocity_x
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_100_t_5.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_200_t_35.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_300_t_155.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_400_t_315.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_500_t_495.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_600_t_720.700000.vtu darcy_velocity_y darcy_velocity_y
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_672_t_900.000000.vtu darcy_velocity_y darcy_velocity_y
    VIS DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_672_t_900.000000.vtu
)

AddTest(
    NAME LARGE_2D_ComponentTransport_Goswami
    PATH Parabolic/ComponentTransport/goswami/
    EXECUTABLE ogs
    EXECUTABLE_ARGS goswami_input.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-1 RELTOL 1e-5
    DIFF_DATA
    ComponentTransport_pcs_0_ts_100_t_150.000000_expected.vtu ComponentTransport_pcs_0_ts_100_t_150.000000.vtu concentration concentration
    ComponentTransport_pcs_0_ts_150_t_350.000000_expected.vtu ComponentTransport_pcs_0_ts_150_t_350.000000.vtu concentration concentration
    ComponentTransport_pcs_0_ts_200_t_750.000000_expected.vtu ComponentTransport_pcs_0_ts_200_t_750.000000.vtu concentration concentration
    ComponentTransport_pcs_0_ts_250_t_1150.000000_expected.vtu ComponentTransport_pcs_0_ts_250_t_1150.000000.vtu concentration concentration
    ComponentTransport_pcs_0_ts_300_t_1950.000000_expected.vtu ComponentTransport_pcs_0_ts_300_t_1950.000000.vtu concentration concentration
    ComponentTransport_pcs_0_ts_350_t_2750.000000_expected.vtu ComponentTransport_pcs_0_ts_350_t_2750.000000.vtu concentration concentration
    ComponentTransport_pcs_0_ts_366_t_3000.000000_expected.vtu ComponentTransport_pcs_0_ts_366_t_3000.000000.vtu concentration concentration
    ComponentTransport_pcs_0_ts_50_t_50.000000_expected.vtu ComponentTransport_pcs_0_ts_50_t_50.000000.vtu concentration concentration
    ComponentTransport_pcs_0_ts_100_t_150.000000_expected.vtu ComponentTransport_pcs_0_ts_100_t_150.000000.vtu pressure pressure
    ComponentTransport_pcs_0_ts_150_t_350.000000_expected.vtu ComponentTransport_pcs_0_ts_150_t_350.000000.vtu pressure pressure
    ComponentTransport_pcs_0_ts_200_t_750.000000_expected.vtu ComponentTransport_pcs_0_ts_200_t_750.000000.vtu pressure pressure
    ComponentTransport_pcs_0_ts_250_t_1150.000000_expected.vtu ComponentTransport_pcs_0_ts_250_t_1150.000000.vtu pressure pressure
    ComponentTransport_pcs_0_ts_300_t_1950.000000_expected.vtu ComponentTransport_pcs_0_ts_300_t_1950.000000.vtu pressure pressure
    ComponentTransport_pcs_0_ts_350_t_2750.000000_expected.vtu ComponentTransport_pcs_0_ts_350_t_2750.000000.vtu pressure pressure
    ComponentTransport_pcs_0_ts_366_t_3000.000000_expected.vtu ComponentTransport_pcs_0_ts_366_t_3000.000000.vtu pressure pressure
    ComponentTransport_pcs_0_ts_50_t_50.000000_expected.vtu ComponentTransport_pcs_0_ts_50_t_50.000000.vtu pressure pressure
    ComponentTransport_pcs_0_ts_100_t_150.000000_expected.vtu ComponentTransport_pcs_0_ts_100_t_150.000000.vtu darcy_velocity_x darcy_velocity_x
    ComponentTransport_pcs_0_ts_150_t_350.000000_expected.vtu ComponentTransport_pcs_0_ts_150_t_350.000000.vtu darcy_velocity_x darcy_velocity_x
    ComponentTransport_pcs_0_ts_200_t_750.000000_expected.vtu ComponentTransport_pcs_0_ts_200_t_750.000000.vtu darcy_velocity_x darcy_velocity_x
    ComponentTransport_pcs_0_ts_250_t_1150.000000_expected.vtu ComponentTransport_pcs_0_ts_250_t_1150.000000.vtu darcy_velocity_x darcy_velocity_x
    ComponentTransport_pcs_0_ts_300_t_1950.000000_expected.vtu ComponentTransport_pcs_0_ts_300_t_1950.000000.vtu darcy_velocity_x darcy_velocity_x
    ComponentTransport_pcs_0_ts_350_t_2750.000000_expected.vtu ComponentTransport_pcs_0_ts_350_t_2750.000000.vtu darcy_velocity_x darcy_velocity_x
    ComponentTransport_pcs_0_ts_366_t_3000.000000_expected.vtu ComponentTransport_pcs_0_ts_366_t_3000.000000.vtu darcy_velocity_x darcy_velocity_x
    ComponentTransport_pcs_0_ts_50_t_50.000000_expected.vtu ComponentTransport_pcs_0_ts_50_t_50.000000.vtu darcy_velocity_x darcy_velocity_x
    ComponentTransport_pcs_0_ts_100_t_150.000000_expected.vtu ComponentTransport_pcs_0_ts_100_t_150.000000.vtu darcy_velocity_y darcy_velocity_y
    ComponentTransport_pcs_0_ts_150_t_350.000000_expected.vtu ComponentTransport_pcs_0_ts_150_t_350.000000.vtu darcy_velocity_y darcy_velocity_y
    ComponentTransport_pcs_0_ts_200_t_750.000000_expected.vtu ComponentTransport_pcs_0_ts_200_t_750.000000.vtu darcy_velocity_y darcy_velocity_y
    ComponentTransport_pcs_0_ts_250_t_1150.000000_expected.vtu ComponentTransport_pcs_0_ts_250_t_1150.000000.vtu darcy_velocity_y darcy_velocity_y
    ComponentTransport_pcs_0_ts_300_t_1950.000000_expected.vtu ComponentTransport_pcs_0_ts_300_t_1950.000000.vtu darcy_velocity_y darcy_velocity_y
    ComponentTransport_pcs_0_ts_350_t_2750.000000_expected.vtu ComponentTransport_pcs_0_ts_350_t_2750.000000.vtu darcy_velocity_y darcy_velocity_y
    ComponentTransport_pcs_0_ts_366_t_3000.000000_expected.vtu ComponentTransport_pcs_0_ts_366_t_3000.000000.vtu darcy_velocity_y darcy_velocity_y
    ComponentTransport_pcs_0_ts_50_t_50.000000_expected.vtu ComponentTransport_pcs_0_ts_50_t_50.000000.vtu darcy_velocity_y darcy_velocity_y
    VIS ComponentTransport_pcs_0_ts_366_t_3000.000000.vtu
)
