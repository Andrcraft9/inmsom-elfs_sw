module math_tools
    implicit none

    ! Math Constants
    real(8), parameter:: lat_extr= 89.99999d0
    real(4), parameter:: Pi = 3.1415926,    & !Pi-constant
                     pip180 = Pi/180.0
    real(8), parameter:: dPi = 3.14159265358979d0,   &  !Pi-constant(double precision)
                     dpip180 = dPi/180.0d0

#ifndef __INTEL_COMPILER

    contains

    function cosd(x)
        real :: x
        real :: cosd

        cosd = cos((x/180.0)*Pi)
    end function

    function sind(x)
        real :: x
        real :: sind

        sind = sin((x/180.0)*Pi)
    end function

    function dcosd(x)
        real*8 :: x
        real*8 :: dcosd

        dcosd = cos((x/180.0)*Pi)
    end function

    function dsind(x)
        real*8 :: x
        real*8 :: dsind

        dsind = sin((x/180.0)*Pi)
    end function

    function dasind(x)
        real*8 :: x
        real*8 :: dasind

        dasind = asin((x/180.0)*Pi)
    end function

    function dacosd(x)
        real*8 :: x
        real*8 :: dacosd

        dacosd = acos((x/180.0)*Pi)
    end function

    function dtand(x)
        real*8 :: x
        real*8 :: dtand

        dtand = tan((x/180.0)*Pi)
    end function

#endif

end module math_tools
