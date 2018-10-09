module math_tools
    implicit none

#ifndef __INTEL_COMPILER

    real*8, parameter:: m_dpi = 3.14159265358979d0
    real, parameter :: m_pi = 3.14159

    contains

    function cosd(x)
        real :: x 
        real :: cosd

        cosd = cos((x/180.0)*m_pi)
    end function

    function sind(x)
        real :: x 
        real :: sind

        sind = sin((x/180.0)*m_pi)
    end function
    
    function dcosd(x)
        real*8 :: x 
        real*8 :: dcosd

        dcosd = cos((x/180.0)*m_pi)
    end function

    function dsind(x)
        real*8 :: x 
        real*8 :: dsind

        dsind = sin((x/180.0)*m_pi)
    end function

    function dasind(x)
        real*8 :: x 
        real*8 :: dasind

        dasind = asin((x/180.0)*m_pi)
    end function

    function dacosd(x)
        real*8 :: x 
        real*8 :: dacosd

        dacosd = acos((x/180.0)*m_pi)
    end function

    function dtand(x)
        real*8 :: x 
        real*8 :: dtand

        dtand = tan((x/180.0)*m_pi)
    end function

#endif

end module math_tools