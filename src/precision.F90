module precision

    use, intrinsic :: iso_fortran_env, only: real32, real64, real128, &
                                             int32, int64

    implicit none

#if INT32
    integer(int32), parameter :: ip = int32
#elif INT64
    integer(int64), parameter :: ip = int64
#else
    integer(int32), parameter :: ip = int32
#endif

#ifdef REAL32
    integer(ip), parameter :: rp = real32
#elif REAL64
    integer(ip), parameter :: rp = real64
#elif REAL128
    integer(ip), parameter :: rp = real128
#else
    integer(ip), parameter :: rp = real64
#endif

    integer(ip), parameter :: lp = int32

end module precision
