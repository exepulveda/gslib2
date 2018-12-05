module kinds_kit
use iso_c_binding

integer, parameter :: fk = c_float !single precision
integer, parameter :: dk = c_double !single precision

integer, parameter :: sik = c_short !short integer
integer, parameter :: gik = c_int !integer
integer, parameter :: lik = c_long !long integer

!integer, parameter :: fk = selected_read_kind(p=4) !single precision
!integer, parameter :: dk = selected_read_kind(p=8) !single precision

contains

end module
