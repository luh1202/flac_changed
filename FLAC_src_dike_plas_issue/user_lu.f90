subroutine write_stress
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

do i = 1,nx-1
    do j = 1,nz-1
        if((j==2).and.(i==44)) then
            open (unit = 1, file = "s111.txt")
                write (1,*) "Here are the s111 ", s11
            close (1)

            open (unit = 1, file = "s122.txt")
                write (1,*) "Here are the s122 ", s22
            close (1)

            open (unit = 1, file = "s133.txt")
                write (1,*) "Here are the s133 ", s33
            close (1)
        end if
    end do
end do

return
end
