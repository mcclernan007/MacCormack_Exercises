program test
implicit none
integer :: idx
integer :: jdx

!exit and cycle
print *,"exit loop behavior\n"
do idx =1,100
	if(idx>10) then
		exit
	end if 
	print *,idx
end do

print *,"cycle loop behavior\n"
do idx =1,10 
	if(mod(idx,2)==0) then
		cycle
	end if 
	print *,idx
end do

print *,"loop labels"
outer_loop: do idx = 1,10,2
	inner_loop: do jdx = 1,10
	if (jdx+idx > 10) then
		cycle outer_loop
	end if
	print*, "idx: ", idx, "jdx",jdx
	end do inner_loop
end do outer_loop

print *,"parallizable loop"


end program test