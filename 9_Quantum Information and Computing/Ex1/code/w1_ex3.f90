!Week 1 - Exercise 3, Faorlin Tommaso (2021857)

module my_mat_mul
    
    implicit none

    integer*4 ii, jj, kk
    
    contains

    !matrix multiplication by row for two generic matrices
    function matmul_byrow(mat_A, mat_B) result(mat_C1)

        real*4, dimension(:,:) :: mat_A, mat_B
        real*4, dimension(size(mat_A,1),size(mat_B,2)) :: mat_C1

        if(size(mat_A,2).NE.size(mat_B,1)) then
            print *, 'Mismatch in matrix dimensions, unable to carry out the multiplication'
        else
            do ii = 1, size(mat_A,1)
                do kk = 1, size(mat_B,2)
                    do jj = 1, size(mat_A,2)
                        mat_C1(ii,jj)=mat_C1(ii,jj)+mat_A(ii,kk)*mat_B(kk,jj)
                    enddo
                enddo
            enddo 
        endif 
        
    end function matmul_byrow

    !matrix multiplication by columns for two generic matrices
    function matmul_bycol(mat_A, mat_B) result(mat_C2)
 
        real*4, dimension(:,:) :: mat_A, mat_B
        real*4, dimension(size(mat_A,1),size(mat_B,2)) :: mat_C2

        if(size(mat_A,2).NE.size(mat_B,1)) then
            print *, 'Mismatch in matrix dimensions, unable to carry out the multiplication'
        else
            do jj = 1, size(mat_A,2)
                do kk = 1, size(mat_B,2)
                    do ii = 1, size(mat_A,1)
                        mat_C2(ii,jj)=mat_C2(ii,jj)+mat_A(ii,kk)*mat_B(kk,jj)
                    enddo
                enddo
            enddo 
        endif
        
    end function matmul_bycol

end module my_mat_mul

program test_performance_squared
    use my_mat_mul
    implicit none

    integer*4 nn, iter_dim
    real*4, dimension(:,:), allocatable :: mat_A, mat_B, mat_C1, mat_C2, mat_C3  
    integer*4, dimension(:), allocatable :: dimensions_vec   
    real*4 t_i, t_f

    dimensions_vec = (/ 100,200,300,400,500,600,700,800,900,1000,11000/)

    do iter_dim=1, size(dimensions_vec)

        if (dimensions_vec(iter_dim).GT.10000) then
            print *, "Matrices are too big, reduce the size."
        else
            nn = dimensions_vec(iter_dim)
            
            allocate(mat_A(nn, nn))
            allocate(mat_B(nn, nn))
            allocate(mat_C1(nn, nn))
            allocate(mat_C2(nn, nn))
            allocate(mat_C3(nn, nn))
            
            call random_number(mat_A)
            call random_number(mat_B)

            call cpu_time(t_i)
            mat_C1 = matmul_bycol(mat_A, mat_B)
            call cpu_time(t_f)
            print *, "Time for moltiplication by columns", nn, t_f-t_i
            open (10, file="matmul_bycols.txt", position="append", action="write")
            write (10,'(I4","F10.7)') nn, t_f-t_i
            close (10)

            call cpu_time(t_i)
            mat_C2 = matmul_byrow(mat_A, mat_B)
            call cpu_time(t_f)
            print *, "Time for moltiplication by rows", nn, t_f-t_i
            open (10, file="matmul_byrows.txt", position="append", action="write")
            write (10,'(I4","F10.7)') nn, t_f-t_i
            close (10)

            call cpu_time(t_i)
            mat_C3 = matmul(mat_A, mat_B)
            call cpu_time(t_f)
            print *, "Time for moltiplication using 'MATMUL' subroutine", nn, t_f-t_i
            open (10, file="matmul_fortran.txt", position="append", action="write")
            write (10,'(I4","F10.7)') nn, t_f-t_i
            close (10)

            deallocate(mat_A)
            deallocate(mat_B)
            deallocate(mat_C1)
            deallocate(mat_C2)
            deallocate(mat_C3)
        
        end if 
    enddo
end program test_performance_squared