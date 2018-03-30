module mod_blocks
    implicit none

contains
    subroutine allocate_blocks2D(blks)
        implicit none

        type(block2D), dimension(:), pointer :: blks
        integer :: i

        allocate(blks(blocks_count))
        do i = 1, blocks_count
            allocate(blks(i)%values(bbnd_start(i) : bbnd_end(i)))
        enddo

    end subroutine

    subroutine deallocate_blocks2D(blks)
        implicit none

        type(block2D), dimension(:), pointer :: blks
        integer :: i

        do i = 1, blocks_count
            deallocate(blks(i))
        enddo
        deallocate(blks)

    end subroutine

end module
