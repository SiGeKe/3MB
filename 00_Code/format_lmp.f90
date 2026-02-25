program format_lammps
    implicit none
    integer::i,N,Npt,dim,ios,narg
    integer, allocatable::pt(:)
    real(8), allocatable::x(:),y(:),z(:),sigma(:)
    real(8)::tmp,tmp2,sigma_avg
    real(8)::xmin,xmax,ymin,ymax,zmin,zmax
    real(8)::boxx,boxy,boxz,volume
    real(8)::rho,box,ibox,hbox,L_target
    logical::use_rho
    character(len=256)::infile,outfile
    character(len=512)::line

    ! ===============
    ! Reading Arguments
    ! ===============

    narg = command_argument_count()
    if ( narg < 2 ) then
        write(*,*) "Usage:"
        write(*,*) "format_lammps infile outfile [rho]"
    end if

    call get_command_argument(1,infile)
    call get_command_argument(2,outfile)
    
    use_rho = .false.
    if ( narg == 3 ) then
        call get_command_argument(3,line)
        read(line,*) rho
        if ( rho <= 0.d0 ) then
            write(*,*) "ERROR: rho must be positive."
            stop
        end if
        use_rho = .true.
    else
        rho = 1.d0
    end if

    ! ===============
    ! Opening File
    ! ===============

    open(10,file=trim(infile),action='read',status='old')

    ! ===============
    ! Reading N
    ! ===============

    read(10,*,iostat=ios) N
    if ( ios /= 0 .or. N <= 0 ) then
        write(*,*) "ERROR: Failed reading particle number."
        stop
    end if

    read(10,*)

    allocate(pt(N),x(N),y(N),z(N))

    ! ===============
    ! Reading Dimensionality
    ! ===============

    read(10,'(A)') line
    read(line,*,iostat=ios) i,tmp,tmp,tmp,tmp

    if ( ios == 0 ) then
        dim = 3
    else
        dim = 2
    end if

    rewind(10)
    read(10,*) N
    read(10,*)

    ! ===============
    ! Reading Data
    ! ===============

    sigma_avg = 0.d0

    do i = 1, N
        if ( dim == 3 ) then
            read(10,*) pt(i),x(i),y(i),z(i),tmp
        else
            read(10,*) pt(i),x(i),y(i),tmp
            z(i) = 0.d0
        end if
    end do

    rewind(10)
    read(10,*)
    read(10,*)

    Npt = maxval(pt)
    allocate(sigma(Npt))

    do i = 1, N
        if ( dim == 3 ) then
            read(10,*) tmp,tmp,tmp,tmp,tmp2
        else
            read(10,*) tmp,tmp,tmp,tmp2
        end if
        sigma(pt(i)) = tmp2
        sigma_avg = sigma_avg + tmp2
    end do

    close(10)

    sigma_avg = sigma_avg / real(N,8)

    ! ===============
    ! Determining Box Size
    ! ===============

    if ( .not. use_rho ) then
        rho = 1.d0
    end if

    L_target = (real(N,8)/rho)**(1.d0/real(dim,8))

    box = L_target * sigma_avg
    ibox = 1.d0/box
    hbox = 0.5d0*box

    ! ===============
    ! Periodic Wrapping
    ! ===============

    do i = 1, N
        x(i) = x(i) - box*anint(x(i)*ibox)
        y(i) = y(i) - box*anint(y(i)*ibox)
        if ( dim == 3 ) then
            z(i) = z(i) - box*anint(z(i)*ibox)
        else
            z(i) = 0.d0
        end if
    end do

    ! ===============
    ! Writing LAMMPS data file
    ! ===============

    open(20,file=trim(outfile),status='replace')

    write(20,*) "LAMMPS data file via format_lammps, version 23 February 2026, timestep = 0"
    write(20,*)
    write(20,*) N, " atoms"
    write(20,*) Npt, " atom types"
    write(20,*)
    write(20,'(F20.16,4X,F20.16,10A)') -hbox, hbox, " xlo xhi"
    write(20,'(F20.16,4X,F20.16,10A)') -hbox, hbox, " ylo yhi"

    if ( dim == 3 ) then
        write(20,'(F20.16,4X,F20.16,10A)') -hbox, hbox, " zlo zhi"
    else
        write(20,'(F20.16,4X,F20.16,10A)') -0.5d0, 0.5d0, " zlo zhi"
    end if

    write(20,*) "0.0 0.0 0.0 xy xz yz"
    write(20,*)
    write(20,*) "Masses"
    write(20,*) 

    do i = 1, Npt
        write(20,*) i, "1.0"
    end do

    write(20,*)
    write(20,*) "Atoms # atomic"
    write(20,*)

    do i = 1, N
        write(20,'(I8,2X,I4,3(2X,F20.16))') i,pt(i),x(i),y(i),z(i)
    end do

    write(20,*)
    write(20,*) "Velocities"
    write(20,*)

    do i = 1, N
        write(20,*) i, "0.0 0.0 0.0"
    end do

    close(20)

end program format_lammps