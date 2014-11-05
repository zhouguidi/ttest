program ttest
    use iso_fortran_env
    implicit none
    include "cdi.inc"

    ! local variables
    ! arguments
    integer :: nargin
    real(8), parameter :: NaN = 0.0/0.0

    type Config
        integer :: nf
        integer :: filetype
        integer :: comptype
        integer :: complev
        integer :: chunks
        logical :: silent
        integer :: tType
        integer :: refDate
        integer :: refTime
        integer :: tUnit
        integer :: calendar
        integer :: nvars
        character(1000) :: vlist
        logical :: exclude
        real(8) :: m
        real(8) :: a
    end type Config

    type FileInfoI
        character(100) :: fileName
        integer :: streamID
        integer :: vlistID
        integer :: taxisID
        integer :: nvars
        integer :: nt
    end type FileInfoI
    
    type FileInfoO
        character(100) :: fileName
        integer :: streamID
        integer :: vlistID
        integer :: taxisID
    end type FileInfoO

    ! configuration and IDs
    character(*), parameter :: VERSION = '1.1.0', VERSION_DATE = '2014-05-05'
    type(Config) :: conf
    type(FileInfoI), dimension(:), allocatable :: fis
    type(FileInfoI) :: fi
    type(FileInfoO) :: fo

    ! temps
    integer :: gridID, nx, ny, ivar, ilev, ifi, itime, stat, nmiss, vdate, vtime
    integer :: chunklen, ichk, tsstart, tsend, n, df, nz, ind
    integer, dimension(:), allocatable :: varIDsIn, varIDsOut
    real(8), dimension(:, :), allocatable :: mvar, stdvar
    logical, dimension(:, :), allocatable :: mask
    real(8), dimension(:, :), allocatable :: t
    real(8) :: timebegin, timeend, miss
    logical :: found, last
    character(1000) :: varname, msg, vlist, temp
    character(:), allocatable :: tok, rem

    call cpu_time(timebegin)
    conf = perseArgs()

    ! loop input files
    do ifi = 1, conf%nf
        ! open the input file
        fis(ifi)%streamID = streamOpenRead(trim(fis(ifi)%fileName))
        if (fis(ifi)%streamID .lt. 0) call errorHandler(cdiStringError(fis(ifi)%streamID))
        fis(ifi)%vlistID = streamInqVlist(fis(ifi)%streamID)
        fis(ifi)%nvars = vlistNvars(fis(ifi)%vlistID)

        ! get t-axis IDs and length
        fis(ifi)%taxisID = vlistInqTaxis(fis(ifi)%vlistID)
        fis(ifi)%nt = 0
        do
            stat = streamInqTimestep(fis(ifi)%streamID, fis(ifi)%nt)
            if (stat .eq. 0) then
                exit
            endif
            fis(ifi)%nt = fis(ifi)%nt + 1
        enddo
    enddo
    fi = fis(1)

    ! loop variables in the first input file, assuming all input files contain the same variables
    vlist = conf%vlist
    if (len_trim(vlist) .eq. 0) then   ! if no varlist specified, process all vars
        allocate(varIDsIn(fi%nvars))
        do ivar = 1, fi%nvars
            varIDsIn(ivar) = ivar - 1
        enddo
    else
        if (conf%exclude) then  ! if vars are going to be excluded, first fill varIDsIn with all the vars
            allocate(varIDsIn(fi%nvars))
            do ivar = 0, fi%nvars - 1
                varIDsIn(ivar + 1) = ivar
            enddo
        endif

        ! loop vars in vlist
        do    
            call strtok(vlist, tok, rem, last = last, deli = ",", allowblk = .false.)
            found = .false.
            ! find the var in vlist in file
            do ivar = 0, fi%nvars - 1
                call vlistInqVarName(fi%vlistID, ivar, varname)
                if (trim(varname) == trim(tok)) then
                    found = .true.
                    ind = ivar
                    exit
                endif
            enddo
            if (found) then   ! the var is present in file
                if (conf%exclude) then  ! the var is to be excluded: remove it from the varIDsIn list
                    do ivar = 1, size(varIDsIn, 1)
                        if (ind == varIDsIn(ivar)) then
                            if (ivar == 1) then
                                varIDsIn = varIDsIn(2 : size(varIDsIn, 1))
                            elseif (ivar == size(varIDsIn, 1)) then
                                varIDsIn = varIDsIn(1 : size(varIDsIn, 1) - 1)
                            else
                                varIDsIn = [varIDsIn(1 : ivar - 1), varIDsIn(ivar + 1 : size(varIDsIn, 1))]
                            endif
                        endif
                    enddo
                else    ! the var is to be included: add it to the varIDsIn list
                    if (allocated(varIDsIn)) then
                        varIDsIn = [varIDsIn, ind]
                    else
                        varIDsIn = [ivar]
                    endif
                endif
            else    ! the var is not in the file
                call errorHandler("variable " // trim(tok) // " not found")
            endif

            vlist = rem
            if (last) exit
        enddo
    endif
    ! output vars
    conf%nvars = size(varIDsIn)
    allocate(varIDsOut(conf%nvars))
    do ivar = 1, conf%nvars
        varIDsOut(ivar) = ivar - 1
    enddo

    ! open output file and initialize
    fo%vlistID = vlistCreate()
    call defOVars(fi, fo, varIDsIn, varIDsOut, conf%nvars) ! define variables
    if (conf%tType .eq. -1) then ! determin taxis from input
        fo%taxisID = taxisCreate(taxisInqType(fi%taxisID))
        if (taxisInqType(fi%taxisID) .eq. TAXIS_RELATIVE) then
            call taxisDefRdate(fo%taxisID, taxisInqRdate(fi%taxisID))
            call taxisDefRtime(fo%taxisID, taxisInqRtime(fi%taxisID))
            call taxisDefCalendar(fo%taxisID, taxisInqCalendar(fi%taxisID))
        endif
    else
        fo%taxisID = taxisCreate(conf%tType)
        call taxisDefRdate(fo%taxisID, conf%refDate)
        call taxisDefRtime(fo%taxisID, conf%refTime)
        call taxisDefCalendar(fo%taxisID, conf%calendar)
        call taxisDefTunit(fo%taxisID, conf%tUnit)
    endif
    call vlistDefTaxis(fo%vlistID, fo%taxisID)

    msg = "Student's-t test version " // trim(VERSION) // " (zhouguidi@gmail.com)"
    stat = vlistDefAttTxt(fo%vlistID, CDI_GLOBAL, "TTEST", len_trim(msg), msg)
    msg = "Helmholtz-Center for Ocean Research"
    stat = vlistDefAttTxt(fo%vlistID, CDI_GLOBAL, "Institute", len_trim(msg), msg)

    fo%streamID = streamOpenWrite(trim(fo%fileName), conf%filetype)
    call assert(fo%streamID .ge. 0, cdiStringError(fo%streamID))
    if (conf%comptype .ne. COMPRESS_NONE) then
        call streamDefCompType(fo%streamID, conf%comptype)
        call streamDefCompLevel(fo%streamID, conf%complev)
    endif
    call streamDefVlist(fo%streamID, fo%vlistID)

    ! compute
    if (conf%nf .eq. 1) then   !single file mode: allowing chunks
        chunklen = fi%nt / conf%chunks
        do ichk = 1, conf%chunks
            tsstart = (ichk - 1) * chunklen
            tsend = ichk * chunklen - 1
            if (ichk .eq. conf%chunks) tsend = fi%nt - 1
            n = tsend - tsstart + 1
            df = n - 1

            stat = streamInqTimestep(fi%streamID, tsend)
            vdate = taxisInqVdate(fi%taxisID)
            vtime = taxisInqVtime(fi%taxisID)
            call taxisDefVdate(fo%taxisID, vdate)
            call taxisDefVtime(fo%taxisID, vtime)
            stat = streamDefTimestep(fo%streamID, ichk - 1)

            do ivar = 1, conf%nvars
                gridID = vlistInqVarGrid(fi%vlistID, varIDsIn(1))
                nx = gridInqXsize(gridID)
                ny = gridInqYsize(gridID)
                allocate(mvar(nx, ny), &
                    stdvar(nx, ny), &
                    mask(nx, ny), &
                    t(nx, ny), &
                    stat = stat)
                call assert(stat .eq. 0, "error allocating memory")

                nz = zaxisInqSize(vlistInqVarZaxis(fi%vlistID, varIDsIn(ivar)))
                do ilev = 0, nz - 1
                    ! return mean and std of the current chunk
                    call getLevDataSingleFile(fi, varIDsIn(ivar), nx, ny, tsstart, tsend, ilev, mvar, stdvar, mask)
                    t = ttest2d(conf%m, conf%a, nx, ny, mvar, stdvar, mask, n, df)
                    call streamWriteVarSlice(fo%streamID, varIDsOut(ivar), ilev, t, 0)
                enddo !ilev
                deallocate(mvar, stdvar, mask, t)
            enddo !ivar
        enddo !ichk
    else   ! multi-file mode: disallowing chunks
        ! define the time stemp of the output file to the last timestep of the last input file
        stat = streamInqTimestep(fis(conf%nf)%streamID, fis(conf%nf)%nt - 1)
        vdate = taxisInqVdate(fis(conf%nf)%taxisID)
        vtime = taxisInqVtime(fis(conf%nf)%taxisID)
        call taxisDefVdate(fo%taxisID, vdate)
        call taxisDefVtime(fo%taxisID, vtime)
        stat = streamDefTimestep(fo%streamID, 0)

        ! loop vars
        do ivar = 1, conf%nvars
            ! get grid and level info from the first input file, assuming all inputs have the same grid and level
            gridID = vlistInqVarGrid(fi%vlistID, varIDsIn(1))
            nx = gridInqXsize(gridID)
            ny = gridInqYsize(gridID)
            allocate(mvar(nx, ny), &
                stdvar(nx, ny), &
                mask(nx, ny), &
                t(nx, ny), &
                stat = stat)
            call assert(stat .eq. 0, "error allocating memory")

            nz = zaxisInqSize(vlistInqVarZaxis(fi%vlistID, varIDsIn(ivar)))
            do ilev = 0, nz - 1
                ! return the mean and std of the whole record spreading several files
                call getLevDataMultiFile(conf%nf, fis, varIDsIn(ivar), nx, ny, ilev, mvar, stdvar, mask, n, df)
                t = ttest2d(conf%m, conf%a, nx, ny, mvar, stdvar, mask, n, df)
                call streamWriteVarSlice(fo%streamID, varIDsOut(ivar), ilev, t, 0)
            enddo !ilev
            deallocate(mvar, stdvar, mask, t)
        enddo !ivar
    endif

    ! finialize
    do ifi = 1, conf%nf
        call streamClose(fis(ifi)%streamID)
    enddo
    call streamClose(fo%streamID)
    call vlistDestroy(fo%vlistID)
    call taxisDestroy(fo%taxisID)
    if (allocated(fis)) deallocate(fis)
    if (allocated(varIDsin)) deallocate(varIDsIn)
    if (allocated(varIDsOut)) deallocate(varIDsOut)
    if (allocated(rem)) deallocate(rem)
    if (allocated(tok)) deallocate(tok)

    ! message
    call cpu_time(timeend)
    if (.not. conf%silent) then
        msg = "ttest: processed"
        write(temp, '(I20)')conf%nvars
        msg = trim(msg) // " " // trim(adjustl(temp)) // " variables over"
        write(temp, '(I20)')n
        msg = trim(msg) // " " // trim(adjustl(temp)) // " timesteps ("
        write(temp, '(F20.3)')timeend - timebegin
        msg = trim(msg) // trim(adjustl(temp)) // "s)"
        write(*, '(A)')trim(msg)
    endif

contains
    subroutine defOVars(fi, fo, varIDsIn, varIDsOut, nvars)
        type(FileInfoI), intent(in) :: fi
        type(FileInfoO), intent(inout) :: fo
        integer, intent(in) :: nvars
        integer, dimension(nvars), intent(in) :: varIDsIn
        integer, dimension(nvars), intent(out) :: varIDsOut
        
        integer :: ivar 
        character(1000) :: names

        do ivar = 1, nvars
            varIDsOut(ivar) = vlistDefVar(fo%vlistID, vlistInqVarGrid(fi%vlistID, varIDsIn(ivar)), &
                vlistInqVarZaxis(fi%vlistID, varIDsIn(ivar)), TIME_VARIABLE)
            call vlistDefVarCode(fo%vlistID, varIDsOut(ivar), vlistInqVarCode(fi%vlistID, varIDsIn(ivar)))
            call vlistInqVarName(fi%vlistID, varIDsIn(ivar), names)
            call vlistDefVarName(fo%vlistID, varIDsOut(ivar), trim(adjustl(names)))
            call vlistInqVarLongName(fi%vlistID, varIDsIn(ivar), names)
            call vlistDefVarLongName(fo%vlistID, varIDsOut(ivar), trim(adjustl(names)))
            call vlistInqVarStdName(fi%vlistID, varIDsIn(ivar), names)
            call vlistDefVarStdName(fo%vlistID, varIDsOut(ivar), trim(adjustl(names)))
            call vlistDefVarUnits(fo%vlistID, varIDsOut(ivar), "")
            call vlistDefVarDataType(fo%vlistID, varIDsOut(ivar), DATATYPE_FLT32)
        enddo
    end subroutine defOVars

    subroutine getLevDataSingleFile(fi, varID, nx, ny, tsfirst, tslast, ilev, mvar, stdvar, mask)
        type(FileInfoI), intent(in) :: fi
        integer, intent(in) :: nx, ny, varID, tsfirst, tslast, ilev
        real(8), dimension(nx, ny), intent(out) :: mvar, stdvar
        logical, dimension(nx, ny), intent(out) :: mask

        real(8) :: miss
        integer :: nmiss, itime, n
        real(8), dimension(nx, ny) :: var, var2, varsum

        miss = vlistInqVarMissval(fi%vlistID, varID)
        n = tslast - tsfirst + 1

        var2 = 0
        varsum = 0
        mvar = 0
        mask = .true.
        do itime = tsfirst, tslast
            stat = streamInqTimestep(fi%streamID, itime)
            call streamReadVarSlice(fi%streamID, varID, ilev, var, nmiss)
            varsum = varsum + var   ! varsum = sum_i(var_i)
            var2 = var2 + var ** 2  ! var2 = sum_i(var_i ** 2)
        enddo
        mvar = varsum / n           ! mvar = varsum / n
        ! n * std**2 = sum_i((var_i - mvar)**2) =  sum_i(var_i ** 2) + n * mvar **2 - 2 * sum_i(var_i) * mvar
        stdvar = sqrt((var2 + n * mvar ** 2 - 2 * varsum * mvar) / n)  

        if (nmiss .ne. 0) then
            where(var .eq. miss)
                mask = .false.
                mvar = miss
                stdvar = miss
            endwhere
        endif
    end subroutine getLevDataSingleFile

    subroutine getLevDataMultiFile(nf, fis, varID, nx, ny, ilev, mvar, stdvar, mask, n, df)
        integer, intent(in) :: nf
        type(FileInfoI), dimension(nf), intent(in) :: fis
        integer, intent(in) :: nx, ny, varID, ilev
        real(8), dimension(nx, ny), intent(out) :: mvar, stdvar
        logical, dimension(nx, ny), intent(out) :: mask
        integer, intent(out) :: n, df

        real(8) :: miss
        integer :: nmiss, ifi, itime
        real(8), dimension(nx, ny) :: var, var2, varsum

        miss = vlistInqVarMissval(fi%vlistID, varID)

        var2 = 0
        varsum = 0
        mvar = 0
        mask = .true.
        n = 0
        do ifi = 1, nf
            do itime = 0, fis(ifi)%nt - 1
                stat = streamInqTimestep(fis(ifi)%streamID, itime)
                call streamReadVarSlice(fis(ifi)%streamID, varID, ilev, var, nmiss)
                varsum = varsum + var   ! varsum = sum_i(var_i)
                var2 = var2 + var ** 2  ! var2 = sum_i(var_i ** 2)
            enddo
            n = n + fis(ifi)%nt
        enddo
        df = n - 1
        mvar = varsum / n           ! mvar = varsum / n
        ! n * std**2 = sum_i((var_i - mvar)**2) =  sum_i(var_i ** 2) + n * mvar **2 - 2 * sum_i(var_i) * mvar
        stdvar = sqrt((var2 + n * mvar ** 2 - 2 * varsum * mvar) / n)  

        if (nmiss .ne. 0) then
            where(var .eq. miss)
                mask = .false.
                mvar = miss
                stdvar = miss
            endwhere
        endif
    end subroutine getLevDataMultiFile

    function ttest2d(m, a, nx, ny, mvar, stdvar, mask, n, df) result(t)
        real(8), intent(in) :: m, a
        integer, intent(in) :: nx, ny, n, df
        real(8), dimension(nx, ny), intent(in) :: mvar, stdvar
        logical, dimension(nx, ny), intent(in) :: mask
        real(8), dimension(nx, ny) :: t

        real(8), dimension(nx, ny) :: ts, p

        t = 0
        ts = 0
        p = 0
        where(mask)
            ts = (mvar - m) * sqrt(n * 1.0_8) / stdvar
            p = betai(df / 2.0_8, 0.5_8, df / (df + ts ** 2))
            where(p .lt. (1 - a))
                t = 1
            elsewhere
                t = 0
            endwhere
        endwhere
    end function ttest2d

    elemental function betai(a, b, x)
        real(8), intent(in) :: x, a, b
        real(8) :: betai

        real(8) :: bt

        if (x .lt. 0 .or. x .gt. 1) then
            betai = NaN
        else
            if (x .eq. 0 .or. x .eq. 1) then
                bt = 0
            else
                bt = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a * log(x) + b * log(1 - x))
            endif
            if (x .lt. (a + 1) / (a + b + 2)) then
                betai = bt * betacf(a, b, x) / a
            else
                betai = 1 - bt * betacf(b, a, 1 - x) / b
            endif
        endif
    end function betai

    elemental function betacf(a, b, x)
        real(8), intent(in) :: a, b, x
        real(8) :: betacf

        integer, parameter :: MAXIT = 100
        real(8), parameter :: EPS = epsilon(x), FPMIN = tiny(x) / EPS
        real(8) :: aa, c, d, del, h, qab, qam, qap
        integer :: m, m2

        qab = a + b
        qap = a + 1
        qam = a - 1
        c = 1
        d = 1 - qab * x / qap
        if (abs(d) .lt. FPMIN) d = FPMIN
        d = 1 / d
        h = d
        do m = 1, MAXIT
            m2 = 2 * m
            aa = m * (b - m) * x / ((qam + m2) * (a + m2))
            d = 1 + aa * d
            if (abs(d) .lt. FPMIN) d = FPMIN
            c = 1 + aa / c
            if (abs(c) .lt. FPMIN) c = FPMIN
            d = 1 / d
            h = h * d * c
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
            d = 1 + aa * d
            if (abs(d) < FPMIN) d = FPMIN
            c = 1 + aa / c
            if (abs(c) < FPMIN) c = FPMIN
            d = 1 / d
            del = d * c
            h = h * del
            if (abs(del - 1) <= EPS) exit
        enddo
        if (m > MAXIT) then
            betacf = NaN
        else
            betacf = h
        endif
    end function betacf

    subroutine strtok(str, tok, rem, last, stat, deli, allowblk)
        character(*), intent(in) :: str
        character(:), allocatable, intent(out) :: tok, rem
        logical, optional, intent(out) :: last
        integer, optional, intent(out) :: stat
        character(1), optional, intent(in) :: deli
        logical, optional, intent(in) :: allowblk

        character(1) :: deli_
        logical :: allowblk_
        integer :: ind

        if (present(stat)) stat = 0
        rem = str ! arg is trimed str
        deli_ = " "
        if (present(deli)) deli_ = deli
        allowblk_ = .true.
        if (present(allowblk)) allowblk_ = allowblk

        if (.not. allocated(tok)) allocate(character(len_trim(str)) :: tok)
        if (.not. allocated(rem)) allocate(character(len_trim(str)) :: rem)

        ind = index(trim(rem), deli_)
        if (ind .eq. 0) then ! no commas found: stop
            tok = rem
            rem = ""
        elseif (ind .eq. 1) then ! comma at pos 1
            if (allowblk_) then ! allowed
                tok = ""
                rem = rem(2 : len_trim(rem))
            else
                if (present(stat)) then
                    stat = 1
                else
                    call errorHandler("blank token detected")
                endif
            endif
        elseif (ind .eq. len_trim(rem)) then ! comma at the end
            if (allowblk_) then
                tok = ""
                rem = ""
            else
                if (present(stat)) then
                    stat = 2
                else
                    call errorHandler("blank token detected")
                endif
            endif
        else !good comma
            tok = rem(1 : ind - 1)
            rem = rem(ind + 1 : len_trim(rem))
        endif

        if (present(last)) then
            if (len_trim(rem) .eq. 0) then
                last = .true.
            else
                last = .false.
            endif
        endif
    end subroutine strtok

    subroutine printUsage
        write(*,'(3A)')"TTEST: Student's t-test (version ", VERSION, ")"
        write(*,'(A)')"Usage:"
        write(*,'(A)')"   ttest [options] INPUT OUTPUT"
        write(*,'(A)')"Options:"
        write(*,'(A)')"   -v[variable list]: varaible list to compute, comma sepereted"
        write(*,'(A)')"   -x[variable list]: varaible list to exclude, comma sepereted"
        write(*,'(A)')"   -c[1..]: number of chunks to divide the input time range into"
        write(*,'(A)')"   -m[mean]: test the records against this mean value (default 0)"
        write(*,'(A)')"   -p[significance level]: significance level (<1, default 0.95)"
        write(*,'(A)')"   -nc: output file type netCDF"
        write(*,'(A)')"   -nc2: output file type netCDF version 2 (64-bit)"
        write(*,'(A)')"   -nc4: output file type netCDF-4 (HDF5)"
        write(*,'(A)')"   -nc4c: output file type netCDF-4 classic"
        write(*,'(A)')"   -z: compress in ZIP format"
        write(*,'(A)')"   -gz: compress in GZIP format"
        write(*,'(A)')"   -bz2: compress in BZIP2 format"
        write(*,'(A)')"   -jpg: compress in JPEG format"
        write(*,'(A)')"   -l[0..9]: compress level (0 for none, 9 for most, default 1)"
        write(*,'(A)')"   -a: absolute time axis"
        write(*,'(A)')"   -r[YYYYMMDD[hhmmss]]: relative time axis with reference date & time"
        write(*,'(A)')"   -u[unit]: time axis unit. supported units are:"
        write(*,'(A)')"       d/day/days"
        write(*,'(A)')"       h/hour/hours"
        write(*,'(A)')"       m/min/minute/minutes"
        write(*,'(A)')"       s/sec/second/seconds"
        write(*,'(A)')"   -C[calendar]: calendar. supported calendars are:"
        write(*,'(A)')"       std/standard"
        write(*,'(A)')"       prol/proleptic/proleptic_gregorian"
        write(*,'(A)')"       360/360_days"
        write(*,'(A)')"       365/365_days"
        write(*,'(A)')"       366/366_days"
        write(*,'(A)')"   -s: silent mode"
        write(*,'(A)')"   -h: print this message and exit (must be the only argument)"
        write(*,'(A)')"   -V: print version and exit (must be the only argument)"
        write(*,'(A)')"Author:"
        write(*,'(A)')"   Guidi Zhou"
        write(*,'(3X, A)')VERSION_DATE
        write(*,'(A)')"   Helmholtz-Center for Ocean Research"
        write(*,'(A)')"   Kiel, Germany"
        write(*,'(A)')"   zhouguidi@gmail.com"
        stop
    end subroutine printUsage

    subroutine printVersion
        write(*,'(5A)')"ttest version ", VERSION, " (", VERSION_DATE, ")"
        stop
    end subroutine printVersion

    function perseArgs() result(args)
        type(Config) :: args

        integer :: iarg, nargin, ind, stat
        character(100) :: arg, tok
        logical :: optcpresent

        ! default values
        args%nf = 0
        args%filetype = FILETYPE_NC4C
        args%comptype = COMPRESS_NONE
        args%complev = 1
        args%chunks = 1
        args%silent = .false.
        args%tType = -1 ! default is to determin from input
        args%refDate = 19700101
        args%refTime = 0
        args%tUnit = TUNIT_DAY
        args%calendar = CALENDAR_PROLEPTIC
        args%vlist = ""
        args%exclude = .false.
        args%m = 0
        args%a = 0.95

        nargin = command_argument_count()

        select case(nargin)
        case (0) ! no args
            call printUsage
            stop
        case (1) ! one arg
            call get_command_argument(1, arg)
            arg = adjustl(arg)
            select case(trim(arg))
            case ('-V')
                call printVersion
                stop
            case ('-h')
                call printUsage
                stop
            case default
                call errorHandler("invalid input argument")
            end select
        case default
            ! the last arg is the output file
            call get_command_argument(nargin, fo%fileName)
            fo%fileName = adjustl(fo%fileName)
            nargin = nargin - 1
            
            ! the options
            do iarg = 1, nargin
                call get_command_argument(iarg, arg)
                arg = adjustl(arg)

                print*,nargin,iarg,trim(arg)
                if (arg(1 : 1) .ne. "-") then   ! not starting with -, it's a file name 
                    print*,'file name'
                    args%nf = args%nf + 1
                    if (.not. allocated(fis)) then
                        allocate(fis(1))
                        fis(1)%filename = arg
                    else
                        fis = (/fis(:), fis(1)/)
                        fis(args%nf)%filename = arg
                    endif
                else
                    print*,'arg'
                    ! starting with -, it's an arg
                    arg = arg(2 : len_trim(arg)) ! now without -
                    select case(arg(1 : 1)) ! swtich options
                    case ("v")
                        if (len_trim(arg) .eq. 1) call errorHandler("invalid input argument") ! only one char (v)
                        args%vlist = arg(2 : len_trim(arg)) ! now the actual variable list
                    case ("x")
                        if (len_trim(arg) .eq. 1) call errorHandler("invalid input argument") ! only one char (v)
                        args%vlist = arg(2 : len_trim(arg)) ! now the actual variable list
                        args%exclude = .true.
                    case ("p")
                        if (len_trim(arg) .eq. 1) call errorHandler("invalid input argument") ! only one char (p)
                        arg = arg(2 : len_trim(arg))
                        read(arg, *, iostat = stat) args%a
                        call assert(stat .eq. 0, "unable to convert string to floating point")
                        call assert(args%a .gt. 0 .and. args%a .lt. 1, "significance level exceeds rang (0..1)")
                    case ("m")
                        if (len_trim(arg) .eq. 1) call errorHandler("invalid input argument") ! only one char (m)
                        arg = arg(2 : len_trim(arg))
                        read(arg, *, iostat = stat) args%m
                        call assert(stat .eq. 0, "unable to convert string to floating point")
                    case ("c")
                        if (len_trim(arg) .eq. 1) call errorHandler("invalid input argument") ! only one char (c)
                        tok = arg(2 : len_trim(arg))
                        read(tok, *, iostat = stat)args%chunks
                        if (stat .ne. 0) call errorHandler("invalid chuank number") ! can't read chunk number
                        optcpresent = .true.
                    case ("n")
                        select case(trim(arg))
                        case ("nc")
                            args%filetype = FILETYPE_NC
                        case ("nc2")
                            args%filetype = FILETYPE_NC2
                        case ("nc4")
                            args%filetype = FILETYPE_NC4
                        case ("nc4c")
                            args%filetype = FILETYPE_NC4C
                        case default
                            call errorHandler("invalid output file format")
                        end select
                    case ("z")
                        select case(trim(arg))
                        case ("z", "zip")
                            args%comptype = COMPRESS_ZIP
                        case default
                            call errorHandler("invalid compress type")
                        end select
                    case ("g")
                        select case(trim(arg))
                        case ("g", "gz", "gzip")
                            args%comptype = COMPRESS_GZIP
                        case default
                            call errorHandler("invalid compress type")
                        end select
                    case ("b")
                        select case(trim(arg))
                        case ("b", "bz", "bz2")
                            args%comptype = COMPRESS_BZIP2
                        case default
                            call errorHandler("invalid compress type")
                        end select
                    case ("j")
                        select case(trim(arg))
                        case ("j", "jpg", "jpeg")
                            args%comptype = COMPRESS_JPEG
                        case default
                            call errorHandler("invalid compress type")
                        end select
                    case ("l")
                        if (len_trim(arg) .eq. 1) call errorHandler("invalid input argument") ! only one char (l)
                        tok = arg(2 : len_trim(arg))
                        read(tok, *, iostat = stat)args%complev
                        if (stat .ne. 0) call errorHandler("invalid compress level") ! can't read compress level
                        if (args%complev .lt. 0 .or. args%complev .gt. 9) &
                            call errorHandler("invalid compress level") ! compress level wrong
                    case ("a")
                        if (len_trim(arg) .eq. 1) then
                            args%tType = TAXIS_ABSOLUTE
                        else
                            call errorHandler("unknown argument")
                        endif
                    case ("r")
                        args%tType = TAXIS_RELATIVE
                        if (len_trim(arg) .eq. 1) then ! default reference date & time
                            args%refDate = 19700101
                            args%refTime = 0
                        else
                            arg = arg(2 : len_trim(arg))
                            if (len_trim(arg) .eq. 8) then ! only ref date is given
                                read(arg, *, iostat = stat)args%refDate
                                if (stat .ne. 0) call errorHandler("invalid reference date")
                            elseif (len_trim(arg) .eq. 14) then ! both ref date and time are givem
                                tok = arg(1 : 8)
                                read(tok, *, iostat = stat)args%refDate
                                if (stat .ne. 0) call errorHandler("invalid reference date")
                                tok = arg(9 : 14)
                                read(tok, *, iostat = stat)args%refTime
                                if (stat .ne. 0) call errorHandler("invalid reference time")
                            else ! incorrect ref date & time
                                call errorHandler("invalid reference date & time")
                            endif
                        endif
                    case ("u")
                        if (len_trim(arg) .eq. 1) call errorHandler("invalid input argument") ! only one char (u)
                        arg = arg(2 : len_trim(arg))
                        select case(trim(arg))
                        case ("d", "day", "days")
                            args%tUnit = TUNIT_DAY
                        case ("h", "hour", "hours")
                            args%tUnit = TUNIT_HOUR
                        case ("m", "min", "minute", "minutes")
                            args%tUnit = TUNIT_MINUTE
                        case ("s", "sec", "second", "seconds")
                            args%tUnit = TUNIT_SECOND
                        case default
                            call errorHandler("invalid time unit")
                        end select
                    case ("C")
                        if (len_trim(arg) .eq. 1) call errorHandler("invalid input argument") ! only one char (C)
                        arg = arg(2 : len_trim(arg))
                        select case(trim(arg))
                        case ("std", "standard")
                            args%calendar = CALENDAR_STANDARD
                        case ("prol", "proleptic", "proleptic_gregorian")
                            args%calendar = CALENDAR_PROLEPTIC
                        case ("360", "360_days")
                            args%calendar = CALENDAR_360DAYS
                        case ("365", "365_days")
                            args%calendar = CALENDAR_365DAYS
                        case ("366", "366_days")
                            args%calendar = CALENDAR_366DAYS
                        case default
                            call errorHandler("invalid calendar")
                        end select
                    case ("s")
                        if (len_trim(arg) .eq. 1) then
                            args%silent = .true.
                        else
                            call errorHandler("unknown argument")
                        endif
                    case ("h")
                        if (len_trim(arg) .eq. 1) then
                            call errorHandler("-h must be used as the only input argument")
                        else
                            call errorHandler("unknown argument")
                        end if
                    case ("V")
                        if (len_trim(arg) .eq. 1) then
                            call errorHandler("-V must be used as the only input argument")
                        else
                            call errorHandler("unknown argument")
                        end if
                    case default
                        call errorHandler("unknown argument")
                    end select
                endif
            enddo
        end select

        if (args%nf .gt. 1 .and. optcpresent) then
            call errorHandler("-c must be used only when there is only one input file")
        endif
        if (args%nf .eq. 0) then
            call errorHandler("no input file specified")
        endif
    end function perseArgs

    subroutine assert(mask, msg)
        ! ASSERT: make sure some condition is satisfied, otherwise abort the
        ! program.
        logical, intent(in) :: mask
        character(*), intent(in) :: msg

        if (.not. mask) call errorHandler(msg)
    end subroutine assert

    subroutine errorHandler(msg)
        character(*), intent(in) :: msg

        write(0, '(A, A)')"ERORR(ttest): ", trim(adjustl(msg))
        stop
    end subroutine errorHandler
end program ttest
