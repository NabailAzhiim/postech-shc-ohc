PROGRAM shcfinal
 
IMPLICIT NONE

! Natural constants
REAL,PARAMETER:: hbar = 1.05457182E-34, lambda = 0.272114079526 
REAL,PARAMETER:: pi = 3.141592653589793238, echarge = 1.60217663E-19, kboltz = 8.6173303E-5

! Calculation parameters
REAL,PARAMETER:: trueEf = 0.8320*13.6057039763, Efmin = -6.00, Efmax = 4.00, lparam = 3.165E-8
COMPLEX,PARAMETER:: ii = (0.0,1.0), i0 = (0.0,0.0), i1= (1.0,0.0)
COMPLEX,PARAMETER,DIMENSION(1,1):: Lsx = (0.0,0.0), Lsy = (0.0,0.0), Lsz = (0.0,0.0)
INTEGER,PARAMETER:: nkmesh = 3, nefmesh = 500
REAL,PARAMETER:: dkb = 2*pi*SQRT(2.0)/(nkmesh-1)
REAL*16,PARAMETER:: tol1 = 0.001, tol2 = 0.0001

! Calculation arrays
REAL*16,ALLOCATABLE:: Emesh(:,:,:,:), omegamesh(:,:,:,:)
REAL*16:: shcresult(nefmesh)

! Temporary variables
INTEGER:: i, j, k, i2, j2
REAL:: kb1, kb2, kb3, Efvalue
REAL*16:: Enk, Emk
COMPLEX*16,DIMENSION(18,1):: umk, unk
REAL*16:: shcEf, Env, omeganv

! Calculation time variables
INTEGER:: tstart, tfinisih, count_rate
REAL:: telapsed

! Date and time variables
CHARACTER(8):: dates
CHARACTER(10):: times
CHARACTER(30):: filename

! ZGEEV variables
INTEGER,PARAMETER:: N = 18 ! matrix size
INTEGER,PARAMETER:: LDA = N, LDVL = N, LDVR = N
INTEGER,PARAMETER:: LWMAX = 1000
INTEGER:: INFO, LWORK
DOUBLE PRECISION:: RWORK(2*N)
COMPLEX*16:: A(LDA,N), VL(LDVL,N), VR(LDVR,N), W(N), WORK(LWMAX)

! External subroutines
EXTERNAL:: ZGEEV,SLASRT

ALLOCATE(Emesh(nkmesh,nkmesh,nkmesh,N))
ALLOCATE(omegamesh(nkmesh,nkmesh,nkmesh,N))

CALL SYSTEM_CLOCK(count_rate=count_rate)
CALL SYSTEM_CLOCK(tstart)
CALL DATE_AND_TIME(DATE=dates,TIME=times)

WRITE(*,*) "===========================SHC CALCULATION==========================="
WRITE(*,*) "Calculation started"
WRITE(*,*) "N k-mesh : ", nkmesh
WRITE(*,*) "N Ef mesh = ", nefmesh

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(A,W,VL,VR,WORK,LWORK,INFO,kb1,kb2,kb3)
DO i=1,nkmesh
    DO j=1,nkmesh
        DO k=1,nkmesh 
            ! k-point
            kb1 = 2*pi*SQRT(2.0)*(i-1)/(nkmesh-1)
            kb2 = 2*pi*SQRT(2.0)*(j-1)/(nkmesh-1)
            kb3 = 2*pi*SQRT(2.0)*(k-1)/(nkmesh-1)

            A = Htot(kr(1,kb1,kb2,kb3),kr(2,kb1,kb2,kb3),kr(3,kb1,kb2,kb3))

            LWORK = -1
            CALL ZGEEV('Vectors','Vectors',N,A,LDA,W,VL,LDVL,VR,LDVR,WORK,LWORK,RWORK,INFO)

            LWORK = MIN(LWMAX,INT(WORK(1)))
            CALL ZGEEV('Vectors','Vectors',N,A,LDA,W,VL,LDVL,VR,LDVR,WORK,LWORK,RWORK,INFO)

            IF (INFO.GT.0) THEN
                WRITE(*,*) 'The algorithm failed to compute eigenvalues'
                STOP
            END IF

            DO i2=1,N
                Emesh(i,j,k,i2) = REALPART(W(i2)) - trueEf
                omegamesh(i,j,k,i2) = omegan(i2,kb1,kb2,kb3,dkb,W,VR)
            END DO
        END DO
    END DO
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(Efvalue,Env,omeganv) REDUCTION(+:shcEf)
DO i2=1,nefmesh
    ! Fermi energy value
    Efvalue = Efmin + (Efmax-Efmin)*(i2-1)/(nefmesh-1)
    shcEf = 0

    DO i=1,nkmesh
        DO j=1,nkmesh
            DO k=1,nkmesh
                ! summation for all energy bands n
                DO j2=1,18
                    Env = Emesh(i,j,k,j2)
                    omeganv = omegamesh(i,j,k,j2)
                    shcEf = shcEf + 2*fdirac(Env,Efvalue)*omeganv*SQRT(2.0)*(echarge**2/hbar)/lparam
                END DO
            END DO
        END DO
    END DO

    shcresult(i2) = shcEf
END DO
!$OMP END PARALLEL DO

DEALLOCATE(Emesh)
DEALLOCATE(omegamesh)

CALL SYSTEM_CLOCK(tfinisih)
telapsed = REAL(tfinisih - tstart)/REAL(count_rate)
WRITE(*,*) "Calculation finished. Execution time = ",telapsed," s"

filename = "shclog_" // dates(7:8) // "-" // dates(5:6) // "-" // dates(1:4) // "_" // &
            & times(1:2) // "." // times(3:4) // "." // times(5:6) // &
            & ".txt"

OPEN(12,FILE=filename,STATUS="NEW",POSITION="APPEND",ACTION="WRITE")

WRITE(12,*) "=================================================================="
WRITE(12,*) "N k-mesh = ", nkmesh
WRITE(12,*) "N Ef mesh = ", nefmesh
WRITE(12,*) "Calculation time = ", telapsed, " s"
WRITE(12,*) "   Ef               SHC"
DO i2=1,nefmesh
    WRITE(12,*) Efmin + (Efmax-Efmin)*(i2-1)/(nefmesh-1), shcresult(i2)
END DO

CLOSE(12)

STOP

CONTAINS

! L operator for p orbital (3x3)
FUNCTION Lp(in) RESULT(C)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: in
    COMPLEX,DIMENSION(3,3):: C

    IF (in == 1) THEN
        C(1,1) = i0 ; C(1,2) = i0 ; C(1,3) = i0
        C(2,1) = i0 ; C(2,2) = i0 ; C(2,3) = -ii*hbar
        C(3,1) = i0 ; C(3,2) = ii*hbar ; C(3,3) = i0
    ELSEIF (in == 2) THEN
        C(1,1) = i0 ; C(1,2) = i0 ; C(1,3) = ii*hbar
        C(2,1) = i0 ; C(2,2) = i0 ; C(2,3) = i0
        C(3,1) = -ii*hbar ; C(3,2) = i0 ; C(3,3) = i0
    ELSEIF (in == 3) THEN
        C(1,1) = i0 ; C(1,2) = -ii*hbar ; C(1,3) = i0
        C(2,1) = ii*hbar ; C(2,2) = i0 ; C(2,3) = i0
        C(3,1) = i0 ; C(3,2) = i0 ; C(3,3) = i0
    ELSE 
        C = i0
    END IF

END FUNCTION Lp

! L operator for d orbital (5x5)
FUNCTION Ld(in) RESULT(C)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: in
    COMPLEX,DIMENSION(5,5):: C

    IF (in == 1) THEN
        C(1,1) = i0 ; C(1,2) = i0 ; C(1,3) = i0 ; C(1,4) = SQRT(REAL(3))*ii*hbar ; C(1,5) = i0
        C(2,1) = i0 ; C(2,2) = i0 ; C(2,3) = i0 ; C(2,4) = ii*hbar ; C(2,5) = i0
        C(3,1) = i0 ; C(3,2) = i0 ; C(3,3) = i0 ; C(3,4) = i0 ; C(3,5) = -ii*hbar
        C(4,1) = -SQRT(REAL(3))*ii*hbar ; C(4,2) = -ii*hbar ; C(4,3) = i0 ; C(4,4) = i0 ; C(4,5) = i0
        C(5,1) = i0 ; C(5,2) = i0 ; C(5,3) = ii*hbar ; C(5,4) = i0 ; C(5,5) = i0
    ELSEIF (in == 2) THEN
        C(1,1) = i0 ; C(1,2) = i0 ; C(1,3) = i0 ; C(1,4) = i0 ; C(1,5) = -SQRT(REAL(3))*ii*hbar
        C(2,1) = i0 ; C(2,2) = i0 ; C(2,3) = i0 ; C(2,4) = i0 ; C(2,5) = ii*hbar
        C(3,1) = i0 ; C(3,2) = i0 ; C(3,3) = i0 ; C(3,4) = ii*hbar ; C(3,5) = i0
        C(4,1) = i0 ; C(4,2) = i0 ; C(4,3) = -ii*hbar ; C(4,4) = i0 ; C(4,5) = i0
        C(5,1) = SQRT(REAL(3))*ii*hbar ; C(5,2) = -ii*hbar ; C(5,3) = i0 ; C(5,4) = i0 ; C(5,5) = i0
    ELSEIF (in == 3) THEN
        C(1,1) = i0 ; C(1,2) = i0 ; C(1,3) = i0 ; C(1,4) = i0 ; C(1,5) = i0
        C(2,1) = i0 ; C(2,2) = i0 ; C(2,3) = -2*ii*hbar ; C(2,4) = i0 ; C(2,5) = i0
        C(3,1) = i0 ; C(3,2) = 2*ii*hbar ; C(3,3) = i0 ; C(3,4) = i0 ; C(3,5) = i0
        C(4,1) = i0 ; C(4,2) = i0 ; C(4,3) = i0 ; C(4,4) = i0 ; C(4,5) = ii*hbar
        C(5,1) = i0 ; C(5,2) = i0 ; C(5,3) = i0 ; C(5,4) = -ii*hbar ; C(5,5) = i0
    ELSE 
        C = i0
    END IF

END FUNCTION Ld

! total L operator (18x18)
FUNCTION Lop(in) RESULT(C)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: in
    COMPLEX,DIMENSION(18,18):: C

    IF (in == 1) THEN
        C = KronProd(DirSum(Lsy,DirSum(Lp(1),Ld(1))),Identity(2))
    ELSE IF (in == 2) THEN
        C = KronProd(DirSum(Lsy,DirSum(Lp(2),Ld(2))),Identity(2))
    ELSE IF (in == 3) THEN
        C = KronProd(DirSum(Lsz,DirSum(Lp(3),Ld(3))),Identity(2))
    ELSE
        C = i0
    END IF

END FUNCTION Lop

! total S operator (18x18)
FUNCTION Sop(in) RESULT(C)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: in
    COMPLEX,DIMENSION(2,2):: S
    COMPLEX,DIMENSION(18,18):: C

    IF (in == 1) THEN
        S(1,1) = i0 ; S(1,2) = i1*hbar/2
        S(2,1) = i1*hbar/2 ; S(2,2) = i0
        C = KronProd(Identity(9),S)
    ELSE IF (in == 2) THEN
        S(1,1) = i0 ; S(1,2) = -ii*hbar/2
        S(2,1) = ii*hbar/2 ; S(2,2) = i0
        C = KronProd(Identity(9),S)
    ELSE IF (in == 3) THEN
        S(1,1) = i1*hbar/2 ; S(1,2) = 0
        S(2,1) = 0 ; S(2,2) = -i1*hbar/2
        C = KronProd(Identity(9),S)
    ELSE
        S = i0; C = i0
    END IF

END FUNCTION Sop

! tight binding hamiltonian (9x9)
FUNCTION Htb(x,y,z) RESULT(C)

    IMPLICIT NONE
    REAL,INTENT(IN):: x, y, z 
    COMPLEX,DIMENSION(9,9):: C

    C(1,1) = 13.586792047772942*i1 + 2*i1*(- 0.339870485327974*Cos(x) - &
            & 0.339870485327974*Cos(y) + &
            & 0.311026392898218*Cos(x)*Cos(y) - &
            & 6.126920614607415*Cos(x/2.)*Cos(y/2.)*Cos(z/2.) - &
            & 0.339870485327974*Cos(z) + &
            & 0.311026392898218*Cos(x)*Cos(z) + &
            & 0.311026392898218*Cos(y)*Cos(z))
    C(2,1) = (0.,6.249642352430717)*Cos(y/2.)*Cos(z/2.)* &
            &  Sin(x/2.) + (0.,2.177184750287526)*Sin(x) - &
            & (0.,0.10967581520682693)*Cos(y)*Sin(x) - &
            & (0.,0.10967581520682693)*Cos(z)*Sin(x) 
    C(1,2) = CONJG(C(2,1))
    C(2,2) = 22.928468397900524*i1 + 2*i1*(3.0890390307791518*Cos(x) - & 
            & 0.11537636971902399*Cos(y) + &
            & 0.4634102774327779*Cos(x)*Cos(y) + &
            & 5.243638312466019*Cos(x/2.)*Cos(y/2.)*Cos(z/2.) - &
            & 0.11537636971902399*Cos(z) + &
            & 0.4634102774327779*Cos(x)*Cos(z) - &
            & 0.22367777337037198*Cos(y)*Cos(z))
    C(3,1) = (0.,6.249642352430717)*Cos(x/2.)*Cos(z/2.)* &
            &  Sin(y/2.) + (0.,2.177184750287526)*Sin(y) - &
            & (0.,0.10967581520682693)*Cos(x)*Sin(y) - &
            & (0.,0.10967581520682693)*Cos(z)*Sin(y)
    C(1,3) = CONJG(C(3,1))
    C(3,2) = - (5.276836230168193,0.)*Cos(z/2.)*Sin(x/2.)* &
            & Sin(y/2.) - (1.3741761016062994,0.)*Sin(x)*Sin(y)
    C(2,3) = CONJG(C(3,2))
    C(3,3) = 22.928468397900524*i1 + 2*i1*(- 0.11537636971902399*Cos(x) + &
            & 3.0890390307791518*Cos(y) + &
            & 0.4634102774327779*Cos(x)*Cos(y) + &
            & 5.243638312466019*Cos(x/2.)*Cos(y/2.)*Cos(z/2.) - &
            & 0.11537636971902399*Cos(z) - &
            & 0.22367777337037198*Cos(x)*Cos(z) + &
            & 0.4634102774327779*Cos(y)*Cos(z))
    C(4,1) = (0.,6.249642352430717)*Cos(x/2.)*Cos(y/2.)* &
            &  Sin(z/2.) + (0.,2.177184750287526)*Sin(z) - &
            & (0.,0.10967581520682693)*Cos(x)*Sin(z) - &
            & (0.,0.10967581520682693)*Cos(y)*Sin(z)
    C(1,4) = CONJG(C(4,1))
    C(4,2) = - (5.276836230168193,0.)*Cos(y/2.)*Sin(x/2.)* &
            & Sin(z/2.) - (1.3741761016062994,0.)*Sin(x)*Sin(z)
    C(2,4) = CONJG(C(4,2))
    C(4,3) = - (5.276836230168193,0.)*Cos(x/2.)*Sin(y/2.)* &
            & Sin(z/2.) - (1.3741761016062994,0.)*Sin(y)*Sin(z)
    C(3,4) = CONJG(C(4,3))
    C(4,4) = 22.928468397900524*i1 + 2*i1*(- 0.11537636971902399*Cos(x) - &
            & 0.11537636971902399*Cos(y) - &
            & 0.22367777337037198*Cos(x)*Cos(y) + &
            & 5.243638312466019*Cos(x/2.)*Cos(y/2.)*Cos(z/2.) + &
            & 3.0890390307791518*Cos(z) + &
            & 0.4634102774327779*Cos(x)*Cos(z) + & 
            & 0.4634102774327779*Cos(y)*Cos(z))
    C(5,1) = 1.5375806063616628*Cos(x) + &
            & 1.5375806063616628*Cos(y) + &
            & 3.0751612127233257*Cos(x)*Cos(y) - &
            & 3.0751612127233257*Cos(z) - &
            & 1.5375806063616628*Cos(x)*Cos(z) - &
            & 1.5375806063616628*Cos(y)*Cos(z)
    C(1,5) = CONJG(C(5,1))
    C(5,2) = - (0.,0.8794727050280321)*Cos(y/2.)*Cos(z/2.)* &
            &  Sin(x/2.) + (0.,1.43540176949965)*Sin(x) + &
            & (0.,0.0796592763081164)*Cos(y)*Sin(x) - &
            & (0.,0.08848710536066326)*Cos(z)*Sin(x)
    C(2,5) = CONJG(C(5,2))
    C(5,3) = - (0.,0.8794727050280321)*Cos(x/2.)*Cos(z/2.)*Sin(y/2.) + &
            & (2.7755575615628914e-17,0.)*Cos(z/2.)*Sin(x/2.)*Sin(y/2.) + &
            & (0.,1.43540176949965)*Sin(y) + &
            & (0.,0.0796592763081164)*Cos(x)*Sin(y) - &
            & (0.,0.08848710536066326)*Cos(z)*Sin(y)
    C(3,5) = CONJG(C(5,3))
    C(5,4) = (0.,1.7589454100560642)*Cos(x/2.)*Cos(y/2.)* &
            &  Sin(z/2.) - (0.,2.8708035389993)*Sin(z) + &
            & (0.,0.008827829052546872)*Cos(x)*Sin(z) + &
            & (0.,0.008827829052546872)*Cos(y)*Sin(z)
    C(4,5) = CONJG(C(5,4))
    C(5,5) = 11.662809448484358*i1 + 2*i1*(- 0.15609143886810173*Cos(x) - &
            & 0.15609143886810173*Cos(y) + &
            & 0.11129465852613399*Cos(x)*Cos(y) + &
            & 1.9294702332256901*Cos(x/2.)*Cos(y/2.)*Cos(z/2.) - &
            & 0.9949851317868189*Cos(z) - &
            & 0.036055115537194996*Cos(x)*Cos(z) - &
            & 0.036055115537194996*Cos(y)*Cos(z))
    C(6,1) = - 1.1874791785029553*Cos(x) + &
            & 1.1874791785029553*Cos(y) + &
            & 0.23495073248014414*Cos(x)*Cos(z) - &
            & 0.23495073248014414*Cos(y)*Cos(z)
    C(1,6) = CONJG(C(6,1))
    C(6,2) = (0.,1.523291408978588)*Cos(y/2.)*Cos(z/2.)* &
            &  Sin(x/2.) - (0.,2.486188794047664)*Sin(x) + &
            & (0.,0.056184803579637656)*Cos(y)*Sin(x) - &
            & (0.,0.04089455514009384)*Cos(z)*Sin(x)
    C(2,6) = CONJG(C(6,2))
    C(6,3) = - (0.,1.523291408978588)*Cos(x/2.)*Cos(z/2.)* &
            &  Sin(y/2.) + (0.,2.486188794047664)*Sin(y) - &
            & (0.,0.056184803579637656)*Cos(x)*Sin(y) + &
            & (0.,0.04089455514009384)*Cos(z)*Sin(y)
    C(3,6) = CONJG(C(6,3))
    C(6,4) = - (0.,0.09707935871973149)*Cos(x)*Sin(z) + &
            & (0.,0.09707935871973149)*Cos(y)*Sin(z)
    C(4,6) = CONJG(C(6,4))
    C(6,5) = 0.9686709988562012*Cos(x) - &
            & 0.9686709988562012*Cos(y) + &
            & 0.17014486344098706*Cos(x)*Cos(z) - &
            & 0.17014486344098706*Cos(y)*Cos(z)
    C(5,6) = CONJG(C(6,5))
    C(6,6) = 11.662809448484358*i1 + 2*i1*(- 0.7153539008139131*Cos(x) - &
            & 0.7153539008139131*Cos(y) - &
            & 0.085171706891638*Cos(x)*Cos(y) + &
            & 1.9294702332256901*Cos(x/2.)*Cos(y/2.)*Cos(z/2.) + &
            & 0.12353979210480398*Cos(z) + &
            & 0.062178067171691*Cos(x)*Cos(z) + &
            & 0.062178067171691*Cos(y)*Cos(z))
    C(7,1) = (7.101780615869233,0.)*Cos(z/2.)*Sin(x/2.)* &
            & Sin(y/2.) - (0.4699014649602882,0.)*Sin(x)*Sin(y)
    C(1,7) = CONJG(C(7,1))
    C(7,2) = - (0.,6.0723174587185165)*Cos(x/2.)*Cos(z/2.)* &
            &  Sin(y/2.) - (0.,0.20408555964449998)*Sin(y) - &
            & (0.,0.13797391385982533)*Cos(x)*Sin(y) + &
            & (0.,0.056184803579637656)*Cos(z)*Sin(y)
    C(2,7) = CONJG(C(7,2))
    C(7,3) = - (0.,6.0723174587185165)*Cos(y/2.)*Cos(z/2.)* &
            &  Sin(x/2.) - (0.,0.20408555964449998)*Sin(x) - &
            & (0.,0.13797391385982533)*Cos(y)*Sin(x) + &
            & (0.,0.056184803579637656)*Cos(z)*Sin(x)
    C(3,7) = CONJG(C(7,3))
    C(7,4) = (0.,7.595608867697108)*Sin(x/2.)*Sin(y/2.)*Sin(z/2.)
    C(4,7) = CONJG(C(7,4))
    C(7,5) = (2.072112021124339,0.)*Cos(z/2.)*Sin(x/2.)*Sin(y/2.) + &
            & (0.2743055693148323,0.)*Sin(x)*Sin(y)
    C(5,7) = CONJG(C(7,5))
    C(7,6) = 0.
    C(6,7) = CONJG(C(7,6))
    C(7,7) = 13.172090190575318*i1 + 2*i1*(- 0.171976098260432*Cos(x) - &
            & 0.171976098260432*Cos(y) + &
            & 0.269665052810266*Cos(x)*Cos(y) - &
            & 1.4591966340093117*Cos(x/2.)*Cos(y/2.)*Cos(z/2.) + &
            & 0.12353979210480398*Cos(z) - &
            & 0.026531122753785*Cos(x)*Cos(z) - &
            & 0.026531122753785*Cos(y)*Cos(z))
    C(8,1) = (7.101780615869233,0.)*Cos(x/2.)*Sin(y/2.)* &
            & Sin(z/2.) - (0.4699014649602882,0.)*Sin(y)*Sin(z)
    C(1,8) = CONJG(C(8,1))
    C(8,2) = (0.,7.595608867697108)*Sin(x/2.)*Sin(y/2.)*Sin(z/2.)
    C(2,8) = CONJG(C(8,2))
    C(8,3) = - (0.,6.0723174587185165)*Cos(x/2.)*Cos(y/2.)* &
            &  Sin(z/2.) - (0.,0.20408555964449998)*Sin(z) + &
            & (0.,0.056184803579637656)*Cos(x)*Sin(z) - &
            & (0.,0.13797391385982533)*Cos(y)*Sin(z)
    C(3,8) = CONJG(C(8,3))
    C(8,4) = - (0.,6.0723174587185165)*Cos(x/2.)*Cos(z/2.)* &
            &  Sin(y/2.) - (0.,0.20408555964449998)*Sin(y) + &
            & (0.,0.056184803579637656)*Cos(x)*Sin(y) - &
            & (0.,0.13797391385982533)*Cos(z)*Sin(y)
    C(4,8) = CONJG(C(8,4))
    C(8,5) = - (1.0360560105621694,0.)*Cos(x/2.)*Sin(y/2.)* &
            & Sin(z/2.) - (0.13715278465741615,0.)*Sin(y)*Sin(z)
    C(5,8) = CONJG(C(8,5))
    C(8,6) = (1.794501649780795,0.)*Cos(x/2.)*Sin(y/2.)* &
            & Sin(z/2.) - (0.,5.551115123125783e-17)*Sin(x/2.)*Sin(y/2.)* &
            & Sin(z/2.) + (0.23755559142619792,0.)*Sin(y)*Sin(z)
    C(6,8) = CONJG(C(8,6))
    C(8,7) = (4.982832084689213,0.)*Cos(y/2.)*Sin(x/2.)* &
            & Sin(z/2.) + (0.11728116827570598,0.)*Sin(x)*Sin(z)
    C(7,8) = CONJG(C(8,7))
    C(8,8) = 13.172090190575318*i1 + 2*i1*(0.12353979210480398*Cos(x) - &
            & 0.171976098260432*Cos(y) - &
            & 0.026531122753785*Cos(x)*Cos(y) - &
            & 1.4591966340093117*Cos(x/2.)*Cos(y/2.)*Cos(z/2.) - &
            & 0.171976098260432*Cos(z) - &
            & 0.026531122753785*Cos(x)*Cos(z) + &
            & 0.269665052810266*Cos(y)*Cos(z))
    C(9,1) = (7.101780615869233,0.)*Cos(y/2.)*Sin(x/2.)* &
            & Sin(z/2.) - (0.4699014649602882,0.)*Sin(x)*Sin(z)
    C(1,9) = CONJG(C(9,1))
    C(9,2) = - (0.,6.0723174587185165)*Cos(x/2.)*Cos(y/2.)* &
            &  Sin(z/2.) - (0.,0.20408555964449998)*Sin(z) - &
            & (0.,0.13797391385982533)*Cos(x)*Sin(z) + &
            & (0.,0.056184803579637656)*Cos(y)*Sin(z)
    C(2,9) = CONJG(C(9,2))
    C(9,3) = (0.,7.595608867697108)*Sin(x/2.)*Sin(y/2.)*Sin(z/2.)
    C(3,9) = CONJG(C(9,3))
    C(9,4) = - (0.,6.0723174587185165)*Cos(y/2.)*Cos(z/2.)* &
            &  Sin(x/2.) - (0.,0.20408555964449998)*Sin(x) + &
            & (0.,0.056184803579637656)*Cos(y)*Sin(x) - &
            & (0.,0.13797391385982533)*Cos(z)*Sin(x)
    C(4,9) = CONJG(C(9,4))
    C(9,5) = - (1.0360560105621694,0.)*Cos(y/2.)*Sin(x/2.)* &
            & Sin(z/2.) - (0.13715278465741615,0.)*Sin(x)*Sin(z)
    C(5,9) = CONJG(C(9,5))
    C(9,6) = - (1.794501649780795,0.)*Cos(y/2.)*Sin(x/2.)* &
            & Sin(z/2.) - (0.23755559142619792,0.)*Sin(x)*Sin(z)
    C(6,9) = CONJG(C(9,6))
    C(9,7) = (4.982832084689213,0.)*Cos(x/2.)*Sin(y/2.)* &
            & Sin(z/2.) + (0.11728116827570598,0.)*Sin(y)*Sin(z)
    C(7,9) = CONJG(C(9,7))
    C(9,8) = (4.982832084689213,0.)*Cos(z/2.)*Sin(x/2.)* &
            & Sin(y/2.) + (0.11728116827570598,0.)*Sin(x)*Sin(y)
    C(8,9) = CONJG(C(9,8))
    C(9,9) = 13.172090190575318*i1 + 2*i1*(- 0.171976098260432*Cos(x) + &
            & 0.12353979210480398*Cos(y) - &
            & 0.026531122753785*Cos(x)*Cos(y) - &
            & 1.4591966340093117*Cos(x/2.)*Cos(y/2.)*Cos(z/2.) - &
            & 0.171976098260432*Cos(z) + &
            & 0.269665052810266*Cos(x)*Cos(z) - &
            & 0.026531122753785*Cos(y)*Cos(z))

END FUNCTION Htb

! total hamiltonian (18x18)
FUNCTION Htot(x,y,z) RESULT(C)

    IMPLICIT NONE
    REAL,INTENT(IN):: x, y, z
    COMPLEX,DIMENSION(18,18):: Hsoc, C

    Hsoc = (2*lambda)*(MATMUL(Lop(1)/hbar,Sop(1)/hbar)+MATMUL(Lop(2)/hbar,Sop(2)/hbar)+MATMUL(Lop(3)/hbar,Sop(3)/hbar))
    C = Hsoc + KronProd(Htb(x,y,z),Identity(2))
    
END FUNCTION Htot

! k transormation (kb1,kb2,kb3 -> kx,ky,kz)
FUNCTION kr(ind,kb1,kb2,kb3) RESULT (C)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: ind
    REAL,INTENT(IN):: kb1, kb2, kb3
    REAL:: C

    IF (ind == 1) THEN
        C = (kb1 + kb3)/SQRT(2.0)
    ELSE IF (ind == 2) THEN
        C = (kb1 + kb2)/SQRT(2.0)
    ELSE IF (ind == 3) THEN
        C = (kb2 + kb3)/SQRT(2.0)
    ELSE
        C = 0
    END IF

END FUNCTION kr    

! partial derivative of Hamiltonian operator (e.g. Hg(1,-,-,-) = dH/dk1)
FUNCTION Hg(ind,kb1,kb2,kb3,dk) RESULT(C)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: ind
    REAL,INTENT(IN):: kb1, kb2, kb3, dk
    COMPLEX,DIMENSION(18,18):: C

    IF (ind == 1) THEN
        C = (- Htot(kr(1,kb1+2*dk,kb2,kb3),kr(2,kb1+2*dk,kb2,kb3),kr(3,kb1+2*dk,kb2,kb3)) + &
            & 4*Htot(kr(1,kb1+dk,kb2,kb3),kr(2,kb1+dk,kb2,kb3),kr(3,kb1+dk,kb2,kb3)) - &
            & 3*Htot(kr(1,kb1,kb2,kb3),kr(2,kb1,kb2,kb3),kr(3,kb1,kb2,kb3)))/(2*dk)
    ELSE IF (ind == 2) THEN
        C = (- Htot(kr(1,kb1,kb2+2*dk,kb3),kr(2,kb1,kb2+2*dk,kb3),kr(3,kb1,kb2+2*dk,kb3)) + &
            & 4*Htot(kr(1,kb1,kb2+dk,kb3),kr(2,kb1,kb2+dk,kb3),kr(3,kb1,kb2+dk,kb3)) - &
            & 3*Htot(kr(1,kb1,kb2,kb3),kr(2,kb1,kb2,kb3),kr(3,kb1,kb2,kb3)))/(2*dk)
    ELSE IF (ind == 3) THEN
        C = (- Htot(kr(1,kb1,kb2,kb3+2*dk),kr(2,kb1,kb2,kb3+2*dk),kr(3,kb1,kb2,kb3+2*dk)) + &
            & 4*Htot(kr(1,kb1,kb2,kb3+dk),kr(2,kb1,kb2,kb3+dk),kr(3,kb1,kb2,kb3+dk)) - &
            & 3*Htot(kr(1,kb1,kb2,kb3),kr(2,kb1,kb2,kb3),kr(3,kb1,kb2,kb3)))/(2*dk) 
    ELSE
        C = i0
    END IF

END FUNCTION Hg

! spin current operator (e.g. Js(1,2,-,-,-) = (vx.Sy+Sy.vx)/2)
FUNCTION Js(vind,sind,kb1,kb2,kb3,dk) RESULT(C)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: vind, sind
    REAL,INTENT(IN):: kb1, kb2, kb3, dk
    COMPLEX,DIMENSION(18,18):: v, s 
    COMPLEX,DIMENSION(18,18):: C

    IF (vind == 1) THEN
        v = (Hg(1,kb1,kb2,kb3,dk)-Hg(2,kb1,kb2,kb3,dk)+Hg(3,kb1,kb2,kb3,dk))/(hbar*SQRT(2.0))
    ELSE IF (vind == 2) THEN
        v = (Hg(1,kb1,kb2,kb3,dk)+Hg(2,kb1,kb2,kb3,dk)-Hg(3,kb1,kb2,kb3,dk))/(hbar*SQRT(2.0))
    ELSE IF (vind == 3) THEN
        v = (-Hg(1,kb1,kb2,kb3,dk)+Hg(2,kb1,kb2,kb3,dk)+Hg(3,kb1,kb2,kb3,dk))/(hbar*SQRT(2.0))
    ELSE
        v = 0
    END IF

    IF (sind == 1) THEN
        s = Sop(1)
    ELSE IF (sind == 2) THEN
        s = Sop(2)
    ELSE IF (sind == 3) THEN
        s = Sop(3)
    ELSE
        s = 0
    END IF

    C = (MATMUL(v,s)+MATMUL(s,v))/2

END FUNCTION Js

! omega(m,n)
FUNCTION omegamn(m,n,kb1,kb2,kb3,dk,W,VR) RESULT (C)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: m,n
    REAL,INTENT(IN):: kb1, kb2, kb3, dk
    COMPLEX*16,INTENT(IN):: W(18)
    COMPLEX*16,INTENT(IN):: VR(18,18)
    COMPLEX*16,DIMENSION(18,1):: um, un
    COMPLEX*16,DIMENSION(18,1):: matJ, matv
    COMPLEX*16:: C1, C2
    REAL*16:: Em, En, C
    REAL*16,PARAMETER:: eta = 0.01
    INTEGER:: i

    Em = REALPART(W(m)) - trueEf
    En = REALPART(W(n)) - trueEf
    DO i=1,18
        um(i,1) = VR(i,m)
        un(i,1) = VR(i,n)
    END DO

    matJ = MATMUL(Js(2,3,kb1,kb2,kb3,dk),um)
    matv = MATMUL((Hg(1,kb1,kb2,kb3,dk)-Hg(2,kb1,kb2,kb3,dk)+Hg(3,kb1,kb2,kb3,dk))/(SQRT(2.0)),un)

    C1 = 0
    C2 = 0

    IF (ABS(En-Em) < tol2) THEN
        C = 0
    ELSE
        DO i=1,18
            IF (ISNAN(REALPART(CONJG(un(i,1))*matJ(i,1)))) C1 = 0
            IF (ISNAN(IMAGPART(CONJG(un(i,1))*matJ(i,1)))) C1 = 0
            C1 = C1 + CONJG(un(i,1))*matJ(i,1)
            IF (ISNAN(REALPART(CONJG(um(i,1))*matv(i,1)))) C2 = 0
            IF (ISNAN(IMAGPART(CONJG(un(i,1))*matJ(i,1)))) C2 = 0
            C2 = C2 + CONJG(um(i,1))*matv(i,1)
        END DO
        C = IMAGPART(C1*C2/((En*i1 - Em*i1 + ii*eta)**2))*((dkb/(2*pi))**3)
        IF (ABS(C) > tol1) THEN
            C = 0
        ELSE IF (ISNAN(C)) THEN
            C = 0
        END IF
    END IF

END FUNCTION omegamn

! omega(n)
FUNCTION omegan(n,kb1,kb2,kb3,dk,W,VR) RESULT (C)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: n
    REAL,INTENT(IN):: kb1, kb2, kb3, dk
    COMPLEX*16,INTENT(IN):: W(18)
    COMPLEX*16,INTENT(IN):: VR(18,18)
    REAL*16:: C
    INTEGER:: i

    C = 0
    DO i=1,18
        C = C + omegamn(i,n,kb1,kb2,kb3,dk,W,VR)
    END DO

END FUNCTION omegan

! Fermi-Dirac distribution
FUNCTION fdirac(Ek,Efermi) RESULT (C)

    IMPLICIT NONE
    REAL*16,INTENT(IN):: Ek
    REAL,INTENT(IN):: Efermi
    REAL,PARAMETER:: TK = 0.00001
    REAL*16:: C

    C = 1/(EXP((Ek-Efermi)/(kboltz*TK))+1)

END FUNCTION fdirac


!=============================================TOOLS FUNCTIONS========================================================!
! identity matrix
FUNCTION Identity(n) RESULT (C)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: n
    COMPLEX,DIMENSION(:,:),ALLOCATABLE:: C
    INTEGER:: i,j

    ALLOCATE(C(n,n))
    DO i=1,n
        DO j=1,n
            IF (i == j) THEN
                C(i,j) = i1
            ELSE
                C(i,j) = i0
            END IF
        END DO
    END DO

END FUNCTION Identity

! kronecker product
FUNCTION KronProd(A,B) RESULT(C)
    
    IMPLICIT NONE
    COMPLEX,DIMENSION(:,:),INTENT(IN):: A, B
    COMPLEX,DIMENSION(:,:),ALLOCATABLE:: C
    INTEGER:: i,j,m,n,p,q

    ALLOCATE(C(SIZE(A,1)*SIZE(B,1),SIZE(A,2)*SIZE(B,2)))
    C = (0.0,0.0)
    
    DO i = 1,SIZE(A,1)
        DO j = 1,SIZE(A,2)
            n = (i-1)*SIZE(B,1) + 1
            m = n+SIZE(B,1) - 1
            p = (j-1)*SIZE(B,2) + 1
            q = p+SIZE(B,2) - 1
            C(n:m,p:q) = A(i,j)*B
        END DO
    END DO

END FUNCTION KronProd

! direct sum
FUNCTION DirSum(A,B) RESULT (C)
    
    IMPLICIT NONE
    COMPLEX,DIMENSION(:,:),INTENT(IN):: A, B
    COMPLEX,DIMENSION(:,:),ALLOCATABLE:: C
    INTEGER :: p,q

    ALLOCATE(C(SIZE(A,1)+SIZE(B,1),SIZE(A,2)+SIZE(B,2)))
    C = i0

    p = SIZE(A,1) + SIZE(B,1) 
    q = SIZE(A,2) + SIZE(B,2) 

    C(1:SIZE(A,1),1:SIZE(A,2)) = A
    C(SIZE(A,1)+1:p,SIZE(A,2)+1:q) = B

END FUNCTION DirSum

END PROGRAM shcfinal
