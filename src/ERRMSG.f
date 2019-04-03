c      ===================================================================================
c      print R-messages
c     ===================================================================================
c     THIS IS PART OF THE DOPRI INTEGRATOR

      subroutine dblepr_k(label, nchar, data, ndata)
      integer nchar, ndata
      character*(*) label
      double precision data(ndata)
      integer nc
       nc = nchar
       if(nc .lt. 0) nc = len(label)
c       call dblep0k(label, nc, data, ndata)
      end subroutine


      subroutine intpr_k(label, nchar, data, ndata)
      integer nchar, ndata
      character*(*) label
      integer data(ndata)
      integer nc
       nc = nchar
       if(nc .lt. 0) nc = len(label)
c       call intp0k(label, nc, data, ndata)
      end subroutine

C     just a string
      subroutine rprint(msg)
      character (len=*) msg
c           call dblepr(msg, -1, 0, 0)
      end subroutine


C     printing with one integer and a double
      subroutine rprintid(msg, i1, d1)
      character (len=*) msg
      double precision d1
      integer i1
c        call dblepr(msg, -1, d1, 1)
c        call intpr(' ', -1, i1, 1)
      end subroutine

      subroutine rprintd4(msg, d1, d2, d3, d4)
      character (len=*) msg
      double precision dbl(4), d1, d2, d3, d4
      integer i1
        DBL(1) = d1
        DBL(2) = d2
        DBL(3) = d3
        DBL(4) = d4

c     call dblepr(msg, -1, dbl, 4)
      end subroutine

      subroutine rprintd3(msg, d1, d2, d3)
      character (len=*) msg
      double precision dbl(3), d1, d2, d3
        DBL(1) = d1
        DBL(2) = d2
        DBL(3) = d3

c     call dblepr(msg, -1, dbl, 3)
      end subroutine

C     printing with one double
      subroutine rprintd1(msg, d1)
      character (len=*) msg
      double precision d1, DBL(1)
        DBL(1) = d1
        call dblepr_k(msg, -1, DBL, 1)
      end subroutine

C     printing with two doubles
      subroutine rprintd2(msg, d1, d2)
      character (len=*) msg
      double precision DBL(2), d1, d2
        DBL(1) = d1
        DBL(2) = d2
        call dblepr_k(msg, -1, DBL, 2)
      end subroutine

C     printing with one integer
      subroutine rprinti1(msg, i1)
      character (len=*) msg
      integer IN(1), i1
        IN(1) = i1
        call intpr_k(msg, -1, IN, 1)
      end subroutine

      subroutine rprinti2(msg, i1, i2)
      character (len=*) msg
      INTEGER IN(2), i1, i2
        IN(1) = i1
        IN(2) = i2
        call intpr_k(msg, -1, IN, 2)
      end subroutine

      subroutine rprinti3(msg, i1, i2, i3)
      character (len=*) msg
      INTEGER IN(3), i1, i2, i3
        IN(1) = i1
        IN(2) = i2
        IN(3) = i3
        call intpr_k(msg, -1, IN, 3)
      end subroutine
