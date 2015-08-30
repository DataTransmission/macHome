c<html>
c<head><title>pointers.f</title></head>
c<body>
c<pre>

      program pointer_test
      implicit none
c
c   Demonstrate use of pointers to manipulate sections of arrays
c   Notice that once associated, p2 and p1 are treated like normal
c   arrays in the program
c
c    John Mahaffy 4/17/96
c
c<a name="target"><font color="FF0000">
      real, target :: a(10,10)
c</font></a>
c<a name="point"><font color="FF0000">
      real, pointer :: p2(:,:), p1(:)
c</font></a>
      integer iub, i, j
c
c   Define a generic interface to permit output of information
c   for either rank 1 or rank 2 arrays.  You need such an interface
c   in any case to get some of the associated pointers to pass their
c   contents efficiently (quietly passes a skip factor (stride).
c   Without an interface our compiler allocates space and makes
c   a copy of the contents of the pointer array to guarantee that the
c   sequential elements in the passed array are in adjacent memory
c   locations
c
c<a name="int"><font color="FF0000">
      interface parray
c</font></a>
            subroutine p2array (a,title,aname)
            real a(:,:)
c<a name="char"><font color="FF0000">
            character*(*) title,aname
c</font></a>
            end subroutine
            subroutine p1array (a,title,aname)
            real a(:)
            character*(*) title,aname
            end subroutine
      end interface
c
c  Open a file to contain results printed to the screen
c
      open (9,file='pointers.out')
c
c    Set some initial values
c
      do    i=1,10
         do j=1,10
            a(i,j)=1000*i+j
         enddo
      enddo
c
      call parray(a, 'original array', 'a')
c
c   p2 becomes a rank 2 array associated with part of "a"
c
      p2 => a(3:6,7:9)
      call parray (p2,'p2 => a(3:6,7:9) ','p2')
c
c   add 800 to all elements in p2
c
      p2 = -p2 - 800
      call parray (p2,'results of p2 = -p2 - 800 ','p2')
c
c   What have we done to the array "a" ?
c
      call parray(a, 'Here''s what happened to "a" ', 'a')
c
c   Associate p1 with column 3 of "a"
c
      p1 => a(:,3)
      call parray (p1, 'Associate p1 with column 3: p1 => a(:,3)','p1')
c
c   Modify  p1
c     <a name=1><font color=FF0000>
      iub = ubound(p1,dim=1)
c     </font>
      do 40 i=1,iub
   40 p1(i) = i**2
      call parray (p1,'Modify p1 with the equation: p1(i) = i**2','p1')
c
c   Now associate p1 with row 8 of "a"
c
      p1 => a(8,:)
      call parray (p1, 'Associate p1 with row 8:  p1 => a(8,:)','p1')
c
c   Modify  p1
c
      iub=ubound(p1,dim=1)
      do 60 i=1,iub
   60 p1(i) = i
      call parray (p1, 'Modify p1 with the equation: p1(i) = i ','p1')
c
c   Now what does the array "a" look like ?
c
      call parray(a, '"a" after all changes to p1 and p2 ', 'a')
c
      stop
      end
c=======================================================================
      subroutine p2array (a,title,aname)
      implicit none
c
c    Print out a 2-D array with header information
c
c    John Mahaffy 4/17/96
c
c   INPUT arguments
c
c   a     -   array to be printed
c   title -   title information for the output
c   aname -   name of the array being printed
c
c
c   The use of colons in the following type declaration is
c   necessary so that the subroutine will look for some
c   hidden arguments passed as a result of an interface statement
c   in the calling routine.   These hidden arguments provide
c   detailed information on the shape of the array.
c
      real a(:,:)
c
      character*(*) title,aname
      character*(4) clb1,cub1,clb2,cub2
      character*32 info
      integer  lb1, ub1, lb2, ub2, i
c
c   Build a header giving the array name and shape
c     <a name=2><font color=FF0000>
      lb1 = lbound(a ,dim=1)
c     </font>
      ub1 = ubound(a ,dim=1)
      lb2 = lbound(a ,dim=2)
      ub2 = ubound (a ,dim=2)
      write(clb1,2001) lb1
      write(cub1,2001) ub1
      write(clb2,2001) lb2
      write(cub2,2001) ub2
c
      info= aname//'('//trim(adjustl(clb1))//':'
     &      //trim(adjustl(cub1))//','
     &      //trim(adjustl(clb2))//':'
     &      //trim(adjustl(cub2))//') ='
c
      write(*,2000) repeat('=',80)
      write(*,2000) title
      write(*,2000) repeat('-',80)
      write(*,2000) info
c
      write(9,2000) repeat('=',80)
      write(9,2000) title
      write(9,2000) repeat('-',80)
      write(9,2000) info
c
      do i = lb1,ub1
         write(*,2002) i, a(i,:)
         write( 9,2002) i, a(i,:)
      enddo
 2000 format (a)
 2001 format (i4)
 2002 format('row',i4,10f7.0,(7x,10f7.0))
      end
c=======================================================================
      subroutine p1array (a,title,aname)
      implicit none
c
c    Print out a 1-D array with header information
c
c
c    John Mahaffy 4/17/96
c
c   INPUT arguments
c
c   a     -   array to be printed
c   title -   title information for the output
c   aname -   name of the array being printed
c
c
c   The use of colon  in the following type declaration is
c   necessary so that the subroutine will look for some
c   hidden arguments passed as a result of an interface statement
c   in the calling routine. These hidden arguments provide
c   detailed information on the shape of the array
c
      real a(:)
c
      character*(*) title,aname
      character*(4) clb1,cub1
      character*32 info
      integer  lb1, ub1
c
c   Build a header giving the array name and shape
c
      lb1 = lbound(a ,dim=1)
      ub1 = ubound(a ,dim=1)
      write(clb1,2001) lb1
      write(cub1,2001) ub1
c
      info= aname//'('//trim(adjustl(clb1))//':'
     &      //trim(adjustl(cub1))//')  '
c
      write(*,2000) repeat('=',80)
      write(*,2000) title
      write(*,2000) repeat('-',80)
      write(*,2000) info
c
      write(9,2000) repeat('=',80)
      write(9,2000) title
      write(9,2000) repeat('-',80)
      write(9,2000) info
c
         write(*,2002)  a
         write( 9,2002)  a
 2000 format (a)
 2001 format (i4)
 2002 format(10f7.0)
      end
c</pre>
c</body>
c</html>
