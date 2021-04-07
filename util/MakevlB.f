      program MakevlB

*** write velocities into an unformatted file

      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)

      open(11, file = 'GRID4.2.5B.txt', status = 'old')
      open(12, file = 'Data4.2.5B.txt', status = 'unknown', 
     *         form = 'unformatted')

      read(11,100) idim, jdim
  100 format(30x,i3,20x,i3)
      write(6, 110) idim, jdim
  110 format(2i6)
      kdim = (idim + 1)*(jdim + 1)
      do 150 i = 1,kdim
         read(11, 120) f1, f2, f3, f4, f5, f6
  120    format(32x, 6f8.2)
         write(12) f1, f2, f3, f4, f5, f6
  150 continue      
      close(11)
      close(12)
      stop
      end

